"""
Spesifikasi job yang dikirim ke container ECS (JSON). Field tak dikenal ditolak.
Aturan rentang sama dengan API (/api/pipeline/submit). Lihat docs/ECS_JOB_SPEC.md.
"""

import os
import re
from enum import Enum
from typing import ClassVar, Optional, List

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

_S3_URI_RE = re.compile(r"^s3://([a-z0-9][a-z0-9.\-]{1,61}[a-z0-9])/(.+)$")
_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")   # gzip maupun teks biasa


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Layout(str, Enum):
    PAIRED = "PAIRED"
    SINGLE = "SINGLE"


class Reference(str, Enum):
    GRCh38 = "GRCh38"
    GRCh37 = "GRCh37"
    GRCh38_chr22 = "GRCh38_chr22"


class AnalysisMode(str, Enum):
    WGS = "WGS"
    WES = "WES"
    Targeted = "Targeted"


class AdapterKit(str, Enum):
    truseq_pe = "truseq_pe"
    truseq_se = "truseq_se"
    truseq_pcr_free = "truseq_pcr_free"
    nextera = "nextera"
    nextera_dna_flex = "nextera_dna_flex"


class Preset(str, Enum):
    research = "research"
    clinical = "clinical"


class Sex(str, Enum):
    male = "male"
    female = "female"


def parse_s3_uri(uri: str):
    m = _S3_URI_RE.match(uri or "")
    if not m:
        raise ValueError(f"bukan URI S3 yang valid (s3://bucket/key): {uri!r}")
    return m.group(1), m.group(2)


def normalize_prefix(v: str) -> str:
    bucket, key = parse_s3_uri(v)   # key tidak boleh kosong: hasil tidak ditulis di akar bucket
    return f"s3://{bucket}/{key.rstrip('/')}/"


def check_job_id(v: str) -> str:
    if not _ID_RE.match(v) or ".." in v:
        raise ValueError("jobId hanya huruf/angka/_/-/. (maks 64, diawali huruf atau angka)")
    return v


def sanitize_name(name: str) -> str:
    """Hanya A-Z a-z 0-9 . - _ (sama dengan Sample_Name_Sanitizer di API)."""
    s = re.sub(r"[^A-Za-z0-9._\-]", "_", name)
    s = re.sub(r"_+", "_", s).strip("_").lstrip(".")
    return s or "sample"


def derive_sample_name(r1_key: str, layout: "Layout") -> str:
    name = os.path.basename(r1_key)
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.IGNORECASE)
    if layout == Layout.PAIRED:
        name = re.sub(r"[_.]R?1(?:_\d+)?$", "", name)   # _R1, _R1_001, .R1, _1, _1_001
    return sanitize_name(name)


class InputSpec(_Strict):
    layout: Layout
    r1: str
    r2: Optional[str] = None

    @field_validator("r1", "r2")
    @classmethod
    def _fastq_uri(cls, v):
        if v is None:
            return v
        _, key = parse_s3_uri(v)
        if not key.endswith(_FASTQ_SUFFIXES):
            raise ValueError("file harus berakhiran .fastq.gz, .fq.gz, .fastq, atau .fq (huruf kecil)")
        return v

    @model_validator(mode="after")
    def _layout_matches_files(self):
        if self.layout == Layout.PAIRED and not self.r2:
            raise ValueError("layout PAIRED membutuhkan r2")
        if self.layout == Layout.SINGLE and self.r2:
            raise ValueError("layout SINGLE tidak memakai r2")
        if self.r2 and self.r1 == self.r2:
            raise ValueError("r1 dan r2 tidak boleh file yang sama")
        return self


class OutputSpec(_Strict):
    prefix: str
    includeBam: bool = True      # upload BAM (post-markdup) + .bai

    @field_validator("prefix")
    @classmethod
    def _prefix(cls, v):
        return normalize_prefix(v)


class TrimmingSpec(_Strict):
    preset: Optional[Preset] = None   # menimpa quality/minLen: research=Q20, clinical=Q30 (minLen 36)
    quality: int = Field(20, ge=10, le=40)
    minLen: int = Field(36, ge=20, le=100)
    adapter: bool = True
    adapterKit: Optional[AdapterKit] = None   # kosong = TruSeq SE untuk SINGLE, TruSeq PE untuk PAIRED
    polyG: bool = False

    def effective(self):
        if self.preset == Preset.clinical:
            return 30, 36
        if self.preset == Preset.research:
            return 20, 36
        return self.quality, self.minLen


class AlignmentSpec(_Strict):
    markDuplicates: bool = True


class VariantCallingSpec(_Strict):
    minQual: int = Field(20, ge=0, le=999)
    minDp: int = Field(10, ge=1, le=500)
    minGq: int = Field(20, ge=0, le=99)
    maxDepth: int = Field(500, ge=100, le=10000)
    minBq: int = Field(20, ge=0, le=40)
    minMq: int = Field(20, ge=0, le=60)


class ResourcesSpec(_Strict):
    threads: Optional[int] = Field(None, ge=1, le=64)   # kosong = otomatis dari CPU container (maks 16)


class JobSpec(_Strict):
    """Spesifikasi untuk `worker vcf`."""
    jobId: str
    sampleName: Optional[str] = None     # kosong = diturunkan dari nama file R1
    input: InputSpec
    output: OutputSpec
    reference: Reference
    analysisMode: AnalysisMode = AnalysisMode.WGS
    bed: Optional[str] = None            # wajib untuk WES/Targeted
    sampleSex: Optional[Sex] = None
    trimming: TrimmingSpec = TrimmingSpec()
    alignment: AlignmentSpec = AlignmentSpec()
    variantCalling: VariantCallingSpec = VariantCallingSpec()
    resources: ResourcesSpec = ResourcesSpec()

    @field_validator("jobId")
    @classmethod
    def _job_id(cls, v):
        return check_job_id(v)

    @field_validator("sampleName")
    @classmethod
    def _sample(cls, v):
        if v is not None and not re.match(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$", v):
            raise ValueError("sampleName hanya huruf/angka/_/-/. (maks 64, diawali huruf atau angka)")
        return v

    @field_validator("bed")
    @classmethod
    def _bed(cls, v):
        if v is None:
            return v
        _, key = parse_s3_uri(v)
        if not key.endswith(".bed"):
            raise ValueError("file harus berakhiran .bed")
        return v

    @model_validator(mode="after")
    def _cross(self):
        needs_bed = self.analysisMode in (AnalysisMode.WES, AnalysisMode.Targeted)
        if needs_bed and not self.bed:
            raise ValueError(f"analysisMode {self.analysisMode.value} membutuhkan bed")
        if not needs_bed and self.bed:
            raise ValueError("bed hanya untuk analysisMode WES atau Targeted")
        return self

    @property
    def sample(self) -> str:
        return self.sampleName or derive_sample_name(parse_s3_uri(self.input.r1)[1], self.input.layout)

    @property
    def status_prefix(self) -> str:
        return self.output.prefix


class ReportOutputSpec(_Strict):
    prefix: str

    @field_validator("prefix")
    @classmethod
    def _prefix(cls, v):
        return normalize_prefix(v)


class ReportSpec(_Strict):
    """
    Spesifikasi untuk `worker report` (FastQC, Step 1 wizard).
    Isi `output.prefix` (lokasi bebas) ATAU `analysisPrefix` (hasil ke `<analysisPrefix>qc/`, konsisten dengan step lain).
    """
    jobId: str
    input: InputSpec
    output: Optional[ReportOutputSpec] = None
    analysisPrefix: Optional[str] = None
    resources: ResourcesSpec = ResourcesSpec()

    @field_validator("jobId")
    @classmethod
    def _job_id(cls, v):
        return check_job_id(v)

    @field_validator("analysisPrefix")
    @classmethod
    def _analysis_prefix(cls, v):
        return normalize_prefix(v) if v is not None else v

    @model_validator(mode="after")
    def _one_prefix(self):
        if bool(self.output) == bool(self.analysisPrefix):
            raise ValueError("isi salah satu: output.prefix atau analysisPrefix")
        return self

    @property
    def status_prefix(self) -> str:
        return self.output.prefix if self.output else self.analysisPrefix + "qc/"

    @property
    def report_dir(self) -> str:
        """Subfolder laporan FastQC: `qc/` di bawah output.prefix; langsung di folder qc/ bila memakai analysisPrefix."""
        return "qc/" if self.output else ""


# ── Worker per langkah (wizard): trimming, alignment, variant_calling ─────────
# Susunan S3 per analisis (analysisPrefix = s3://bucket/analyses/<id>/):
#   qc/ (report) | trimmed/ (step 2) | alignment/ (step 3) | variants/ (step 4); tiap folder berisi status.json
STEP_DIRS = {"qc": "qc/", "trimming": "trimmed/", "alignment": "alignment/", "variant_calling": "variants/"}
_DOWNSTREAM = {"trimming": ["alignment", "variant_calling"], "alignment": ["variant_calling"], "variant_calling": []}


class _StepBase(_Strict):
    jobId: str
    analysisPrefix: str
    sampleName: Optional[str] = None
    resources: ResourcesSpec = ResourcesSpec()

    STEP: ClassVar[str] = ""

    @field_validator("jobId")
    @classmethod
    def _job_id(cls, v):
        return check_job_id(v)

    @field_validator("analysisPrefix")
    @classmethod
    def _analysis_prefix(cls, v):
        return normalize_prefix(v)

    @field_validator("sampleName")
    @classmethod
    def _sample(cls, v):
        if v is not None and not re.match(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$", v):
            raise ValueError("sampleName hanya huruf/angka/_/-/. (maks 64, diawali huruf atau angka)")
        return v

    @property
    def status_prefix(self) -> str:
        return self.analysisPrefix + STEP_DIRS[self.STEP]

    def step_prefix(self, step: str) -> str:
        return self.analysisPrefix + STEP_DIRS[step]

    @property
    def downstream_prefixes(self) -> List[str]:
        return [self.step_prefix(s) for s in _DOWNSTREAM[self.STEP]]


def _check_bed(mode: "AnalysisMode", bed: Optional[str]):
    needs_bed = mode in (AnalysisMode.WES, AnalysisMode.Targeted)
    if needs_bed and not bed:
        raise ValueError(f"analysisMode {mode.value} membutuhkan bed")
    if not needs_bed and bed:
        raise ValueError("bed hanya untuk analysisMode WES atau Targeted")
    if bed:
        _, key = parse_s3_uri(bed)
        if not key.endswith(".bed"):
            raise ValueError("bed harus berakhiran .bed")


class TrimmingStepSpec(_StepBase):
    """Step 2: FASTQ mentah (S3) -> <analysisPrefix>trimmed/ (FASTQ ter-gzip + stats.json)."""
    STEP: ClassVar[str] = "trimming"
    input: InputSpec
    trimming: TrimmingSpec = TrimmingSpec()

    @property
    def sample(self) -> str:
        return self.sampleName or derive_sample_name(parse_s3_uri(self.input.r1)[1], self.input.layout)

    @property
    def layout(self) -> Layout:
        return self.input.layout


class AlignmentStepSpec(_StepBase):
    """Step 3: <analysisPrefix>trimmed/ -> <analysisPrefix>alignment/ (BAM markdup + coverage)."""
    STEP: ClassVar[str] = "alignment"
    layout: Layout
    reference: Reference
    analysisMode: AnalysisMode = AnalysisMode.WGS
    bed: Optional[str] = None
    alignment: AlignmentSpec = AlignmentSpec()

    @model_validator(mode="after")
    def _cross(self):
        if not self.sampleName:
            raise ValueError("sampleName wajib (pakai nilai yang sama dengan step trimming)")
        _check_bed(self.analysisMode, self.bed)
        return self

    @property
    def sample(self) -> str:
        return self.sampleName


class VariantStepSpec(_StepBase):
    """Step 4: <analysisPrefix>alignment/ -> <analysisPrefix>variants/ (VCF terfilter)."""
    STEP: ClassVar[str] = "variant_calling"
    reference: Reference
    analysisMode: AnalysisMode = AnalysisMode.WGS
    bed: Optional[str] = None
    sampleSex: Optional[Sex] = None
    variantCalling: VariantCallingSpec = VariantCallingSpec()

    @model_validator(mode="after")
    def _cross(self):
        if not self.sampleName:
            raise ValueError("sampleName wajib (pakai nilai yang sama dengan step sebelumnya)")
        _check_bed(self.analysisMode, self.bed)
        return self

    @property
    def sample(self) -> str:
        return self.sampleName


def validation_problems(exc: ValidationError) -> List[dict]:
    """Ubah error pydantic menjadi [{path, reason}] untuk status.json."""
    out = []
    for e in exc.errors():
        out.append({"path": ".".join(str(p) for p in e.get("loc", ())), "reason": str(e.get("msg", "")).replace("Value error, ", "")})
    return out
