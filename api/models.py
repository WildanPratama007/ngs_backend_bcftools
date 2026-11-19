"""
NGS API — Pydantic models untuk request/response
Fix: enum status sesuai API Contract, format waktu ISO 8601 UTC,
     field opsional null bukan hilang
"""

from pydantic import BaseModel, Field
from typing import Optional
from datetime import datetime, timezone
from enum import Enum


class AnalysisMode(str, Enum):
    WGS      = "WGS"
    WES      = "WES"
    TARGETED = "Targeted"


class EndType(str, Enum):
    SE = "SE"
    PE = "PE"


class ReferenceGenome(str, Enum):
    GRCh38_chr22 = "GRCh38_chr22"
    GRCh38       = "GRCh38"
    GRCh37       = "GRCh37"


class VariantType(str, Enum):
    GERMLINE = "germline"
    SOMATIC  = "somatic"


class TrimPreset(str, Enum):
    """Preset kualitas trimming — ditampilkan di UI"""
    RESEARCH = "research"   # SLIDINGWINDOW Q20, MINLEN 36
    CLINICAL = "clinical"   # SLIDINGWINDOW Q30, MINLEN 36


# JobStatus sesuai API Contract section 2.3
class JobStatus(str, Enum):
    QUEUED            = "QUEUED"
    RUNNING           = "RUNNING"
    CANCEL_REQUESTED  = "CANCEL_REQUESTED"
    SUCCEEDED         = "SUCCEEDED"
    FAILED            = "FAILED"
    CANCELLED         = "CANCELLED"


# Alias untuk backward compat internal
class StageStatus(str, Enum):
    PENDING   = "pending"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    SKIPPED   = "skipped"


# ── Response Models ───────────────────────────────────────────────────────────

class StageInfo(BaseModel):
    status:      StageStatus    = StageStatus.PENDING
    progress:    int            = 0        # 0-100
    started_at:  Optional[datetime] = None
    finished_at: Optional[datetime] = None
    error:       Optional[str]  = None

    model_config = {"json_encoders": {datetime: lambda v: v.isoformat() + "Z" if v else None}}


class CoverageStats(BaseModel):
    """Extended coverage statistics — lebih informatif dari mean saja."""
    mean_depth:     float = 0.0
    median_depth:   float = 0.0
    pct_bases_10x:  float = 0.0   # % bases dengan depth ≥10x
    pct_bases_20x:  float = 0.0   # % bases dengan depth ≥20x
    pct_bases_30x:  float = 0.0   # % bases dengan depth ≥30x
    total_bases:    int   = 0
    covered_bases:  int   = 0
    warning:        Optional[str] = None


class PipelineStatus(BaseModel):
    job_id:        str
    sample_name:   str
    status:        StageStatus
    analysis_mode: AnalysisMode
    end_type:      EndType
    reference:     ReferenceGenome
    variant_type:  VariantType
    coverage:      Optional[CoverageStats] = None   # extended coverage stats
    stages: dict = Field(default_factory=lambda: {
        "fastqc":          StageInfo().model_dump(),
        "trimming":        StageInfo().model_dump(),
        "alignment":       StageInfo().model_dump(),
        "variant_calling": StageInfo().model_dump(),
    })
    created_at:   datetime = Field(default_factory=datetime.now)
    elapsed:      Optional[str] = None
    nextflow_log: Optional[str] = None


class SubmitResponse(BaseModel):
    job_id:   str
    message:  str
    status:   str


class ErrorResponse(BaseModel):
    error: dict  # {"code": "...", "message": "...", "details": {...}}
