"""
NGS API — Pydantic models untuk request/response
"""

from pydantic import BaseModel, Field
from typing import Optional, Literal
from datetime import datetime
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


class StageStatus(str, Enum):
    PENDING   = "pending"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    SKIPPED   = "skipped"


# ── Response Models ───────────────────────────────────────────────────────────

class StageInfo(BaseModel):
    status:   StageStatus = StageStatus.PENDING
    progress: int         = 0       # 0-100
    started_at:  Optional[datetime] = None
    finished_at: Optional[datetime] = None
    error:    Optional[str]         = None


class PipelineStatus(BaseModel):
    job_id:        str
    sample_name:   str
    status:        StageStatus
    analysis_mode: AnalysisMode
    end_type:      EndType
    reference:     ReferenceGenome
    variant_type:  VariantType
    coverage:      Optional[float] = None   # mean coverage setelah alignment
    stages: dict = Field(default_factory=lambda: {
        "fastqc":          StageInfo().model_dump(),
        "trimming":        StageInfo().model_dump(),
        "alignment":       StageInfo().model_dump(),
        "variant_calling": StageInfo().model_dump(),
    })
    created_at:  datetime = Field(default_factory=datetime.now)
    elapsed:     Optional[str] = None
    nextflow_log: Optional[str] = None


class SubmitResponse(BaseModel):
    job_id:     str
    message:    str
    status:     str


class ErrorResponse(BaseModel):
    detail: str
