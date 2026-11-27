"""
Titik masuk worker:

    python -m worker vcf    [--spec-file spec.json] [--validate-only]
    python -m worker report [--spec-file spec.json] [--validate-only]
    python -m worker trimming | alignment | variant_calling [--spec-file spec.json] [--validate-only]   (wizard per langkah)

Sumber spesifikasi (urutan): --spec-file, env JOB_SPEC (JSON), env JOB_SPEC_S3_URI (s3://...).
Exit code: 0 sukses, 1 gagal, 2 spesifikasi tidak valid, 143 dihentikan.
"""

import argparse
import json
import logging
import os
import sys
import traceback
from pathlib import Path
from typing import Optional

from pydantic import ValidationError

from . import config, proc
from .errors import (EXIT_FAILED, EXIT_INVALID_SPEC, EXIT_OK, INTERNAL_ERROR, INVALID_SPEC, WorkerError)
from .s3io import S3
from .spec import (STEP_DIRS, AlignmentStepSpec, JobSpec, ReportSpec, TrimmingStepSpec, VariantStepSpec,
                   normalize_prefix, validation_problems)
from .status import StatusReporter

log = logging.getLogger("worker")

SPEC_MODELS = {"vcf": JobSpec, "report": ReportSpec, "trimming": TrimmingStepSpec,
               "alignment": AlignmentStepSpec, "variant_calling": VariantStepSpec}
STEP_TASKS = ("trimming", "alignment", "variant_calling")


def _load_raw_spec(args, s3: Optional[S3]) -> dict:
    try:
        if args.spec_file:
            return json.loads(Path(args.spec_file).read_text())
        if os.environ.get("JOB_SPEC"):
            return json.loads(os.environ["JOB_SPEC"])
        if os.environ.get("JOB_SPEC_S3_URI"):
            return (s3 or S3()).read_json(os.environ["JOB_SPEC_S3_URI"])
    except (OSError, ValueError) as e:
        raise WorkerError(INVALID_SPEC, f"Spesifikasi tidak bisa dibaca: {e}", EXIT_INVALID_SPEC)
    raise WorkerError(INVALID_SPEC, "Tidak ada spesifikasi: isi JOB_SPEC, JOB_SPEC_S3_URI, atau --spec-file",
                      EXIT_INVALID_SPEC)


def _best_effort_prefix(raw: dict, task: str) -> Optional[str]:
    """Prefix status.json dari spesifikasi yang tidak valid, supaya pemanggil tetap bisa membaca kegagalannya."""
    try:
        if isinstance(raw.get("output"), dict) and raw["output"].get("prefix"):
            return normalize_prefix(raw["output"]["prefix"])
        step_dir = STEP_DIRS["qc" if task == "report" else task]
        return normalize_prefix(raw["analysisPrefix"]) + step_dir
    except Exception:  # noqa: BLE001
        return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="worker")
    parser.add_argument("task", choices=list(SPEC_MODELS))
    parser.add_argument("--spec-file")
    parser.add_argument("--validate-only", action="store_true", help="validasi spesifikasi lalu keluar (tanpa S3)")
    args = parser.parse_args(argv)

    logging.basicConfig(level=os.environ.get("LOG_LEVEL", "INFO"),
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    proc.install_signal_handlers()

    scratch = Path(config.SCRATCH_DIR)
    status: Optional[StatusReporter] = None
    s3: Optional[S3] = None
    try:
        if not args.validate_only:
            s3 = S3()
        raw = _load_raw_spec(args, s3)
        model = SPEC_MODELS[args.task]
        try:
            spec = model.model_validate(raw)
        except ValidationError as e:
            problems = validation_problems(e)
            # Tulis status FAILED bila prefix output bisa dibaca, agar pemanggil tahu tanpa membuka log
            prefix = None if args.validate_only else _best_effort_prefix(raw, args.task)
            if prefix and s3:
                st = StatusReporter(args.task, str(raw.get("jobId", "unknown")), prefix, s3, scratch / "status")
                st.fail(INVALID_SPEC, "Spesifikasi job tidak valid", {"fields": problems})
            log.error("Spesifikasi tidak valid: %s", json.dumps(problems, ensure_ascii=False))
            return EXIT_INVALID_SPEC

        if args.validate_only:
            normalized = spec.model_dump(mode="json")
            if args.task in ("vcf",) + STEP_TASKS:
                normalized["_derived"] = {"sampleName": spec.sample, "statusPrefix": spec.status_prefix}
            print(json.dumps(normalized, indent=2))
            return EXIT_OK

        scratch.mkdir(parents=True, exist_ok=True)
        if not os.access(scratch, os.W_OK):
            raise WorkerError(INTERNAL_ERROR, f"{scratch} tidak bisa ditulis oleh user container "
                                              f"(uid {os.getuid()}). Jalankan task sebagai root atau perbaiki izin volume.")
        status = StatusReporter(args.task, spec.jobId, spec.status_prefix, s3, scratch / "status" / spec.jobId)
        if args.task in STEP_TASKS:
            from . import step_jobs
            step_jobs.prepare(args.task, spec, s3)   # cek input + bersihkan hasil lama/hilir (kegagalan tercatat di status.json)
        status.start()

        if args.task == "vcf":
            from . import vcf_job
            vcf_job.run(spec, status, s3)
        elif args.task == "report":
            from . import report_job
            report_job.run(spec, status, s3)
        else:
            from . import step_jobs
            {"trimming": step_jobs.run_trimming, "alignment": step_jobs.run_alignment,
             "variant_calling": step_jobs.run_variant_calling}[args.task](spec, status, s3)
        return EXIT_OK

    except WorkerError as e:
        if status:
            status.fail(e.code, e.message, e.details)
        else:
            log.error("%s: %s", e.code, e.message)
        return e.exit_code
    except Exception as e:  # noqa: BLE001
        log.error("Error tak terduga:\n%s", traceback.format_exc())
        if status:
            status.fail(INTERNAL_ERROR, f"{type(e).__name__}: {e}")
        return EXIT_FAILED
