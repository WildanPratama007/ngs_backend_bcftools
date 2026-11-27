"""Akses S3: unduh input, unggah hasil, tulis JSON. Kredensial dari IAM task role (rantai kredensial boto3)."""

import json
import logging
import mimetypes
import os
from pathlib import Path
from typing import Iterable, List, Tuple

import boto3
from boto3.s3.transfer import TransferConfig
from botocore.config import Config
from botocore.exceptions import BotoCoreError, ClientError

from . import config
from .errors import DOWNLOAD_FAILED, INPUT_NOT_FOUND, UPLOAD_FAILED, WorkerError
from .spec import parse_s3_uri

log = logging.getLogger("worker.s3")

# Multipart besar dan paralel: FASTQ/BAM berukuran puluhan GB
_TRANSFER = TransferConfig(
    multipart_threshold=64 << 20,
    multipart_chunksize=64 << 20,
    max_concurrency=10,
    use_threads=True,
)


def _client():
    kwargs = {"region_name": config.AWS_REGION,
              "config": Config(retries={"max_attempts": 8, "mode": "standard"},
                               s3={"addressing_style": "path"} if config.S3_ENDPOINT_URL else {})}
    if config.S3_ENDPOINT_URL:
        kwargs["endpoint_url"] = config.S3_ENDPOINT_URL
    return boto3.client("s3", **kwargs)


class S3:
    def __init__(self):
        self.client = _client()

    # ── baca ─────────────────────────────────────────────────────────────────
    def size_of(self, uri: str) -> int:
        bucket, key = parse_s3_uri(uri)
        try:
            return int(self.client.head_object(Bucket=bucket, Key=key)["ContentLength"])
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code", "")
            if code in ("404", "NoSuchKey", "NotFound"):
                raise WorkerError(INPUT_NOT_FOUND, f"File input tidak ditemukan: {uri}", details={"uri": uri})
            if code in ("403", "AccessDenied"):
                raise WorkerError(INPUT_NOT_FOUND, f"Tidak punya akses membaca {uri} (cek IAM task role)", details={"uri": uri})
            raise WorkerError(DOWNLOAD_FAILED, f"Gagal membaca metadata {uri}: {code}", details={"uri": uri})
        except BotoCoreError as e:
            raise WorkerError(DOWNLOAD_FAILED, f"Gagal menghubungi S3 untuk {uri}: {e}", details={"uri": uri})

    def download(self, uri: str, dest: Path) -> int:
        bucket, key = parse_s3_uri(uri)
        dest.parent.mkdir(parents=True, exist_ok=True)
        expected = self.size_of(uri)
        try:
            self.client.download_file(bucket, key, str(dest), Config=_TRANSFER)
        except (ClientError, BotoCoreError, OSError) as e:
            raise WorkerError(DOWNLOAD_FAILED, f"Gagal mengunduh {uri}: {e}", details={"uri": uri})
        got = dest.stat().st_size
        if got != expected:
            raise WorkerError(DOWNLOAD_FAILED, f"Ukuran {uri} tidak cocok: {got} != {expected}", details={"uri": uri})
        return got

    def read_json(self, uri: str) -> dict:
        bucket, key = parse_s3_uri(uri)
        try:
            body = self.client.get_object(Bucket=bucket, Key=key)["Body"].read()
        except (ClientError, BotoCoreError) as e:
            raise WorkerError(INPUT_NOT_FOUND, f"Gagal membaca {uri}: {e}", details={"uri": uri})
        return json.loads(body)

    def list_keys(self, uri_prefix: str) -> List[Tuple[str, int]]:
        bucket, key = parse_s3_uri(uri_prefix)
        out = []
        for page in self.client.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=key):
            for obj in page.get("Contents", []):
                out.append((f"s3://{bucket}/{obj['Key']}", int(obj["Size"])))
        return out

    # ── tulis ────────────────────────────────────────────────────────────────
    def _extra(self, path_or_key: str) -> dict:
        extra = {}
        ctype = mimetypes.guess_type(path_or_key)[0]
        if ctype:
            extra["ContentType"] = ctype
        if config.S3_SSE:
            extra["ServerSideEncryption"] = config.S3_SSE
            if config.S3_SSE == "aws:kms" and config.S3_KMS_KEY_ID:
                extra["SSEKMSKeyId"] = config.S3_KMS_KEY_ID
        return extra

    def upload(self, src: Path, uri: str) -> int:
        bucket, key = parse_s3_uri(uri)
        size = src.stat().st_size
        try:
            self.client.upload_file(str(src), bucket, key, ExtraArgs=self._extra(key), Config=_TRANSFER)
            remote = int(self.client.head_object(Bucket=bucket, Key=key)["ContentLength"])
        except (ClientError, BotoCoreError, OSError) as e:
            raise WorkerError(UPLOAD_FAILED, f"Gagal mengunggah ke {uri}: {e}", details={"uri": uri})
        if remote != size:
            raise WorkerError(UPLOAD_FAILED, f"Ukuran {uri} setelah upload tidak cocok: {remote} != {size}", details={"uri": uri})
        return size

    def delete_prefix(self, uri_prefix: str) -> int:
        """Hapus semua objek di bawah prefix (harus berakhir dengan '/'). Mengembalikan jumlah objek yang dihapus."""
        if not uri_prefix.endswith("/"):
            raise WorkerError(UPLOAD_FAILED, f"delete_prefix menolak prefix tanpa '/': {uri_prefix}")
        bucket, _ = parse_s3_uri(uri_prefix)
        keys = [u.split(f"s3://{bucket}/", 1)[1] for u, _ in self.list_keys(uri_prefix)]
        deleted = 0
        try:
            for i in range(0, len(keys), 1000):
                batch = [{"Key": k} for k in keys[i:i + 1000]]
                self.client.delete_objects(Bucket=bucket, Delete={"Objects": batch, "Quiet": True})
                deleted += len(batch)
        except (ClientError, BotoCoreError) as e:
            raise WorkerError(UPLOAD_FAILED, f"Gagal menghapus hasil lama di {uri_prefix}: {e}", details={"uri": uri_prefix})
        return deleted

    def exists(self, uri: str) -> bool:
        bucket, key = parse_s3_uri(uri)
        try:
            self.client.head_object(Bucket=bucket, Key=key)
            return True
        except ClientError:
            return False

    def put_json(self, uri: str, payload: dict):
        bucket, key = parse_s3_uri(uri)
        extra = {k: v for k, v in self._extra(key).items() if k != "ContentType"}
        try:
            self.client.put_object(Bucket=bucket, Key=key, Body=json.dumps(payload, indent=2).encode("utf-8"),
                                   ContentType="application/json", **extra)
        except (ClientError, BotoCoreError) as e:
            raise WorkerError(UPLOAD_FAILED, f"Gagal menulis {uri}: {e}", details={"uri": uri})


def files_under(root: Path) -> Iterable[Path]:
    for p in sorted(root.rglob("*")):
        if p.is_file():
            yield p
