"""
NGS API — Format error sesuai API Contract bagian 4:

    { "error": { "code": "STATUS_CONFLICT", "message": "...", "details": {...} } }

Semua error (HTTPException, validasi request, error tak terduga) diseragamkan ke format ini.
"""

import logging

from fastapi import FastAPI, HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from starlette.exceptions import HTTPException as StarletteHTTPException

logger = logging.getLogger("ngs.api")

_CODE_BY_STATUS = {
    400: "VALIDATION_FAILED",
    401: "UNAUTHENTICATED",
    403: "FORBIDDEN",
    404: "NOT_FOUND",
    405: "VALIDATION_FAILED",
    409: "STATUS_CONFLICT",
    422: "VALIDATION_FAILED",
}


def error_response(code: str, message: str, http_status: int, details: dict = None):
    """Lempar HTTPException dengan envelope kontrak."""
    body = {"error": {"code": code, "message": message}}
    if details:
        body["error"]["details"] = details
    raise HTTPException(status_code=http_status, detail=body)


def _envelope(code: str, message: str, details: dict = None) -> dict:
    body = {"error": {"code": code, "message": message}}
    if details:
        body["error"]["details"] = details
    return body


def install_handlers(app: FastAPI):

    @app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(request: Request, exc: StarletteHTTPException):
        detail = exc.detail
        if isinstance(detail, dict) and "error" in detail:
            body = detail
        else:
            code = _CODE_BY_STATUS.get(exc.status_code, "SERVICE_UNAVAILABLE" if exc.status_code >= 500 else "VALIDATION_FAILED")
            body = _envelope(code, str(detail))
        return JSONResponse(status_code=exc.status_code, content=body, headers=getattr(exc, "headers", None))

    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(request: Request, exc: RequestValidationError):
        fields = [
            {"path": ".".join(str(p) for p in err.get("loc", ())), "reason": str(err.get("msg", ""))}
            for err in exc.errors()
        ]
        return JSONResponse(
            status_code=400,
            content=_envelope("VALIDATION_FAILED", "Request tidak valid.", {"fields": fields}),
        )

    @app.exception_handler(Exception)
    async def unhandled_exception_handler(request: Request, exc: Exception):
        # Detail internal hanya di log server; klien tidak menerima path atau stack trace
        logger.exception("Error tak terduga di %s %s", request.method, request.url.path)
        return JSONResponse(
            status_code=503,
            content=_envelope("SERVICE_UNAVAILABLE", "Terjadi kesalahan di server."),
        )
