"""Self-contained Matplotlib export helper for the supplementary package.

This is a trimmed copy of the subset of the ``scientific-visualization``
``figure_export`` module that ``make_figures.py`` relies on.  It depends only on
Matplotlib (and Pillow for the optional RGB-TIFF normalisation) -- no external
asset files, no private paths.

Behaviour preserved from the original tool:

* files are written atomically through a same-directory temporary file;
* ``bbox_inches=None`` keeps the exact physical page dimensions of the figure;
* TIFF output uses LZW compression and is flattened to RGB before saving;
* PDF text stays editable (TrueType, ``pdf.fonttype = 42``);
* a JSON manifest (``<name>.export.json``) records settings, versions and the
  provenance dictionary when ``write_manifest=True``.

One deliberate difference: ``overwrite`` defaults to ``True`` so that the
pipeline can be re-run in place, which is what anyone reproducing the
figures actually needs.
"""

from __future__ import annotations

import importlib.metadata
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

SCHEMA_VERSION = "1.0"

SUPPORTED_FORMATS = {"pdf", "svg", "eps", "ps", "png", "tiff", "jpeg", "webp"}
FORMAT_ALIASES = {"jpg": "jpeg", "tif": "tiff"}
VECTOR_FORMATS = {"pdf", "svg", "eps", "ps"}
METADATA_FORMATS = {"pdf", "svg", "png", "eps", "ps"}


class ExportError(ValueError):
    """An expected export validation error."""


def _normalize_formats(formats: Iterable[str]) -> list[str]:
    normalized: list[str] = []
    for item in formats:
        for value in item.split(","):
            raw_name = value.strip().lower().lstrip(".")
            name = FORMAT_ALIASES.get(raw_name, raw_name)
            if not name:
                continue
            if name not in SUPPORTED_FORMATS:
                raise ExportError(
                    f"unsupported format {name!r}; "
                    f"available: {', '.join(sorted(SUPPORTED_FORMATS))}"
                )
            if name not in normalized:
                normalized.append(name)
    if not normalized:
        raise ExportError("at least one output format is required")
    return normalized


def _base_output_path(value: str | os.PathLike[str]) -> Path:
    path = Path(value).expanduser()
    if path.suffix.lower().lstrip(".") in SUPPORTED_FORMATS | set(FORMAT_ALIASES):
        path = path.with_suffix("")
    if not path.name:
        raise ExportError("output must include a base filename")
    return path


def _font_rc(font_mode: str) -> dict[str, Any]:
    if font_mode == "current":
        return {}
    if font_mode == "truetype":
        return {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}
    if font_mode == "paths":
        return {"pdf.fonttype": 3, "ps.fonttype": 3, "svg.fonttype": "path"}
    raise ExportError("font_mode must be 'truetype', 'paths', or 'current'")


def _package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _atomic_savefig(
    fig: Any,
    destination: Path,
    *,
    format_name: str,
    save_kwargs: dict[str, Any],
    tiff_rgb: bool = False,
) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.stem}.",
        suffix=f".{format_name}",
        dir=destination.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        fig.savefig(temporary, format=format_name, **save_kwargs)
        if tiff_rgb and format_name == "tiff":
            from PIL import Image

            with Image.open(temporary) as rendered:
                if "A" in rendered.getbands() and rendered.getchannel(
                    "A"
                ).getextrema() != (255, 255):
                    raise ExportError("RGB TIFF requires an opaque background")
                rgb = rendered.convert("RGB")
                pillow_options = dict(save_kwargs.get("pil_kwargs", {}))
                pillow_options.setdefault(
                    "dpi", rendered.info.get("dpi", (save_kwargs["dpi"],) * 2)
                )
                if rendered.info.get("icc_profile"):
                    pillow_options.setdefault("icc_profile", rendered.info["icc_profile"])
            try:
                rgb.save(temporary, format="TIFF", **pillow_options)
            finally:
                rgb.close()
        if not temporary.exists() or temporary.stat().st_size == 0:
            raise ExportError(f"Matplotlib produced no data for {destination}")
        os.replace(temporary, destination)
    except ExportError:
        raise
    except Exception as exc:  # pragma: no cover - depends on matplotlib internals
        raise ExportError(
            f"failed to export {destination.name} as {format_name}: {exc}"
        ) from exc
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _validate_provenance(provenance: dict[str, Any] | None) -> dict[str, Any]:
    if provenance is None:
        return {}
    if not isinstance(provenance, dict):
        raise ExportError("provenance must be a dictionary")
    try:
        encoded = json.dumps(provenance, sort_keys=True).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ExportError(f"provenance must be JSON serializable: {exc}") from exc
    if len(encoded) > 1_000_000:
        raise ExportError("provenance metadata is limited to 1,000,000 bytes")
    return provenance


def export_figure(
    fig: Any,
    filename: str | os.PathLike[str],
    *,
    formats: Iterable[str] = ("pdf", "png"),
    dpi: float = 300,
    transparent: bool = False,
    tiff_rgb: bool = False,
    bbox_inches: str | None = None,
    pad_inches: float = 0.1,
    facecolor: str = "white",
    edgecolor: str = "none",
    font_mode: str = "truetype",
    overwrite: bool = True,
    mkdir: bool = False,
    metadata: dict[str, Any] | None = None,
    provenance: dict[str, Any] | None = None,
    write_manifest: bool = False,
    savefig_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Export a Matplotlib figure atomically with explicit format settings.

    ``dpi`` is passed to vector backends too because it controls rasterized
    artists and embedded raster images.  ``bbox_inches=None`` preserves the
    figure's physical page dimensions; use ``"tight"`` only when that change is
    intentional.
    """
    if not (float(dpi) > 0 and math.isfinite(float(dpi))):
        raise ExportError("dpi must be finite and greater than zero")
    if not (float(pad_inches) >= 0 and math.isfinite(float(pad_inches))):
        raise ExportError("pad_inches must be finite and non-negative")
    normalized_formats = _normalize_formats(formats)
    if tiff_rgb and ("tiff" not in normalized_formats or transparent):
        raise ExportError("tiff_rgb requires TIFF output and transparent=False")

    base = _base_output_path(filename)
    parent = base.parent
    if mkdir:
        parent.mkdir(parents=True, exist_ok=True)
    if not parent.exists() or not parent.is_dir():
        raise ExportError(f"output parent directory does not exist: {parent}")

    import matplotlib as mpl

    provenance_document = _validate_provenance(provenance)
    user_metadata = dict(metadata or {})
    creator = user_metadata.setdefault(
        "Creator",
        f"Matplotlib {mpl.__version__}",
    )
    if not isinstance(creator, str):
        raise ExportError("metadata Creator must be a string")

    extra_kwargs = dict(savefig_kwargs or {})
    if "format" in extra_kwargs:
        raise ExportError("savefig_kwargs must not override format")

    figure_size = [float(value) for value in fig.get_size_inches()]
    outputs: list[dict[str, Any]] = []
    warnings_list: list[str] = []
    export_rc = _font_rc(font_mode)
    if bbox_inches is None:
        # Prevent a user's global savefig.bbox="tight" from changing page size.
        export_rc["savefig.bbox"] = None

    with mpl.rc_context(export_rc):
        for format_name in normalized_formats:
            destination = parent / f"{base.name}.{format_name}"
            if format_name in {"eps", "ps"}:
                format_metadata = {"Creator": user_metadata["Creator"]}
            elif format_name in METADATA_FORMATS:
                format_metadata = user_metadata
            else:
                format_metadata = None
            if user_metadata and format_name not in METADATA_FORMATS:
                warnings_list.append(
                    f"{format_name}: Matplotlib does not support savefig metadata "
                    "for this format; metadata remains in the optional manifest."
                )
            save_kwargs: dict[str, Any] = {
                "dpi": float(dpi),
                "transparent": transparent,
                "bbox_inches": bbox_inches,
                "pad_inches": float(pad_inches),
                "facecolor": "none" if transparent else facecolor,
                "edgecolor": "none" if transparent else edgecolor,
                **extra_kwargs,
            }
            if format_metadata is not None:
                save_kwargs["metadata"] = format_metadata
            if format_name == "tiff":
                save_kwargs["pil_kwargs"] = {
                    **save_kwargs.get("pil_kwargs", {}),
                    "compression": "tiff_lzw",
                }

            _atomic_savefig(
                fig,
                destination,
                format_name=format_name,
                save_kwargs=save_kwargs,
                tiff_rgb=tiff_rgb,
            )
            outputs.append(
                {
                    "path": str(destination.resolve()),
                    "format": format_name,
                    "kind": (
                        "vector-container"
                        if format_name in VECTOR_FORMATS
                        else "raster"
                    ),
                    "size_bytes": destination.stat().st_size,
                }
            )

    report: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "outputs": outputs,
        "settings": {
            "formats": normalized_formats,
            "dpi": float(dpi),
            "transparent": transparent,
            "tiff_rgb": tiff_rgb,
            "bbox_inches": bbox_inches,
            "pad_inches": float(pad_inches),
            "facecolor": facecolor,
            "font_mode": font_mode,
            "figure_size_inches": figure_size,
        },
        "versions": {
            "matplotlib": mpl.__version__,
            "pillow": _package_version("Pillow"),
        },
        "provenance": provenance_document,
        "warnings": sorted(set(warnings_list)),
        "notice": (
            "Export completed with explicit settings. This report does not "
            "establish visual, accessibility, scientific-integrity, or "
            "journal-submission compliance."
        ),
    }

    if write_manifest:
        manifest_path = parent / f"{base.name}.export.json"
        payload = (
            json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
        ).encode("utf-8")
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{manifest_path.name}.", suffix=".tmp", dir=manifest_path.parent
        )
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        os.chmod(temporary_name, 0o600)
        os.replace(temporary_name, manifest_path)
        report["manifest"] = str(manifest_path.resolve())
    return report
