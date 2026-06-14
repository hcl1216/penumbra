"""
Penumbra -- central path configuration.

Single source of truth for every filesystem path used by the Phase-0 scripts.
NO script may hard-code a drive letter or absolute path; import from here instead.

    from config import PATHS
    slide_dir = PATHS["cosmx_slide"]

Design rules (see CLAUDE.md "Infrastructure / conventions"):
  - Cross-platform: built with pathlib, no hard-coded "C:\\..." anywhere.
  - Windows-local (authoring/CPU) and Colab-Linux (GPU) must both work unchanged.
  - Every path is overridable by an environment variable so the same code runs
    against a local copy, a Colab-mounted Drive, or CI without edits.

Environment variables (override the defaults below):
  PENUMBRA_DRIVE_ROOT   root of the (Colab-mounted) Google Drive holding raw data
  PENUMBRA_COSMX_SLIDE  dir of the CosMx breast multiomic demo slide (raw flat files)
  PENUMBRA_METABRIC     dir of the METABRIC-IMC (Danenberg 2022) + bulk RNA (Curtis 2012)
  PENUMBRA_DATA_DIR     small/derived artifacts (defaults to <repo>/data)
  PENUMBRA_RESULTS_DIR  gate outputs / figures / verdicts (defaults to <repo>/results)

Nothing here reads or writes data; this module only resolves locations.
"""

from __future__ import annotations

import os
from pathlib import Path

# Repo root = directory containing this file. Never a hard-coded drive letter.
REPO_ROOT = Path(__file__).resolve().parent


def _env_path(var: str, default: Path) -> Path:
    """Return Path from env var `var` if set, else `default`. Pure resolution."""
    val = os.environ.get(var)
    return Path(val).expanduser() if val else default


# Raw data lives off-repo (Google Drive / Colab mount); these are placeholders
# until the real locations are recorded. Override via env vars above.
DRIVE_ROOT = _env_path("PENUMBRA_DRIVE_ROOT", REPO_ROOT / "_drive")

PATHS: dict[str, Path] = {
    "repo_root": REPO_ROOT,
    "drive_root": DRIVE_ROOT,
    # Phase-0 blocking input -- raw CosMx breast multiomic demo slide (flat files).
    "cosmx_slide": _env_path("PENUMBRA_COSMX_SLIDE", DRIVE_ROOT / "cosmx_breast_multiomic"),
    # Gate-E / later -- METABRIC-IMC (Danenberg 2022) + METABRIC bulk RNA (Curtis 2012).
    "metabric": _env_path("PENUMBRA_METABRIC", DRIVE_ROOT / "metabric_imc"),
    # Panel target lists used by Step 0 (panel overlap).
    "panels": _env_path("PENUMBRA_PANELS", REPO_ROOT / "data" / "panels"),
    # Small / derived artifacts kept in-repo-adjacent (gitignored).
    "data_dir": _env_path("PENUMBRA_DATA_DIR", REPO_ROOT / "data"),
    # Gate outputs, figures, verdicts (gitignored; mirrored to Drive).
    "results_dir": _env_path("PENUMBRA_RESULTS_DIR", REPO_ROOT / "results"),
}
