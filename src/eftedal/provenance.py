"""Provenance stamps for generated results: every pipeline script records which commit (and whether
the code had uncommitted changes) produced its outputs, in results/PROVENANCE/<script>.txt. Used by
ASHG/build.sh to refuse to rebuild the ASHG 2026 materials from results made by newer code."""

from __future__ import annotations

import datetime as dt
import subprocess
from pathlib import Path

from .data import REPO_ROOT

CODE_PATHS = ("src", "scripts", "R")


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO_ROOT), *args], capture_output=True, text=True).stdout.strip()


def stamp(script: str | Path) -> Path:
    """Write results/PROVENANCE/<script name>.txt: commit, whether src/scripts/R were dirty, time."""
    out = REPO_ROOT / "results" / "PROVENANCE"
    out.mkdir(parents=True, exist_ok=True)
    dirty = bool(_git("status", "--porcelain", "--", *CODE_PATHS))
    path = out / f"{Path(script).stem}.txt"
    path.write_text(f"commit={_git('rev-parse', 'HEAD')}\ndirty={int(dirty)}\n"
                    f"time={dt.datetime.now().isoformat(timespec='seconds')}\n")
    return path
