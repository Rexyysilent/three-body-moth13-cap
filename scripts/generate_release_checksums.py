#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "RELEASE_SHA256SUMS"

EXCLUDED_DIRS = {
    ROOT / ".git",
    ROOT / "proof" / "reproduction" / "run_quick",
    ROOT / "proof" / "reproduction" / "run_full",
}
EXCLUDED_FILES = {
    OUTPUT,
    ROOT / "release-check.log",
}


def excluded(path: Path) -> bool:
    if path in EXCLUDED_FILES or path.suffix in {".zip", ".pyc"}:
        return True
    return any(directory == path or directory in path.parents for directory in EXCLUDED_DIRS)


lines: list[str] = []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file() or excluded(path):
        continue
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    lines.append(f"{digest}  {path.relative_to(ROOT).as_posix()}")

OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
print(f"wrote {OUTPUT} with {len(lines)} entries")
