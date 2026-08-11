#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT.parent / f"{ROOT.name}.zip"
ARCHIVE_HASH = ARCHIVE.with_suffix(ARCHIVE.suffix + ".sha256")
FIXED_TIME = (2026, 8, 11, 0, 0, 0)

subprocess.run(["sha256sum", "-c", "RELEASE_SHA256SUMS"], cwd=ROOT, check=True)

listed = []
for line in (ROOT / "RELEASE_SHA256SUMS").read_text(encoding="utf-8").splitlines():
    _digest, relative = line.split("  ", 1)
    listed.append(Path(relative))
listed.append(Path("RELEASE_SHA256SUMS"))

with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
    for relative in sorted(listed, key=lambda item: item.as_posix()):
        source = ROOT / relative
        info = zipfile.ZipInfo(f"{ROOT.name}/{relative.as_posix()}", FIXED_TIME)
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = (0o100644 & 0xFFFF) << 16
        output.writestr(info, source.read_bytes())

digest = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
ARCHIVE_HASH.write_text(f"{digest}  {ARCHIVE.name}\n", encoding="utf-8", newline="\n")
print(f"wrote {ARCHIVE}")
print(f"sha256 {digest}")
