#!/usr/bin/env python3
"""Validate and assemble 244 single-segment interval logs."""
from decimal import Decimal, getcontext
from pathlib import Path
import sys

getcontext().prec = 100
if len(sys.argv) != 4:
    raise SystemExit("usage: assemble_logs.py JOB_DIR CENTER_OUT OUTER_OUT")
job_root = Path(sys.argv[1])
outputs = {"center": Path(sys.argv[2]), "outer": Path(sys.argv[3])}
expected_tokens = (4, 6, 17, 17, 129)
bound = Decimal("0.003")

for kind, destination in outputs.items():
    records = []
    minimum = Decimal("Infinity")
    for i in range(244):
        path = job_root / kind / f"{i:03d}.log"
        lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        if len(lines) != 6:
            raise SystemExit(f"{path}: expected five record lines plus footer, got {len(lines)}")
        rec = lines[:5]
        lengths = tuple(len(line.split()) for line in rec)
        if lengths != expected_tokens:
            raise SystemExit(f"{path}: token counts {lengths}, expected {expected_tokens}")
        if rec[0] != "REC 0 OK 1":
            raise SystemExit(f"{path}: failed record {rec[0]!r}")
        footer = lines[5].split()
        if footer[:5] != ["BATCH_OK", "1", "COUNT", "1", "GLOBAL_D2"]:
            raise SystemExit(f"{path}: malformed footer")
        d2 = Decimal(rec[1].split()[3])
        footer_d2 = Decimal(footer[5])
        if d2 <= bound or footer_d2 <= bound:
            raise SystemExit(f"{path}: collision bound is not above {bound}")
        minimum = min(minimum, d2, footer_d2)
        rec[0] = f"REC {i} OK 1"
        records.extend(rec)
    records.append("BATCH_OK 1 COUNT 244 GLOBAL_D2 3.00000000000000000000e-03")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(records) + "\n")
    print(f"{kind}: 244 records; raw minimum d^2 lower bound {minimum}; archived bound 0.003")
