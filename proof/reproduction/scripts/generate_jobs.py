#!/usr/bin/env python3
"""Generate one rigorous flow-enclosure job per shooting segment."""
from pathlib import Path
import sys

if len(sys.argv) != 3:
    raise SystemExit("usage: generate_jobs.py MESH OUTPUT_DIR")

mesh_path = Path(sys.argv[1])
out = Path(sys.argv[2])
lines = [line.strip() for line in mesh_path.read_text().splitlines() if line.strip()]
header = lines[0].split()
n = int(header[0])
tau = header[1]
if n != 244:
    raise SystemExit(f"expected 244 segments, got {n}")
den = lines[1].split()
nodes = lines[2:]
if len(den) != n or len(nodes) != n + 1:
    raise SystemExit("malformed mesh")
if den.count("2.000000000000000000000000000000000000000000000000000000000000000000000000000e+02") != 198:
    raise SystemExit("unexpected 1/200 mesh count")
if den.count("3.200000000000000000000000000000000000000000000000000000000000000000000000000e+03") != 30:
    raise SystemExit("unexpected 1/3200 mesh count")
if den.count("2.560000000000000000000000000000000000000000000000000000000000000000000000000e+04") != 16:
    raise SystemExit("unexpected 1/25600 mesh count")

for kind in ("center", "outer"):
    (out / kind).mkdir(parents=True, exist_ok=True)

for i in range(n):
    center_order = 20 if i == 148 else 24
    center_tol = "1e-20" if i in (148, 149) else "1e-25"
    center = f"0 {tau} 0 0 {center_order} {center_tol} 1\n{den[i]} {nodes[i]}\n"
    (out / "center" / f"{i:03d}.in").write_text(center)

    outer_tol = "3e-10" if i >= 217 else "1e-10"
    outer = f"1 {tau} 1.001e-13 1.001e-14 20 {outer_tol} 1\n{den[i]} {nodes[i]}\n"
    (out / "outer" / f"{i:03d}.in").write_text(outer)

print(f"generated {2*n} jobs in {out}")
