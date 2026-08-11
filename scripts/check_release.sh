#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PROOF="$ROOT/proof"
LOG="$ROOT/release-check.log"

# The published CI/container path uses system packages. During local release
# preparation, reuse the workspace dependency cache when it is available.
LOCAL_DEPS="$ROOT/../../../work/deps"
if [[ -d "$LOCAL_DEPS/boost-pkg/root/usr/include" ]]; then
  export CPLUS_INCLUDE_PATH="$LOCAL_DEPS/boost-pkg/root/usr/include:${CPLUS_INCLUDE_PATH:-}"
fi
if [[ -d "$LOCAL_DEPS/mpfr-pkg/root/usr/include" ]]; then
  export CPLUS_INCLUDE_PATH="$LOCAL_DEPS/mpfr-pkg/root/usr/include:$LOCAL_DEPS/mpfr-pkg/root/usr/include/x86_64-linux-gnu:${CPLUS_INCLUDE_PATH:-}"
  export LIBRARY_PATH="$LOCAL_DEPS/mpfr-pkg/root/usr/lib/x86_64-linux-gnu:${LIBRARY_PATH:-}"
  export LD_LIBRARY_PATH="$LOCAL_DEPS/mpfr-pkg/root/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
fi
required=(
  "$ROOT/README.md"
  "$ROOT/AUTHORSHIP.md"
  "$ROOT/DISCLOSURE.md"
  "$ROOT/CITATION.cff"
  "$ROOT/.zenodo.json"
  "$ROOT/LICENSE"
  "$ROOT/LICENSES.md"
  "$ROOT/LICENSES/Apache-2.0.txt"
  "$ROOT/LICENSES/CC-BY-4.0.txt"
  "$ROOT/VERIFICATION.md"
  "$ROOT/RELEASE_CHECKLIST.md"
  "$ROOT/RELEASE_SHA256SUMS"
  "$PROOF/candidate.json"
  "$PROOF/proof_report.md"
  "$PROOF/checksums.sha256"
  "$PROOF/reproduction/logs/final_certificate.txt"
)

for path in "${required[@]}"; do
  [[ -f "$path" ]] || { echo "missing required file: $path" >&2; exit 1; }
done

cd "$ROOT"
sha256sum -c RELEASE_SHA256SUMS

cd "$PROOF"
sha256sum -c checksums.sha256

bash -n reproduction/scripts/verify_existing.sh
bash -n reproduction/scripts/rerun_all.sh
python3 -c "import ast, pathlib; [ast.parse(pathlib.Path(p).read_text()) for p in ('reproduction/scripts/generate_jobs.py','reproduction/scripts/assemble_logs.py')]"
python3 -c "import json, pathlib; [json.loads(pathlib.Path(p).read_text()) for p in ('candidate.json','manifest.json','report_artifact.json')]"
python3 -c "import json, pathlib; json.loads(pathlib.Path('../.zenodo.json').read_text())"

cmp "$ROOT/LICENSE" "$ROOT/LICENSES/Apache-2.0.txt"
if grep -R --line-number --fixed-strings 'REPLACE-WITH' "$ROOT/CITATION.cff" "$ROOT/.zenodo.json" "$ROOT/README.md" "$ROOT/AUTHORSHIP.md"; then
  echo "release metadata placeholder found" >&2
  exit 1
fi

if grep -R --line-number --fixed-strings -- '-ffast-math' reproduction/source reproduction/scripts; then
  echo "forbidden -ffast-math flag found" >&2
  exit 1
fi
grep -q -- '-frounding-math' reproduction/scripts/verify_existing.sh
grep -q -- '-fno-fast-math' reproduction/scripts/verify_existing.sh

set -o pipefail
bash reproduction/scripts/verify_existing.sh | tee "$LOG"
grep -qx 'PROOF_OK 1' reproduction/run_quick/certificate.txt
grep -qx 'certificate verified' <(tail -n 1 "$LOG")
cmp reproduction/run_quick/certificate.txt reproduction/logs/final_certificate.txt

echo "RELEASE_CHECKS_OK"
