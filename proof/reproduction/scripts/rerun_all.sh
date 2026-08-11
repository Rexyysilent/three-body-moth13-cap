#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/run_full/build"
RUN="$ROOT/run_full"
SRC="$ROOT/source"
DATA="$ROOT/data"
JOBS="${JOBS:-$(nproc)}"
mkdir -p "$BUILD" "$RUN"

CXX="${CXX:-g++}"
CXXFLAGS=(-O2 -std=c++17 -frounding-math -fno-fast-math)
CPPFLAGS=()
LDFLAGS=()
LOCAL_DEPS="$ROOT/../../work/deps"
if [[ -d "$LOCAL_DEPS/boost-pkg/root/usr/include" ]]; then
  CPPFLAGS+=("-I$LOCAL_DEPS/boost-pkg/root/usr/include")
fi
if [[ -d "$LOCAL_DEPS/mpfr-pkg/root/usr/include" ]]; then
  CPPFLAGS+=("-I$LOCAL_DEPS/mpfr-pkg/root/usr/include" "-I$LOCAL_DEPS/mpfr-pkg/root/usr/include/x86_64-linux-gnu")
  LDFLAGS+=("-L$LOCAL_DEPS/mpfr-pkg/root/usr/lib/x86_64-linux-gnu")
  export LD_LIBRARY_PATH="$LOCAL_DEPS/mpfr-pkg/root/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
fi

"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/boost_taylor_batch.cpp" -o "$BUILD/outer_flow"
"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/boost_taylor_batch_mpfi.cpp" "${LDFLAGS[@]}" -lmpfi -lmpfr -lgmp -o "$BUILD/center_flow"
"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/krawczyk244.cpp" -o "$BUILD/krawczyk244"
"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/krawczyk_mpfi_terms.cpp" "${LDFLAGS[@]}" -lmpfi -lmpfr -lgmp -o "$BUILD/krawczyk_mpfi_terms"

python3 "$ROOT/scripts/generate_jobs.py" "$DATA/mesh244.certmesh" "$RUN/jobs"

find "$RUN/jobs/center" -name '*.in' -print0 | sort -z |   xargs -0 -P "$JOBS" -I{} bash -c 'f="$1"; exe="$2"; "$exe" < "$f" > "${f%.in}.log" 2> "${f%.in}.err"' _ {} "$BUILD/center_flow"
find "$RUN/jobs/outer" -name '*.in' -print0 | sort -z |   xargs -0 -P "$JOBS" -I{} bash -c 'f="$1"; exe="$2"; "$exe" < "$f" > "${f%.in}.log" 2> "${f%.in}.err"' _ {} "$BUILD/outer_flow"

python3 "$ROOT/scripts/assemble_logs.py" "$RUN/jobs" "$RUN/mpfi_center.log" "$RUN/outer_box.log"

rm -f "$RUN"/preconditioner.cert*
set +e
"$BUILD/krawczyk244" "$DATA/mesh244.certmesh" "$RUN/mpfi_center.log" "$RUN/outer_box.log" "$RUN/preconditioner.cert"
point_rc=$?
set -e
if [[ "$point_rc" -ne 0 && "$point_rc" -ne 10 ]]; then
  echo "point-preconditioner stage failed with exit code $point_rc" >&2
  exit "$point_rc"
fi

"$BUILD/krawczyk_mpfi_terms"   "$DATA/mesh244.certmesh" "$RUN/mpfi_center.log" "$RUN/outer_box.log"   "$RUN/preconditioner.cert.C.bin" "$RUN/preconditioner.cert.point.bin" "$RUN/preconditioner.cert.R.bin"   "$RUN/certificate.txt"

grep -qx 'PROOF_OK 1' "$RUN/certificate.txt"
cat "$RUN/certificate.txt"
echo "full enclosure and certificate rerun passed"
