#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/run_quick/build"
RUN="$ROOT/run_quick"
SRC="$ROOT/source"
DATA="$ROOT/data"
mkdir -p "$BUILD"

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

"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/krawczyk244.cpp" -o "$BUILD/krawczyk244"
"$CXX" "${CXXFLAGS[@]}" "${CPPFLAGS[@]}" "$SRC/krawczyk_mpfi_terms.cpp" "${LDFLAGS[@]}" -lmpfi -lmpfr -lgmp -o "$BUILD/krawczyk_mpfi_terms"

rm -f "$RUN"/preconditioner.cert*
set +e
"$BUILD/krawczyk244" "$DATA/mesh244.certmesh" "$DATA/mpfi_center.log" "$DATA/outer_box.log" "$RUN/preconditioner.cert"
point_rc=$?
set -e
if [[ "$point_rc" -ne 0 && "$point_rc" -ne 10 ]]; then
  echo "point-preconditioner stage failed with exit code $point_rc" >&2
  exit "$point_rc"
fi

"$BUILD/krawczyk_mpfi_terms"   "$DATA/mesh244.certmesh" "$DATA/mpfi_center.log" "$DATA/outer_box.log"   "$RUN/preconditioner.cert.C.bin" "$RUN/preconditioner.cert.point.bin" "$RUN/preconditioner.cert.R.bin"   "$RUN/certificate.txt"

grep -qx 'PROOF_OK 1' "$RUN/certificate.txt"
cat "$RUN/certificate.txt"
echo "certificate verified"
