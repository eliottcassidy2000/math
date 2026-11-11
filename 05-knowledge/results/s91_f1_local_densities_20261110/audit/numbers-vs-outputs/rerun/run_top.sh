#!/usr/bin/env bash
# Re-run of the top-level fast scripts in an isolated copy (auditor numbers-vs-outputs).
# At most 3 threads.
export NUMBA_NUM_THREADS=3
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3
cd "$(dirname "$0")/top" || exit 2
run() {
  name="$1"; shift
  start=$(date +%s)
  python3 "$@" > "$name.rerun.out" 2> "$name.rerun.err"
  code=$?
  end=$(date +%s)
  echo "$name exit=$code seconds=$((end-start))" >> ../rerun_status_top.txt
}
: > ../rerun_status_top.txt
run constants constants.py
run mod4_symbol mod4_symbol.py
run box_identity box_identity.py
run ap_constants ap_constants.py
run regions_check regions_check.py
run twisted_spectrum twisted_spectrum.py 16
echo DONE >> ../rerun_status_top.txt
