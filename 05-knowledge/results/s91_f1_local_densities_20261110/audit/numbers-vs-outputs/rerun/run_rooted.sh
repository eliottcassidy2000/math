#!/usr/bin/env bash
# Re-run of the rooted_spectrum scripts in an isolated copy (auditor numbers-vs-outputs).
# At most 3 threads.
export NUMBA_NUM_THREADS=3
export OMP_NUM_THREADS=3
export MKL_NUM_THREADS=3
export OPENBLAS_NUM_THREADS=3
cd "$(dirname "$0")/rooted_spectrum" || exit 2
run() {
  name="$1"; shift
  start=$(date +%s)
  python3 "$@" > "$name.rerun.out" 2> "$name.rerun.err"
  code=$?
  end=$(date +%s)
  echo "$name exit=$code seconds=$((end-start))" >> ../rerun_status_rooted.txt
}
: > ../rerun_status_rooted.txt
run source_supply source_supply.py
run head_parity_check head_parity_check.py
run source_witnesses source_witnesses.py
run line_density_order9 line_density_order10.py 8
run line_density_order10 line_density_order10.py 9
echo DONE >> ../rerun_status_rooted.txt
