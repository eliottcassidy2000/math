#!/bin/sh
set -eu
cd "$(dirname "$0")"
NS='4096 4097 5997 8191 8192 8193 11001 11002 11003 16383 16384 16385 18274 18275 18276 32768 65536 131072 262144'
python scan.py $NS > scan.normal.log
python -O scan.py $NS --suffix .optimized > scan.optimized.log
python literal_verify.py $NS > literal.normal.log
python -O literal_verify.py $NS --suffix .optimized > literal.optimized.log
python mutation_controls.py > mutations.log
python -O mutation_controls.py --suffix .optimized > mutations.optimized.log
for N in $NS; do
    cmp "results/receipt_N$N.json" "results/receipt_N$N.optimized.json"
    cmp "results/literal_receipt_N$N.json" "results/literal_receipt_N$N.optimized.json"
done
cmp results/mutation_receipt.json results/mutation_receipt.optimized.json
python mechanisms.py > mechanisms.log
python -O mechanisms.py --suffix .optimized > mechanisms.optimized.log
cmp results/N4096_clock_border_mechanism.json results/N4096_clock_border_mechanism.optimized.json
cmp results/resource_port_diagnostics.json results/resource_port_diagnostics.optimized.json
python summarize.py > summary.log
cmp results/N4096_depth_boundary_stability_hostile.json results/N4096_depth_boundary_stability_hostile.optimized.json
python pair_check.py 65536 > pair_birth_N65536.log
python -O pair_check.py 65536 --suffix .optimized > pair_birth_N65536.optimized.log
cmp results/pair_birth_receipt_N65536.json results/pair_birth_receipt_N65536.optimized.json
python boundary_invoice.py > boundary_invoice.log
python -O boundary_invoice.py --suffix .optimized > boundary_invoice.optimized.log
cmp results/boundary_invoice_receipt.json results/boundary_invoice_receipt.optimized.json
