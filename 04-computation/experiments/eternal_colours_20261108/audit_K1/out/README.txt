audit K1 (2026-11-08) outputs. All produced by the auditor's programs in ../ (see ../aud_run.sh).

census_vals_n.txt      n <= 9: "n h SCclasses strongSCclasses allclasses" from aud_census.c on THM-4922's census
aud_census.log         per-order summary: classes, self-converse classes, field mismatches (0), orbit-stabilizer sums
anti_n_full.txt        aud_sc.c, all sigma-anti-invariant labellings of order n: "h count strongcount"
anti_n_red.txt         the same with within-pair arcs fixed (reduced model); anti_10_red.txt is order 10
aud_compare_le9.out    census vs anti-invariant model, orders 1..9 (identical incl. strong flags)
aud_compare_o10.out    order 10: auditor vs packet sc_even10.out vs THM-4922 audit aud_o10_byH.txt
o11_full.txt           order 11, ALL self-converse classes (reduced model, 2^25 codes, no pruning): "h count strongcount"
o11.log                worker logs
o12_lt625.txt          order 12, all self-converse tournaments with H < 625 (aud_sc12.c): "h count strongcount"
o12.log                its log (27,784 parents; 1,410,552 order-11 intermediates; 8,036 split checks)
pair_n.txt, pair_10_lt625.txt, pair_validation.log
                       aud_scpair.c validation: identical to anti_n_red.txt (n = 4, 6, 8, 10; and below 625 at n = 10)
aud_lists.out          aud_lists.py on the auditor's data with --no12 (the packet's setting)
aud_lists_with12.out   aud_lists.py including order 12 below 625 (decides 377)
aud_witness.out        aud_witness.c on the packet's sc_witness.out
aud_275.out            aud_275.py: 275 at order 13, 121 at order 10
aud_find11.out         aud_find11.c: explicit strong self-converse 11-tournaments with H = 1551 and 1627
aud_find11_check.out   their independent re-check (H, backtracking self-converseness, strongness)
pair12_lt625.txt       order 12 below 625 recomputed with the parametrized aud_scpair.c (identical to o12_lt625.txt)
aud_refine.out         aud_refine.py (independent refinement table)
pk_refine.out          the packet's refine_dichotomy.py rerun on the regenerated census (byte-identical to the packet's)
