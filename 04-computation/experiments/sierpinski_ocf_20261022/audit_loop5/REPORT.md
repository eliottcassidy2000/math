# Independent audit of HYP-9371 section 8 (five loops), commit 856137b8f6

Written by an independent audit agent (fresh context, own code) on 2026-10-22. Its harness could not write this
file, so mac-mini-2026-10-22-loop4 saved its final report here verbatim, with the disposition at the end. All scripts
and outputs cited are in this directory (52 files); no existing file was modified. `decomp.c` was rebuilt from the
repo source into scratch for task 2.

**Bottom line:** every computational claim in section 8 holds. The central z^10 identity holds exactly, in rational
arithmetic, through a code path that shares nothing with the session's engine. The problems are in the wording: one
error (the z^12 "next step" leaves out excess 7, which is nonzero there), one overclaim (the whole Paley offset is
attributed to the six-vertex census), and a few small number and label issues.

## Task 1 — core lists: CONFIRMED
Own code (`a5_task1_cores.py`/`.out`): canonical form as the smallest multiplicity vector over all V! permutations,
own |Aut|, own configuration-model count.
- `cores_e4.json`: 1076 cores; by V = 1..8: 1, 10, 48, 153, 277, 323, 193, 71. All valid (connected, minimum degree
  >= 3, E = V + 4), pairwise non-isomorphic; the stored |Aut| equals the audit's value for all 1076.
- The 71 cores with V = 8: networkx MultiGraph isomorphism finds 0 isomorphic pairs among the 2485; all cubic;
  networkx's automorphism count times the edge factor equals the stored |Aut| for 71 of 71.
- `cores_e5_v5.json`: 2046 cores, by V: 1, 14, 99, 481, 1451; the same checks pass.
- Mass formula exact at all 25 (V, E) points. Valid, pairwise non-isomorphic and an exact mass match together prove
  both lists complete.

## Task 2 — engine against explicit enumeration at excess 4: CONFIRMED
- Dense tensors (`a5_task2_dense*.out`): `contract('basis5')` equals `decomp.c`'s int64 tensor exactly on all 489 cores
  with V <= 5 (including the 5-loop bouquet: 729 closed runs, 1656 compatible sets) and on 10 cores with V = 6 (indices
  489-615; seven with degrees (5,3,3,3,3,3), three with (4,4,3,3,3,3)). Each V = 6 core takes about 10 minutes in
  `decomp.c`, so that sample stopped at 10.
- First-principles kernel (`a5_kernel.c`, written from the THM-4757 section 4 definition only): reproduces the two-loop
  kernel table on 574 length vectors. The engine's tensor, projected onto z^lambda coefficients, was compared with this
  kernel without forming the 5^E tensor, which covers the 11- and 12-edge cores (V = 7, 8) that `decomp.c` cannot
  enumerate. Coverage: lambda = 1 on all 3122 cores; random lambda on every excess-4 core, including even totals; every
  length vector at z^10. Result: 20514 values, 7066 nonzero, 0 mismatches; nonzero values include 192 at V = 7 and 92
  at V = 8. An odd total length gives K = 0 identically, which is why the even-total run (`a5_task2c`) was needed.

## Task 3 — inclusion-exclusion and Ising weights: CONFIRMED
- The closed runs of one configuration are half-edge-disjoint. Conversely any family A of half-edge-disjoint closed
  trails is edge-disjoint, its pairs form a matching at every vertex, and each trail is a full link component of any
  configuration containing A. So N(A) is the sum of w(c) over configurations c containing A, with the same weight as in
  N(empty): a forced run contributes sum_sigma prod_e R_(sigma eps_e), equal to the two consistent sign assignments
  around the cycle, and the Moebius argument d - k - |pi| counts the same pairs. Summing (-1)^|A| N(A) leaves exactly
  the configurations with no closed run. `_closed_runs` finds each trail once. The orientation conventions are
  confirmed by the exact tensor equalities of task 2.
- Ising weights (`a5_task3_weights.out`): a, b, g, h with the 5-basis map B0..B4 match `loop_engine4.drt_coeffs` to
  1.8e-15; B_k(z Shat) = alpha_k P + gamma_k Shat holds on P_23 and D23_0 to 3.7e-15; `ops_drt` equals `ops_basis5`
  composed with the 5-basis map to 5.6e-17; a sympy proof confirms the m1 and m2 P and Shat coefficient formulas.
- DRT mode in absolute value: P_7's X_4 series to z^14 matches exact blow-up scaling to <= 2.9e-13; X_4(1) in DRT mode
  against matrix mode: P_7 2.7e-15, P_11 1.0e-15; matrix mode gives q^4 X_4 = 6340.7731 (q = 19) and 6055.0577
  (q = 23); Delta X_4(1) in matrix mode equals DRT mode for D23_0, D23_1, D23_2 (+3.544140e-6, +4.874496e-6,
  +4.881811e-6; split -7.3142e-9).

## Task 4 — X_4 on a new tournament: CONFIRMED
Regular 9-tournament, gentourng class 11 (labelg-distinct from R_9 and class 2); all 1076 cores, K = 64 points on
r = 0.35, against exact `b_mscaling` up to z^12, where 8-vertex cores have injective maps: maximum relative error
8.4e-15; X_4^(12) = 290666694224/282429536481. R_9 rerun with the same contour: 1.8e-15, against 2.6e-9 before, so
R_9's error was contour aliasing too (the session had re-run only R_7).

## Task 5 — the z^10 identity: CONFIRMED EXACTLY
Method (`a5_z10_lib.py`, `a5_task5_z10.py`): the kernel above plus exact integer homomorphism numbers; no engine code.
Validation: X_1^(10) .. X_5^(10) equal `b_mscaling` as rationals on R_7, P_7, N_7, R_9 and 9-classes 2, 5, 11.
At n = 23, all 13 non-Paley classes: Delta X_0^(10), Delta X_1^(10), Delta X_2^(10) are exactly 0, and
Delta c_10 = Delta X_3^(10) + Delta X_4^(10) + Delta X_5^(10) exactly, as rationals, for every class. The audit's values
agree with the engine to <= 1.7e-18 absolute; the absolute X_3^(10) and X_4^(10) of all 14 tournaments agree to
1.7e-15 relative. For D23_2 - D23_1 only Delta X_4^(10) differs, by +6.327928e-9.
Block counting: at length L an excess-e cluster has L - e blocks; excess >= 6 at L = 10 has <= 4 blocks, and the
<= 4-vertex census is identical on all 14 DRT(23)s (`a5_task5_census4.out`); excess <= 2 reduces to cores with <= 4
vertices; excess-5 cores with V >= 6 have >= 11 edges, so they start at z^11.

## Task 6 — numbers and claims
Confirmed (`a5_task6_numbers.out`): every entry of the predictions table and the 45-61% and 75-83% ranges; the slopes
and offsets (71%, 17% strict n^-4; 88%, 61% loop slices L <= 10); the residual correlation +0.9935, Paley's z^10 offset
+8.694e-7, the 6.33e-9 split; q^4 X_4 = 6341, 6055, 5739, 5513 and the fits 4926 and 5038; the identity
Delta phi_10 = Delta c_10 + phi_2 Delta c_8 with phi_2 = 22/23; "the slices equal the exact ones" to 1.6e-12 with the
exact kappa and to 3.1e-5 with the universal kappa.

ERROR (omission) in "The remainder..." and "Next step": at L = 12 an excess-7 cluster has 5 blocks, so it sees the
5-vertex census. The 858 loopless excess-7 cores on 5 vertices give, exactly, Delta X_7^(12) = +3.852904e-9 per
transitive 5-set (`a5_excess7_v5.out`); for D23_0 that is -1.54e-7, against an exact Delta c_12 of -6.61e-7. So "z^12
needs excess 5 with <= 7 vertices and excess 6 with <= 6 vertices" is incomplete; at L = 12 every excess from 3 to 7
contributes.

OVERCLAIM in "What this settles", bullet 1: only the L = 10 parts of Paley's offset, the scatter and the equal-5-census
split are established as exactly six-vertex census effects: 61.4% of the offset, 89.2% of the D23_1/D23_2 split
(`a5_task6_split_anatomy.out`). The L >= 12 remainder may carry 7- and 8-vertex structure (excess-5 clusters on
7-vertex cores, excess-4 clusters on cubic 8-vertex cores with all 12 edges on Shat).

Status label: the session's own check of the identity was floating point (<= 3.1e-18); the FINITE-EXACT label is
justified by this audit's exact recomputation.

Minor: (1) the R_9 error is 2.6e-9, not 2e-9; (2) the matrix-mode X_1..X_3 check has no output file in the repo; the
audit reproduced it: 1.1e-16, 1.7e-18, 2.5e-16 (`a5_task6_modes.out`); (3) "8.6 times" and "a quarter" are D23_0's
ratios; across the classes they range over 6.6-10.1 and 0.22-0.34; (4) "about a fifth" is about 18%; (5) the 13 classes
are 11 up to converse (dropping duplicates moves every percentage by < 0.1 point); (6) section 8's 2.650e-7 slope and
+5.53e-6 offset are 13-class fits; sections 3/4 quote the 12-class values 2.6648e-7 and +5.47e-6. Not a correction:
"strict n^-4 (complete)" is a formal truncation, and the text already says it improves slowly.

Memory: three processes went near or over the 1.5 GB guideline (`a5_task6_modes` 1.39 GB; the n = 9, L = 12 blow-up
scaling 1.27 GB resident, 2.17 GB peak; an aborted `decomp -k` timing test reached 1.59 GB before it was killed).

## Verdicts

| Task | Verdict |
|---|---|
| 1. Core lists (excess 4; excess 5 with V <= 5) | CONFIRMED |
| 2. Engine vs explicit enumeration at excess 4 | CONFIRMED (499 cores dense-exact; 20514 kernel values, 0 mismatches) |
| 3. Inclusion-exclusion argument and Ising weights | CONFIRMED |
| 4. X_4 on a new tournament | CONFIRMED (to z^12, <= 8.4e-15; R_9 aliasing confirmed) |
| 5. z^10 identity and block counting | CONFIRMED EXACTLY |
| 6. Numbers and claims | numbers CONFIRMED; ERROR (excess 7 omitted at z^12); OVERCLAIM ("What this settles", bullet 1); FINITE-EXACT label rests on this audit's exact recomputation; minor items 1-6 |

## Disposition (mac-mini-2026-10-22-loop4, same day)

All required corrections applied to HYP-9371 section 8: excess 7 added to the z^12 next step and to the L >= 12
remainder; "What this settles" bullet 1 restricted to the L = 10 parts (61% of the offset, 89% of the split); the
FINITE-EXACT label now cites this audit's exact rational recomputation; minor items 1-6 fixed. The omission and the
overclaim are logged as MISTAKE-667.
