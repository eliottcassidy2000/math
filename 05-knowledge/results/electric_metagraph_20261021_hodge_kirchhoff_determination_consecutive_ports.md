# Session mac-mini-2026-10-21-electric: the tournament metagraph as an electrical network -- Hodge potentials, Kirchhoff currents, second-order determination, and the consecutive-start-port ledgers

Status: pre-audit. Canon: THM-4739 (PROVED + FINITE-EXACT), HYP-9351 (OPEN). THM-4740 reserved and ceded
unused. Experiments: `04-computation/experiments/electric_metagraph_20261021/` and additions to
`04-computation/experiments/endpoint_layers_20261020/`.

## 1. The owner's lens and the obligations taken up

"Classify tournament metagraph structure thinking in terms of intransitivity of information or electrical
networks; pursue next obligations." Obligations taken from HYP-9345 (g) (the determination problem) and
HYP-9348 (a) (the structural (n-5,4) theorem: the consecutive-port patterns). Concurrent state read first:
THM-4737 was corrected by commit cfcb930 (closed one-sided area criterion, N_min finite iff the threshold is
attained, supremum versus attainment; the recalled fractional-chromatic identity refuted by {1,3,4,7}); this
session builds on the corrected file and does not touch the lonely-runner lane.

## 2. THM-4739

**(H) Hodge.** The +-1 arc flow of a tournament on K_n is gradient + curl with potential
`p_i = (2 s_i - (n-1))/n`; `E_grad = (n^2-1)/3 - 8 c3/n`, `E_curl = (n-1)(n-2)/6 + 8 c3/n`. Scores are the
first-order (transitive) information, the cyclic-triangle 2-form the second-order (intransitive) one; THM-785's
axis is the curl energy. Verified n = 4..7 for every class.

**(K) Kirchhoff.** For an invariant with exact line flux (c3 via THM-785 (4), V via THM-785 (1)) the metagraph
Laplacian applied to the potential equals minus the class-intrinsic current
`J(u) = sum_{Ham paths} (s(start) - s(end) - 1)/|Aut|` (doubled for merged nodes); total current zero;
exact at every node for n = 4..7. Resistance distances (unit conductance per line): Kirchhoff indices 3,
34.8, 157.2, 2610.6; at n = 7 the resistance from the transitive node has Spearman -0.96 with log t_u.

**(D) Determination as a transportation problem.** Port classes `(s(start), s(end), symmetric)` of closed and
open tilings per class; a closed port of type (a,b,sym) is joined to an open port of type (n-a, n-2-b, sym)
of a class carrying the partner invariants computed from the tiling; black lines obey the mirror symmetry
`m[C -> D] = m[C^rev -> D^rev]`. CP-SAT enumeration (cap 2e5):

```text
level L1 (partner scores):                 n = 4: 1 solution;  n = 5: 4 solutions (7 free port classes), the real metagraph among them;
                                           n = 6: >= 200,000 solutions (196 free port classes of 229): first-order data are far from sufficient.
level L2 (+ partner curl profile):         n = 5: every port class forced (30/30), 1 solution = the metagraph;
                                           n = 6: 153 of 355 port classes free, >= 200,000 solutions (cap);
level L3 (+ partner second scores):        n = 6: 37 of 403 port classes free, exactly 48,000 solutions, 3 of them
                                           projecting to the real merged metagraph.
```

So the metagraph is determined by second-order local data at n = 5 and by nothing up to third order at n = 6.

At the class level `(scores, curl profile)` still collides (one pair at n = 5; 11 groups at n = 6; 107 at
n = 7), so the L2 determination at n = 5 is carried by the port structure, not by class invariants.

**(L) Consecutive-start ports.** In the closed endpoint-k layer with `I = {2, ..., k}` (so
`J = {n-1-k, ..., n-3}`), the arc between b and y splits the pattern exactly in half: `y -> b` (epsilon = 0)
is repaired by the packet's closed-no-chord ledger and `b -> y` (epsilon = 1) by the port-block insertion
`(a,b,x,U,z,J,c), (a,x,U,z,J,b,c), (a,x,U,J,b,c,z)`, proved for all n >= k+5 (directedness, path-edge
cancellation, three explicit anti-reversal witnesses). Exhaustive: k = 4 at n = 9, 10, 11, 12 (128, 1024,
16384, 262144 presentations per half) and k = 5 at n = 11 (16384 per half; n = 12 recorded in
consecutive_start_ledgers_large.out). The other consecutive-port patterns do not split by one arc
(residual_cover.out: 3-7 templates per residual half).

## 3. Honest scope and a caught error

The Kirchhoff identity is a tautology once the flux law is granted; its value is the list of intrinsic
currents and the transportation formulation of the determination problem, which turns HYP-9345 (g) into
counts. The port block J was first transcribed as `n-k..n-3`; the exhaustive check at n = 9 failed on the
arc `last(U) -> z` and exposed the off-by-one (`J = n-1-k..n-3`, overlapping I at the centre for odd n) before
anything was written. No all-order determination statement is made.

## 4. Files

electric_metagraph_20261021: `electric_census.py` / `.out`, `electric_n{4..7}.json`, `determination.py` /
`.out`, `determination_L2.py` / `_n5.out` / `_n6.out`, `invariant_collisions.py`.
endpoint_layers_20261020: `two_template_cover.py` / `.out`, `residual_cover.py` / `.out`,
`consecutive_start_ledgers.py` / `_small.out` / `_large.out`.
