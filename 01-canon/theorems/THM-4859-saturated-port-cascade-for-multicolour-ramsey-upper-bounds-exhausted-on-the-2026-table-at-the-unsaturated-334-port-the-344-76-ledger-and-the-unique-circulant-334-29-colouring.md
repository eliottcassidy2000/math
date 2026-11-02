---
id: THM-4859
title: "The saturated-port cascade for multicolour Ramsey upper bounds: if R(k) equals its recursive bound P(k) and the nested colour-i ports are saturated to depth t-2, the number of colour-i cliques K_t is the forced integer n_0 n_1 ... n_(t-1)/t! (t = 2 is the classical parity strictness, t = 3 is Boza's mod-3 criterion), each such condition yields exactly one unit; on the April 2026 table the cascade reproduces Boza's R(3,3,6) <= 91, R(3,4,5) <= 157, R(4,4,4) <= 229 and is otherwise exhausted, dying for R(3,4,4) = 77 and R(3,3,5) = 57 at the unsaturated (3,3,4;29) port (30 < 34), so the targets 73, 137, 86, 198 require six, eleven and twenty-one units of direct exclusion; the (3,4,4;76) saturation ledger (Paley17 ports, |N1(r) ∩ N2(v)| = 8, e1(N2(v)) = e1(N3(v)), odd colour-1 cross count); and the cyclic (3,3,4;29) colouring is unique up to multipliers and the colour swap, with split (116,116,174)"
status: >
  PROVED (klein-2026-10-30, elementary; not independently audited): (A) the depth-t cascade, generalising CITED
  Boza 2026 Theorem 2.1 (t = 3) and the classical parity strictness (t = 2) of Greenwood-Gleason's recursion; (C) the
  (3,4,4;76) ledger. FINITE-EXACT given CITED inputs (Radziszowski DS1.18, April 2026, sections 2.1 and 6.1; Boza
  arXiv:2603.10851; Codish-Frank-Itzhakov-Miller 2016 for R(3,3,4) = 30; McKay's database for the unique (4,4;17) graph
  Paley17): (B) the scan over all three-colour tuples with entries 3..7 and four-colour tuples with entries 3..5, and
  the target arithmetic; (D) the circulant census over 3^14 difference colourings. NO NEW RAMSEY BOUND is claimed; the
  content is the exact location of the barrier and the exact requirements of the owner's targets. The "R(3,6,6) to 86"
  of the directive is read as R(3,3,6) (338 <= R(3,6,6) by DS1.18).
session: klein-2026-10-30 (codex/session-ramsey-20261030)
source: 05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/README.md
depends_on:
  - 05-knowledge/results/ramsey_saturated_ports_20261030/PROOFS.md
related:
  - 05-knowledge/hypotheses/HYP-9450-rigidity-of-334-29-ports-exclusion-of-344-76-and-335-56-and-what-the-ramsey-targets-73-137-86-198-would-require.md
  - 01-canon/theorems/THM-4840-tutte-camion-dichotomy-for-arithmetic-square-class-banks-hamming-15-is-the-equal-parity-degree-space-of-k6-rank-6-prime-banks-unit-cayley-tournaments-c3-t5-p7-and-fano-chromatic-numbers.md
scripts:
  - 05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/scripts/ramsey_cascade.py (+ data/ramsey_cascade.out, receipts/)
  - 05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/scripts/circ_probe_334.py (+ data/circ_probe_334.out, receipts/)
  - 05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/scripts/port_probe_334.py (+ data/, solver-dependent samples)
---

# THM-4859 — the saturated-port cascade and where it dies

Notation. For `k = (k_1, …, k_r)` and valid upper bounds `U` (exact values where known), `P(k) = 2 − r + Σ_i
U(k^{(i)})`, `k^{(i)}` lowering `k_i` by one, so `R(k) ≤ P(k)`. An `(k; n)` colouring is an `r`-colouring of `K_n` with no
colour-`i` `K_{k_i}`.

## (A) The cascade

Assume `R(k) = P(k)`. Then `R(k^{(i)}) = U(k^{(i)})` for all `i`, and an `(k; P−1)` colouring is saturated: every
vertex has colour-`i` degree `n_1^{(i)} = U(k^{(i)}) − 1`. Fix `i` with `k_i ≥ 4`, write `k^{(i,j)}` for `k_i` lowered by
`j`, and suppose `U(k^{(i,j)}) = P(k^{(i,j)})` for `1 ≤ j ≤ t−2`. Then the nested colour-`i` ports have forced sizes
`n_j = U(k^{(i,j)}) − 1` and forced colour-`i` degrees `n_{j+1}`, so the colour-`i` cliques `K_t` number exactly
`n_0 n_1 ⋯ n_{t−1}/t!`, `n_0 = P − 1`. If this is not an integer, `R(k) ≤ P(k) − 1`. Proof: a saturated colouring
restricted to a saturated port is again saturated, by the same forcing one level down (`R(k^{(i,1)}) = P(k^{(i,1)})`
forces its terms to attain their bounds); count ordered `t`-cliques by choosing a vertex, then a colour-`i` neighbour in
its port, and so on. `t = 2` is the classical parity rule, `t = 3` is Boza's Theorem 2.1. ∎

The method yields one unit per application and cannot be iterated: on `P − 2` vertices saturation fails.

## (B) Exhaustion on the 2026 table

With DS1.18's values and bounds (`R(3,3,4) = 30`, `R(3,3,5) ≤ 57`, `R(3,4,4) ≤ 77`, `R(3,3,6) ≤ 91`, `R(3,4,5) ≤ 157`,
`R(4,4,4) ≤ 229`, `R(3,3,3,3) ≤ 62`, `R(3,3,3,4) ≤ 149`, …) the scan reproduces Boza: `73892/3`, `57148/3` and
`173014/3`, `147934/3` are the non-integral forced triangle counts behind `R(3,3,6) ≤ 91`, `R(3,4,5) ≤ 157`,
`R(4,4,4) ≤ 229`. No other condition fires on a tuple with a non-recursive bound. **For `R(3,4,4) = 77` and
`R(3,3,5) = 57` the first colour-4 (resp. colour-5) port is a `(3,3,4;29)` colouring with `U = 30 < P = 34`:
unsaturated, so no count is forced and the cascade stops.** Consequently, by recursion arithmetic alone: `R(3,4,4) ≤ 73`
has no recursive route (`77 = 18 + 30 + 30 − 1`, all terms exact); `R(3,3,6) ≤ 86` requires `R(3,3,5) ≤ 51`;
`R(3,4,5) ≤ 137` requires `R(3,3,5) + R(3,4,4) ≤ 113`; `R(4,4,4) ≤ 198` requires `R(3,4,4) ≤ 66`.

## (C) The `(3,4,4;76)` ledger

In a `(3,4,4;76)` colouring every vertex has degrees `(17, 29, 29)`; each `N_1(v)` is Paley17 in colours 2/3 (68 + 68
edges, 8-regular each); `|N_1(r) ∩ N_2(v)| = |N_1(r) ∩ N_3(v)| = 8` for `r ∈ N_1(v)`; `e_1(N_1(v), N_2(v)) = e_1(N_1(v),
N_3(v)) = 136`; `e_1(N_2(v)) = e_1(N_3(v)) = e^*(v)` with `2e^*(v) + e_1(N_2(v), N_3(v)) = 357`, so the colour-1 cross
count between the two 29-ports is odd and at least 125 (in-port colour-1 degrees are at most 8); the mean colour-1
codegree along colour-2 edges is `136/29`. Nothing forces `e_2(N_2(v))`, which a mod-3 contradiction would need. ∎

## (D) The cyclic `(3,3,4;29)` colouring

Among the `3^14` colourings of the differences of `Z/29` exactly 14 are `(3,3,4)` colourings, all with class sizes
`(4, 4, 6)`, hence edge counts `(116, 116, 174)` and degrees `(8, 8, 12)`; they form one orbit under the multiplier
group and the swap of the two triangle-free colours (`D_0 = {1,3,7,12}`, `D_1 = {2,5,6,14}`, `D_2 = {4,8,9,10,11,13}`).
Both triangle-free colours sit at their degree cap 8. ∎ (Presumably Kalbfleisch's 1966 construction.)

Scope. The theorem explains why the four targets are unreachable by recursion-type arguments and what each would cost;
it proves no bound. The enabling datum for any further progress is the set of `(3,3,4;29)` colourings (HYP-9450).
