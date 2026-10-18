---
id: HYP-9287
title: "Does the class number decide whether the natural order ranks the Paley tournament optimally? REFUTED in the threshold form (natural order optimal iff L(1, chi_(-p)) exceeds a slowly rising threshold): rotation (Beatty) orders x -> {alpha x} beat the natural order at high-L primes 167, 263, 383, 439, 503, 599, 647, 887, 911, 983, while it stays unbeaten only at 47, 71, 191, 239, 311, 359, 431, 479 (full rotation family, audit); surviving question: is the natural order optimal at only finitely many primes (a largest anomaly), or on a thin infinite set governed by L(1, chi_(-p)) against its small-conductor twists?"
status: >
  REFUTED as first stated (independent S40 audit, 04-computation/experiments/s40_audit_20261017/REPORT.md, MISTAKE-626).
  Both pre-registered tests were genuine (pushed before the runs; reruns byte-identical) and passed as executed, but the
  pre-registered search (12-30 restarts of insertion local search) was too weak for its "unbeaten" claims: certified
  rotation orders beat the natural order at 167, 263, 383 (and 439), so P1 has 4 of 10 high-L primes beaten, the failure
  its own rule names; in the second test 503, 599, 647 and 887 (audit) and 911 (S40 Farey sweep) are beaten, so Q1 fails (5 of 7; at most one allowed). P2 and Q2
  (low L always beaten, certified witnesses) stand. PROVED (finite, via Heegner-Baker-Stark): never optimal for
  h(-p) = 1, p > 7 (11, 19 by exact MAS; 43, 67, 163 by witnesses of excess 183, 393, 1577). FINITE-EXACT: optimal at
  3, 7, 23; not optimal at 11, 19, 31, 43. OPEN (CONJECTURE, revised): the natural order is optimal only on a thin set of
  primes inside L+, characterised by L(1, chi_(-p)) dominating the twisted values that control the rotation orders;
  whether that set is finite was tested by a pre-registered exact scan (THM-4666 1c, the Farey sweep over all rotation
  orders): the finite-list prediction FAILED, and 46 of the 245 primes with L >= 2.4 in (479, 20000] keep the natural order
  as their best rotation order (set R+), all with L >= 3.58, with a boundary rising from about 3.4 (p ~ 1000) to
  4.35-4.47 (p ~ 16000-20000). R+ looks infinite and thin; optimality among all orders is open even on R+.
session: opus-2026-10-17-S40
source: 05-knowledge/results/binary_information_class_groups_largest_anomalies_20261017.md
scripts:
  - 04-computation/experiments/paley_ranking_local_search_20261017.py (runs '131 30 150', pre-registered '132:500 12 80' and '501:1000 6 40'; + .out, .json)
  - 04-computation/experiments/paley_rotation_orders_20261017.py (rotation scans m <= 30 at 47, 71; m <= 12 for 167-479; m <= 6 for 503-983)
  - 04-computation/experiments/paley_mas_exact_20261017.py (exact MAS, p <= 31); paley_mas_milp_20261017.py (plain 3-dicycle MIP in HiGHS: optimal at p = 23 in 72 s; at p = 31 the 1200 s limit stops it at incumbent 283, bound 307, against MAS 285)
  - 04-computation/experiments/paley_best_orders_vs_affine_20261017.py (structure of the witnesses)
  - audit: 04-computation/experiments/s40_audit_20261017/audit_rotation_witnesses.py, audit_search.py (rotation certificates, 480-restart annealing at 47, 71)
related: [THM-4666, HYP-9286, THM-871, THM-4644, HYP-9277]
---

# HYP-9287 — class numbers and optimal rankings of Paley tournaments

## The original statement and its refutation

For primes `p ≡ 3 (mod 4)` write `L = L(1, χ_(-p)) = π h(-p)/sqrt(p)` (`p > 3`; `L(1, χ_(-3)) = π/(3 sqrt 3) = 0.60`). The
natural order `0 < 1 < ... < p-1` has excess `p h(-p)` (THM-4666). The first form of this hypothesis said: the natural
order is a maximum acyclic ordering iff `L >= τ(p)`, a slowly rising threshold. Necessary conditions are proved: the
natural order must be insertion-locally optimal, i.e. `p ∈ L+` (THM-4666 (5)), a set of relative density 0 (Kalmynin
Thm 1, via Baker-Montgomery).

**Certificates against the threshold form** (rotation orders: list `x ∈ {0, ..., p-1}` by `(x·num) mod den`; every excess
recomputed in pure Python by the audit or in the S40 scan):

| p | L | natural `p h` | rotation order | excess | source |
|---:|---:|---:|---|---:|---|
| 167 | 2.674 | 1837 | `7599/19912` | 1861 (best known 1863) | audit; S40 scan `k/(mp)` = 701/(11·167) |
| 263 | 2.518 | 3419 | `57115/72812` | 4079 (best known 4089) | audit; S40 scan |
| 383 | 2.729 | 6511 | `38609/72644` | 7519 (best known 7529) | audit; S40 scan |
| 439 | 2.249 | 6585 | `27601/171580` | 9313 | audit (local search had 7353) |
| 503 | 2.942 | 10563 | `9521/197064` | 11657 | audit |
| 599 | 3.209 | 14975 | `12767/298528` | 15507 | audit (the S40 scan with m <= 6 missed it) |
| 647 | 2.841 | 14881 | `10901/258330` | 17661 | audit |
| 887 | 3.059 | 25723 | `348751/880400` | 28235 | audit; equals the Farey-sweep maximum (`α = 532/1343`) |
| 911 | 3.227 | 28241 | `358/1263` | 30327 | S40 Farey sweep (exact maximum over rotation orders; not audited) |
| 983 | 2.705 | 26541 | `376249/642500` | 34225 | audit (S40 scan: 32103 at `1727/(6·983)`) |

**Unbeaten by any known order** (rotation family, insertion and window-DP polishing, annealing with 480 restarts at 47
and 71): 3, 7, 23 (exact optima) and 47, 71, 191, 239, 311, 359, 431, 479 (audit, `L >= 1.97`). The coarse S40 scans
(m <= 6) also left 599, 719, 839 unbeaten, but the audit's full Farey-midpoint family beats 599, so coarse scans are not
evidence; 719 and 839 await the full family. A pure threshold in `L` fails: 503 (`L = 2.942`) is beaten while 191
(`L = 2.955`) is not, and the boundary rises quickly with `p` (2.73-2.955 for `p <= 500`, above 3.21 at 599). The local-search optima recorded earlier (best excess
`0.75-0.96 p^(3/2)`) were search-limited; rotation orders reach `0.96-1.04 p^(3/2)` (263: 4089, 383: 7529, 463: 9553,
983: 32103).

## Revised conjecture

The natural order is a maximum acyclic ordering of `P_p` only on a thin set of primes `p ∈ L+`, and among those exactly
when `L(1, χ_(-p))` dominates the twisted values `L(1, χ_(-p) ψ)` (`ψ` of small conductor `m`) that govern the rotation
orders `x -> {kx/(mp)}`. Mechanism (PROVED as an identity, heuristic as a criterion): in multiplier coordinates a rotation
order sorts first by a block label and then by `<cx>`, so its excess is `χ(c)[p h(-p) - 2 Σ χ(t - s)]` over the pairs whose
block labels are inverted (for `m = 2` the blocks are the Eisenstein parity sets `floor(d u/p) mod 2`); these crossing
sums are character analogues of Dedekind sums (Berndt), i.e. combinations of twisted `L`-values. NUMERICAL: the twist hint seen for `p <= 1000` does not generalise. On the 245 primes with `L >= 2.4` in
`(479, 20000]`, a logistic boundary in `(L(1, χ_(-p)), log p)` classifies `R+` with leave-one-out accuracy 0.963. A threshold on
`L` alone reaches 0.935. Adding `L(1, χ_(-8p))/L` or `L(1, χ_(-5p))/L` does not help (0.955, 0.951;
paley_rplus_twist_separation_20261017.py). The exact criterion is the Farey prefix positivity of THM-4666 1c;
its analytic content beyond a slowly rising `L`-threshold is open.

## Structure of the earlier witnesses (paley_best_orders_vs_affine_20261017.py)

At primes where local search beat the natural order, the best order found is either close to an affine order
(normalised Kendall distance 0.11-0.25: 743: 0.114, 103: 0.131, 127: 0.134, 863: 0.166) or generic (0.41-0.46, near the
random level 0.33-0.45 for `p = 43-991`). These are local-search witnesses only; rotation orders are better at many
primes, so this describes the search, not the optimum.

## The rotation family exactly (THM-4666 1c, pre-registered scan to 20000)

The Farey sweep computes the best rotation order in `O(p^2)`.

* **Agreement.** It reproduces every audited certificate, and it improves 911 to 30327 (`α = 358/1263`).
* **Survivors below 480.** `R+ ∩ [3, 479]` is 3, 7, 23, 47, 71, 191, 239, 311, 359, 431, 479; 11, 19, 31, 43 are beaten by
  other orders.
* **Survivors in `(479, 20000]`, among the 245 primes with `L >= 2.4`:** 719, 839, 1151, 1319, 1511, 1559, 2351, 2999, 3191, 3671, 3911, 4079, 4679, 4751, 4871, 5591, 5711, 5879, 6551, 6599, 8111, 8951, 9239, 9719, 10391, 10559, 10631, 10799, 11399, 11471, 12239, 12479, 12911, 13751, 14159, 14759, 14951, 15791, 16631, 17471, 18119, 18191, 18311, 18959, 19031, 19319.

Boundary by range (minimum `L` of the survivors / maximum `L` of the beaten):

| `p` range | minimum `L`, survivors | maximum `L`, beaten |
|---|---:|---:|
| `[480, 1000)` | 3.579 | 3.227 |
| `[1000, 2000)` | 3.797 | 3.424 |
| `[2000, 4000)` | 3.837 | 3.881 |
| `[4000, 8000)` | 4.096 | 4.155 |
| `[8000, 12000)` | 4.149 | 4.224 |
| `[12000, 16000)` | 4.472 | 4.252 |
| `[16000, 20000]` | 4.350 | 4.422 |

The overlapping ranges show that `L(1, χ_(-p))` alone does not decide even rotation-optimality. The exact criterion is
the prefix positivity of the Farey jump sequence.

**Full census to 20000** (all 1136 primes `p ≡ 3 (mod 4)`; `paley_rotation_farey_sweep_3_20000_20261017.out`):

* **`R+` has 57 members:** 3, 7, 23, 47, 71, 191, 239, 311, 359, 431, 479 and the 46 above.
* **Low-`L` members.** The only ones with `L < 2.4` are 3, 7, 23, 47 (`L(1, χ_(-3)) = 0.60`; the script prints 1.81 at
  `p = 3`, the same unit-factor display bug).
* **`R+ ⊂ L+`** at every member (observed, not proved).
* **Per 2000-block:** 17, 5, 7, 2, 4, 6, 4, 4, 2, 6. The count does not thin out, so `R+` looks infinite.
* **A necessary condition.** The first `(p-1)/2` Farey jumps are `-2kχ(k)`, so `R+` requires
  `Σ_(k<=K) kχ(k) >= 0` for `K <= (p-1)/2`. This is necessary but far from sufficient: 118 of the 245 high-`L` primes
  satisfy it and are still beaten.

**Beyond rotation orders on `R+`** (`paley_rplus_anneal_20261017.py`, the audit's polishing and annealing routines):

| p | natural | after insertion + window-DP polishing | annealing best, 48 random restarts |
|---:|---:|---:|---:|
| 719 | 22289 | 22289 | 16025 |
| 839 | 27687 | 27687 | 20327 |
| 1151 | 47191 | 47191 | 32839 |

So the natural order is a local optimum in a larger neighbourhood. Annealing from random starts reaches only 0.70-0.73
of its excess. That is consistent with optimality on `R+`, but only weak evidence: by the calibration rule above,
"unbeaten" is not a certificate.

## Tests

* An exact MAS for `P_47` or `P_71` (the first unbeaten primes beyond the DP range), e.g. the exact DP modulo
  `Aut(P_p)` with meet-in-the-middle of Chindelevitch-Harutyunyan (their App. E certified `p = 43`); a plain MIP already
  leaves a gap of 22 at `p = 31` after 1200 s.
* The finite-list question: run the full rotation family on every `p ≡ 3 (mod 4)` up to 2000 with `L >= 2.9`
  (pre-register first) and record the largest prime where the natural order survives.
* Any future "unbeaten" claim must include the full rotation family (all Farey midpoints of order `p - 1`) and be
  calibrated on known beaten primes (199, 223, 463).
* An explicit formula for the rotation-order excess in terms of twisted `L`-values, and a test of the revised criterion
  on `p <= 2000` (pre-register first).
