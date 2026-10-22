---
id: HYP-9362
title: "(a) One-flip conjecture: for every odd n, every circulant tournament maximising the number H of Hamiltonian paths is the rotational tournament R_n or differs from it in one difference class; sharper, for odd n >= 13 the H-maximal circulant is R_n itself, the unique maximiser of the additive energy (refines HYP-480). (b) Paley peak: for every prime q = 3 mod 4 the Paley tournament is a strict local maximum of H on the circulant class graph; more generally, for odd n >= 17 every energy minimiser is. (c) The rotational/Paley ratio H(R_q)/H(P_q) increases to a finite limit greater than 1. (d) NUMEROLOGY to be explained or refuted: the H-values {7, 21}, proved never to occur, are the point counts |PG(2,2)|, |PG(2,4)| of the only projective planes that extend to 3-designs, whose extensions S(3,4,8) and S(3,6,22) < S(5,8,24) are the combinatorial skeletons of E_8 and the Leech lattice. (e) Packing law: for odd n >= 17 the energy minimiser maximises alpha_1 and alpha_2 and R_n maximises every alpha_k, 3 <= k <= floor(n/3), in the OCF H = sum 2^k alpha_k."
status: >
  OPEN. Audited in two passes 2026-10-22 (05-knowledge/results/lattice_bridges_audit_20261021.md,
  lattice_bridges_audit_pass2_20261022.md): every number confirmed; statuses and prior art corrected here.
  FINITE-EXACT evidence:
  - (a): complete multiplier-class tables for every odd n <= 25 (THM-4751 (3); n = 25 new, 208 classes; reproduced by the
    audit). R_n is moreover a strict local maximum at n = 27, 29, 31 against all its flip-neighbours (rot*_local_max.out).
  - (b): Paley at q = 7, 11, 19, 23 (tables) and q = 31 (two independent H computations). All one-flip neighbours of P_q
    are isomorphic, so (b) is a single inequality per q. Energy minimisers: strict local maxima for 17 <= n <= 25.
  - (c): ratios 1.0098, 1.0159, 1.0239 at q = 19, 23, 31 (and 1.0205 at the prime power 27, where P_27 is not a circulant).
  - (d): the classical facts on both sides are cited below. No map is known; recorded as NUMEROLOGY (MISTAKE-559 kept the
    NUMEROLOGY verdict on {7, 21}).
  - (e): n = 17, 19, 21 over all classes (alpha_all.out, alpha_all_21.out); the split H(R_n) - H(energy minimum) for
    n = 9..21 (THM-4751 (7)); independently verified by the second audit pass.
session: mac-mini-2026-10-21-fanopsp
source: 01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md
related:
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
  - 05-knowledge/results/lattices_codes_tournaments_20261021.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
  - 01-canon/theorems/THM-137-paley-crossover-mechanism.md
  - 01-canon/theorems/THM-140-energy-H-correlation-reversal.md
---

# HYP-9362 — the one-flip conjecture, the Paley peak, the packing law, and the {7, 21} → {8, 24} numerology

## (a) One-flip conjecture

**Prior art.** This refines HYP-480 ("interval universal `H`-maximiser", exhaustive for primes `p ≤ 19` and `p = 23`, March 2026) and THM-137's open question 2. It extends them to all odd `n` and adds "within one flip".

**Data.** On odd `n ≤ 25` the `H`-maximal circulant is unique. It is:
* `R_n` itself for `n ∈ {3, 5, 13, 15, 17, 19, 21, 23, 25}`;
* a one-flip neighbour of `R_n` for `n ∈ {7, 9, 11}` (`P_7`, `Z_9{2,3,4,8}`, `P_11`).

These are the two extremes of the additive energy `E(S)` (THM-4751 (6)): the energy minimum wins for `n ≤ 11`, and `R_n`, the energy maximum, for `n ≥ 13`.

**Consequence.** By THM-4751 (2), Paley is at flip distance `B_q ≥ 2` from `R_q` for `q ≥ 19`. So (a), in the form "every maximiser is within one flip of `R_n`", would settle the circulant form of HYP-9357's "only if": `P_q` is never the circulant maximum for `q ≥ 19`.

**Pre-registered partial test (done).** `R_n` is a strict local maximum for `n = 25, 27, 29, 31`.
* The best neighbour always flips `{m, m + 1}`.
* Its ratio is 0.99585, 0.99613, 0.99641, 0.99667, so `1 − ratio ≈ 0.104/n`.

## (b) Paley peak

`P_q` is a strict local maximum at `q = 7, 11, 19, 23, 31`. The margins `H(P_q^(1))/H(P_q)` (`paley_flip_margin.out`) are:

| `q` | margin |
|---|---|
| 7 | 0.92593 |
| 11 | 0.97825 |
| 19 | 0.99246 |
| 23 | 0.99542 |
| 31 | 0.99822 |

`(1 − margin)·q²` is 3.63, 2.63, 2.72, 2.42, 1.72 at `q = 7, 11, 19, 23, 31`. So the decay is at least quadratic, with no stable constant. The rotational peak's margin decays only like `1/n` ((a)).

**Generalization.** For `17 ≤ n ≤ 25`, every energy minimiser is a strict local maximum of `H`; at `n = 19, 23` this is Paley. The lists use the hypercube convention: a flip that lands in the same multiplier class counts as a tie (THM-4751 (3)).

**Heuristic.** The doubly regular symmetry of the Paley tournament makes its subtournament statistics uniform. Breaking one difference class costs `q` arcs of that uniformity at once.

## (c) Rotational/Paley ratio

**Data.** The ratio grows slowly with `q`: 1.0098, 1.0159, 1.0239 at `q = 19, 23, 31`. (At `q = 27` it is 1.0205, but 27 is a prime power and `P_27` is not a circulant.)

**Normalised values `H·2^(q−1)/q!`.**
* Rotational: 2.552, 2.598, 2.632, 2.658.
* Paley: 2.527, 2.557, 2.579, 2.596.
* Wormald's blow-up constant, about 2.8559, is the best known lower-bound constant for max `H`.

**Conflict with earlier conjectures.** HYP-9274 (A) conjectures that the Paley value tends to `e`, and HYP-9028 that the rotational value does too. Together they predict ratio → 1. So (c) amounts to the rotational limit exceeding `e`, i.e. to HYP-9028 failing for `R_n`.

The data lean that way: the rotational values rise by 0.0138, then 0.0122, at `n = 27, 29, 31`. No prediction of the limit is made.

## (d) The {7, 21} → {8, 24} numerology

**Facts (classical; citations in the results note):**

* Hughes (1965): a projective plane of order `n` extends to a 3-design only if `n ∈ {2, 4, 10}`. There is no plane of order 10 (Lam–Thiel–Swiercz 1989).
* The extensions are `PG(2,2)` (7 points) → `S(3,4,8)` (one step) and `PG(2,4)` (21 points) → `S(3,6,22) ⊂ S(4,7,23) ⊂ S(5,8,24)` (three steps). There is no `S(4,5,9)`.
* `S(3,4,8)` is the set of 14 weight-4 words of the `[8,4,4]` code that gives `E_8`. The octads of `S(5,8,24)` are the 759 weight-8 words of the Golay code from which the Leech lattice is built. Both chains were rederived by the audit from the codes of `P_7` and `P_23`.
* THM-1370 proves that 7 and 21 are never `H` values, at any `n`. That they are the only odd values missed is THM-1370's completeness CONJECTURE: FINITE-EXACT below 80,265 (THM-4749); MISTAKE-559 (1).
* Two Paley routes lead to `{8, 24}` as well: `P_7` and `P_23` give the codes of length 8 and 24 (THM-481).

**Prior art.** Earlier chains `7 → 21 → 24`:
* opus S80, `fano_hamming_forbidden.out`;
* opus S88, `golay_baer_bridge_88.py`;
* kind-pasteur S100, `projective_geometry_deep.out`: "PG(2,4) + 3 = 24";
* the 2026-05-02 reflection `sphere-packing-eight-twentyfour.md`;
* concurrently, opus S52, `lattices_modular_repunits_20261021.md` §5, also typed NUMEROLOGY.

New here is only the use of Hughes and Lam as the selection principle.

**The question.** Is there a structural reason that the two forbidden `H` values are the point counts of the only extendable projective planes?

* Both are values of `Φ_3`: `7 = Φ_3(2)` and `21 = Φ_3(4)`. A Φ₃ reading was already examined in HYP-1317. There `Φ_3(6) = 43` is attained, and 6 is not a prime power. The extendable-plane reading picks out exactly `q = 2, 4`.
* Repo history for the null hypothesis:
  * `H = 13 = Φ_3(3)` is attained (first at `n = 5`), and so are `31 = Φ_3(5)` and `43 = Φ_3(6)` (first at `n = 6`; audit witnesses).
  * 31 is rare in the `H` spectrum (S50/S51).
  * 31 is also `|PG(4,2)|` (Goormaghtigh: `31 = 11111_2 = 111_5`).

**What would refute (d) as a meaningful alignment.** An independent characterization of the forbidden values that makes no reference to the projective plane structure. THM-1370's proof is of this kind already: 7 and 21 sit in the strong-floor crevasses `(5, 9)` and `(15, 25)`. So by this criterion (d) is NUMEROLOGY unless a map is found.

## (e) Packing law

**Data.** In `H = Σ_k 2^k α_k`, over all classes at `n = 17, 19, 21` (THM-4751 (7)):
* an energy minimiser is the unique maximiser of `α_1` and `α_2`;
* `R_n` is the unique maximiser of every `α_k` with `3 ≤ k ≤ ⌊n/3⌋`, except the perfect triangle factors `α_7` at `n = 21`.

At `n = 13, 15`, `R_n` maximises `α_3` and `α_4` (and at `n = 15` the triangle factors `α_5` are maximised elsewhere). The `H`-switch at `n = 13` comes through `α_4`, when four disjoint triangles first fit.

**Conjecture.** The pattern of `n = 17, 19, 21` holds for every odd `n ≥ 17`, with `α_(n/3)` excepted when `3 | n`. It would imply (a) in the sharper form as soon as the packing terms dominate, which they do from `n = 13` on in all the data.

**Tests.**
* `n = 23, 25` need a lower-memory version of `alpha_decomp.c`.
