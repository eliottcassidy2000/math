---
id: HYP-9362
title: "(a) One-flip conjecture: for every odd n the circulant tournament maximising the number H of Hamiltonian paths is the rotational tournament R_n or differs from it in one difference class. (b) Paley peak: for every prime q = 3 mod 4 the Paley tournament is a strict local maximum of H on the circulant class graph, H(P_q) > H(P_q with one difference class reversed). (c) The rotational/Paley ratio H(R_q)/H(P_q) increases to a finite limit greater than 1. (d) Dictionary to be explained or refuted: the forbidden Hamiltonian-path counts {7, 21} = {|PG(2,2)|, |PG(2,4)|} are exactly the orders of the projective planes that extend to 3-designs, whose extensions S(3,4,8) and S(3,6,22) < S(5,8,24) are the combinatorial skeletons of E_8 and the Leech lattice; is there a tournament mechanism (for example through the OCF independence polynomial) behind this alignment?"
status: >
  OPEN.
  FINITE-EXACT evidence:
  - (a): every odd n <= 23 (complete multiplier-class tables; THM-4751 (3)).
  - (b): q = 7, 11, 19, 23 (tables) and q = 31 (one H computation). All one-flip neighbours of P_q are isomorphic, so
    (b) is a single inequality per q.
  - (c): ratios 1.0098, 1.0159, 1.0205, 1.0239 at q = 19, 23, 27, 31 (HYP-9357).
  - (d): the classical facts on both sides are cited in the results note. No map is known, and the coincidence involves
    two values; it is recorded as DICTIONARY.
session: mac-mini-2026-10-21-fanopsp
source: 01-canon/theorems/THM-4751-circulant-h-landscape-paley-peak-at-class-number-distance.md
related:
  - 05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-and-the-kantor-galois-exceptional-primes-3-7-11.md
  - 05-knowledge/results/lattices_codes_tournaments_20261021.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
---

# HYP-9362 — the one-flip conjecture, the Paley peak, and the {7, 21} → {8, 24} dictionary

## (a) One-flip conjecture

On odd `n ≤ 23` the `H`-maximal circulant is:

* `R_n` itself for `n ∈ {3, 5, 13, 15, 17, 19, 21, 23}`;
* a one-flip neighbour of `R_n` for `n ∈ {7, 9, 11}` (`P_7`, `Z_9{2,3,4,8}`, `P_11`).

**Consequence.** By THM-4751 (2), Paley is at flip distance `B_q ≥ 2` from `R_q` for `q ≥ 19`. So (a) would settle the circulant form of HYP-9357's "only if": `P_q` is never the circulant maximum for `q ≥ 19`.

**Test.** Full circulant maxima at `n = 25..31` are expensive (hundreds of classes, minutes each). A cheaper partial test is whether `R_n` is a strict local maximum for `n = 25, 27, 29, 31`; that takes about `n/2` computations each.

## (b) Paley peak

`P_q` is a strict local maximum at `q = 7, 11, 19, 23, 31`. The margins `H(P_q^(1))/H(P_q)` (`paley_flip_margin.out`) are:

| `q` | margin |
|---|---|
| 7 | 0.92593 |
| 11 | 0.97825 |
| 19 | 0.99246 |
| 23 | 0.99542 |
| 31 | 0.99822 |

`1 − margin` decays roughly like `2.5/q²` over this range.

Heuristic: the Paley tournament's doubly regular symmetry makes its subtournament statistics uniform. Breaking one difference class costs `q` arcs of that uniformity at once.

## (c) Rotational/Paley ratio

The ratio grows slowly with `q`: 1.0098, 1.0159, 1.0205, 1.0239. Wormald's blow-up constant for the asymptotic maximum (about 2.8559) and the Paley normalised value `H·2^(q−1)/q!` (2.55, 2.60, 2.63, 2.66) suggest a finite limit. No prediction is made.

## (d) The {7, 21} → {8, 24} dictionary

**Facts (classical; citations in the results note):**

* Hughes (1965): a projective plane of order `n` extends to a 3-design only if `n ∈ {2, 4, 10}`. There is no plane of order 10 (Lam–Thiel–Swiercz 1989).
* The extensions are `PG(2,2)` (7 points) → `S(3,4,8)` and `PG(2,4)` (21 points) → `S(3,6,22) ⊂ S(4,7,23) ⊂ S(5,8,24)`.
* `S(3,4,8)` is the set of 14 weight-4 words of the `[8,4,4]` code that gives `E_8`. The octads of `S(5,8,24)` are the 759 weight-8 words of the Golay code that gives the Leech lattice.
* THM-1370: `7` and `21` are the only values never taken by `H`, at any `n`.
* Two Paley routes lead to `{8, 24}` as well: `P_7` and `P_23` give the codes of length 8 and 24 (THM-4751 (5)).

**The question.** Is there a structural reason that the two forbidden `H` values are the point counts of the only extendable projective planes?

* Both are values of `Φ_3`: `7 = Φ_3(2)` and `21 = Φ_3(4)`. A Φ₃ reading was already examined in HYP-1317. There `Φ_3(6) = 43` is attained, and `q = 6` is not a prime power. The extendable-plane reading picks out exactly `q = 2, 4`.
* Repo history for the null hypothesis:
  * `H = 13 = Φ_3(3)` and `31 = Φ_3(5)` are attained.
  * 31 is rare in the `H` spectrum (S50/S51).
  * 31 is also `|PG(4,2)|` (Goormaghtigh: `31 = 11111_2 = 111_5`).

**What would refute (d) as a meaningful alignment:** an independent characterization of the forbidden values, for example through the OCF independence polynomial, that makes no reference to the projective plane structure. THM-1370's proof is of this kind already, so (d) is at best a dictionary unless a map is found.
