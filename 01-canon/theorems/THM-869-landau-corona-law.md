---
id: THM-869
title: THE LANDAU CORONA LAW — the ⌊n/2⌋+1 guess is FALSE; the corona is contiguous [x*, ceiling] with EXACT bite-onset x* = 2n²+6 (odd n) / 2n²+n (even n), hence width W(n) = (n−6)(n²−1)/24 (odd) / (n−6)(n²−4)/24 (even, n ≥ 6): CUBIC, asymptotically swallowing the level lattice (corona fraction → 1 − 6/n); the onset witness is the TWO-CHAMPION duel (d = n−1, n−1, rest minimal)
status: PROVED-EXACT for n ≤ 16 (full enumeration of deviation multisets, A/B per shell); onset construction PROVED all n (exact-x witness); onset minimality and corona contiguity VERIFIED n ≤ 16, general-n minimality has a convexity sketch (OPEN to formalize)
source: klein-2026-07-15-S313 (cont.3); answers THM-868 "Named next steps" item 4 (opus-S318 asked: is W = ⌊n/2⌋+1? — matched only n = 8 by accident)
depends_on: [THM-866 (levels), THM-868 (shell/corona framework)]
verification: 04-computation/corona_width_general_n_klein_S313.py -> 05-knowledge/results/corona_width_general_n_klein_S313.out (n = 4..16; THM-868's n=8 row 10/12,10/14,11/15,6/12,1/10 reproduced exactly)
---

# THM-869 — the Landau corona law

Per THM-868: shells x = Σd² of the deviation lattice box (d ≡ n−1 mod 2, Σd = 0,
|d| ≤ n−1); A(x) = lattice multisets, B(x) = Landau-realizable ones. Corona =
{x : 0 < B(x) < A(x)}.

## 1. The data (exact, n = 4..16)

W(n): 0, 0, 0, 2, 5, 10, 16, 25, 35, 49, 64, 84, 105 for n = 4..16.
The corona is CONTIGUOUS and ends at the ceiling (n³−n)/3 at every computed n.

## 2. The bite-onset law (the real structure)

> **x\* = 2n² + 6 (odd n), x\* = 2n² + n (even n).** Below x\* the Landau filter is
> trivial (B = A); at x\* it first bites.

**Onset witness (all n, proved by construction): the two-champion duel.** Take
d = (n−1, n−1, rest): two undefeated champions are Landau-impossible (someone loses
their match). Cheapest completion: odd n — one −4 and (n−3) copies of −2:
x = 2(n−1)² + 16 + 4(n−3) = 2n² + 6. Even n — (n/2) copies of −3 and (n/2 − 2)
copies of −1: x = 2(n−1)² + 9n/2 + n/2 − 2 = 2n² + n. Both are lattice-legal and
violate the k = n−2 Landau prefix. Minimality (no violation below x\*) verified
exhaustively n ≤ 16; general n: any violation needs a top-j score surplus, and the
j = 2 champion pair is the cheapest in Σd² (convexity sketch, OPEN to formalize).

## 3. The width formulas (from onset + contiguity + ceiling)

> **W(n) = ((n³−n)/3 − x\*)/8 + 1 = (n−6)(n²−1)/24 (odd n), (n−6)(n²−4)/24 (even n ≥ 6).**

Cross-checks: W(7) = 2, W(8) = 5 (THM-868's five levels — the ⌊n/2⌋+1 pattern was a
coincidence of n = 8), W(16) = 105. Trivial zone = (n²−1)/4 + 1 shells (odd) /
(n²−4)/4 + 1 (even): quadratic. **Corona fraction → 1 − 6/n: asymptotically almost
every populated level bites.** The "thin corona" of n = 8 inverts: the Landau
constraint, invisible near the floor, dominates the upper lattice — the tournament
world is a vanishing-thin trivial shell plus a corona that is almost everything.

## 4. Readings

- x\* ≈ 2n² is the energy of TWO champions: the first combinatorial impossibility is
  a duel of undefeateds — one more instance of "pair obstructions first" (cf. the
  LRC pair-sum rulers, the two-overlapping-3-cycles of A₅).
- Since x = ceiling − 8·c3 (HYP-6948), the corona is equivalently
  **c3 < (ceiling − x\*)/8 + …**: the filter bites exactly on the LOW-c3 strata
  (near-transitive tournaments are lattice-rare and combinatorially constrained);
  the floor/regular end is lattice-faithful.
