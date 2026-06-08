---
id: HYP-2360
name: triangle-recursion-temperature-ladder-and-pfaffian-translator
status: MIXED — (1)(4) PROVED+verified; (2)(3) verified; the "temperature ladder" (5) is the conjecture
date: 2026-06-07
session: claudebox-2026-06-07-S715
depends_on:
  - THM-291   # Mode-B n->n+2 multilinear recursion (overlap + boundary wiring + apex)
  - THM-337   # base-path staircase order-3 recurrence H=3H(k-1)+H(k-2)+H(k-3)
  - THM-329   # A038375 max-H new terms
  - THM-326   # staircase HP = independence polynomial (H = I(Omega,2))
  - THM-174   # det(I+2A)=Pf(S)^2
  - THM-435   # Pfaffian = the even/odd seam (S713)
provisional_id: true
---

# HYP-2360: The triangular IE recursion is one recursion at three "temperatures"; the Pfaffian is the translator

Engages the user's 7-piece tiling decomposition (3 corner (n-1)-triangles {A,B,C} +, 3 edge (n-2)-overlaps
{D,E,F} -, 1 center (n-3) {G} +) and the request to use the Pfaffian to translate across domains.

## (1) PROVED: the 7-piece IE is the triangular partition of unity

With `Tri(k) = k(k+1)/2` cells in a side-k triangular grid:
```
   3*Tri(n-1) - 3*Tri(n-2) + Tri(n-3) = Tri(n)      (verified n=3..39, identically)
```
Per-cell net multiplicity is 1 everywhere: a corner cell lies in 1 big triangle; an edge cell in 2 big
minus 1 overlap; an interior cell in 3 big minus 3 overlaps plus 1 center (`3-3+1`). The coefficient
vector `(1,-3,3,-1)` is `(x-1)^3` — the **third finite difference** `Delta^3`. So the user's
construction is the Mobius reconstruction of the staircase triangle by its three corner subtriangles.

## (2) PROVED/verified: RECURSIVE TRUTH A — tile-additive invariants obey (x-1)^3 (are quadratic in n)

If `F(T) = sum over cells of a local value` (tile-additive), then `F(n) - 3F(n-1) + 3F(n-2) - F(n-3) = 0`,
i.e. `Delta^3 F = 0`, i.e. **F is a quadratic polynomial in n**. Verified `Delta^3 = 0` for `#tiles =
C(n-1,2)`, `#arcs = C(n,2)`, `#vertices = n`. This is the exact algebraic content of the user's IE: the
geometry of the triangle dictates that everything *additive over cells* grows quadratically.

## (3) verified: RECURSIVE TRUTH B — H is multiplicative, so it breaks (x-1)^3 but keeps the leading 3

`H = I(Omega, 2)` (THM-326) is a PRODUCT over the conflict graph, not a sum over cells, so it is NOT
tile-additive: `Delta^3` of A038375 `= [-2,8,12,94,214,1896,...] != 0`. The base-path staircase family
(THM-337) instead satisfies, with H `= [1,5,17,57,193,653,2209,7473,25281]` (verified by direct count),
```
   H(k) = 3 H(k-1) + H(k-2) + H(k-3)        char poly x^3 - 3x^2 - x - 1,  lambda ~ 3.383.
```
The **leading coefficient 3 survives** (the three corners), but the lower-order signs flip from the
additive IE's alternating `-3,+1` to **positive** `+1,+1`. Additive geometry -> multiplicative
combinatorics turns subtraction into addition; the backbone (the 3) is invariant, the corrections change
sign.

## (4) PROVED/verified: the PFAFFIAN realizes the user's "negative (n-2) pieces"

The user's instinct that the size-(n-2) pieces `{D,E,F}` are *subtracted* is exactly the Pfaffian's
minor expansion:
```
   Pf(S) = sum_{j} (-1)^{j-1} S_{1j} Pf(S_{1-hat,j-hat})       (verified 2000/2000 random n=6 tournaments)
```
The minors `S_{1-hat,j-hat}` are **(n-2)-vertex deleted subtournaments**, entering with **alternating
sign**. This is Mode-B (THM-291, the n->n-2 step: overlap + boundary) made into an *exact signed
recursion*. So among the three step-sizes, the Pfaffian is the n->n-2 (Mode-B) translator, sitting
between the additive n->n-1 IE and the multiplicative H, and its sign structure is permutation-parity
(neither uniformly alternating nor uniformly positive).

## (5) THE CONJECTURE — the temperature ladder

The same triangular recursion appears at three "temperatures", differing only in the SIGN structure of
the sub-leading (n-2, n-3) corrections, with the leading 3 (three corners) invariant:

| temperature | object | step | recursion / char poly | growth |
|---|---|---|---|---|
| additive (cold) | area / linear coeff | n->n-1 | `(x-1)^3`, signs `+3,-3,+1` | quadratic |
| Pfaffian (mid) | `Pf(S)`, `det=Pf^2` | n->n-2 | signed minor expansion (parity signs) | `sqrt(det)` |
| multiplicative (hot) | `H = I(Omega,2)` | n->n-2 (THM-337) | `x^3-3x^2-x-1`, signs `+3,+1,+1` | exponential `~3.383^k` |

**Conjecture:** these are one recursion analytically continued in a "sign/temperature" parameter; the
leading 3 is the geometric backbone (3 corners), and the additive->Pfaffian->multiplicative passage is
the same additive-vs-multiplicative (flat-vs-peaked, S714; 21-vs-22, S712) axis the project keeps hitting.

## The Pfaffian as cross-domain translator (the creative request)

`Pf(S)` is the single object visible in all five domains, which is why it translates between them:
- **topology**: Euler class — Chern-Gauss-Bonnet `chi = integral Pf(curvature)`; and THM-120 `|Pf|`
  separates the C-phase (beta_1) from the S-phase (beta_3).
- **geometry**: `|Pf| = sqrt(det(MM*))` = oriented volume / the autocorrelation determinant (S714).
- **graphs**: `Pf(S) = signed perfect-matching sum` (FKT), and mod 2 `= #PM(K_n) = (n-1)!!` (odd, THM-435).
- **tournaments**: `Pf(S_T)` odd; `det(I+2A)=Pf^2` (THM-174); the n->n-2 recursion above.
- **algebras**: `Pf` = Berezin integral `integral exp(1/2 theta^T S theta) d-theta` (top exterior form);
  the doubling tower R->C->H->O->S (n=2,3,5,9,17) is the Cayley-Dickson shadow of the even/odd seam.

So the Pfaffian is the Rosetta object: Euler characteristic (topology) = curvature integral (geometry)
= signed matchings (graphs) = skew-tournament parity (tournaments) = Berezin top form (algebras).

## Next
- prove the temperature-ladder as a literal one-parameter family of recursions (sign-twisted transfer);
- find the n->n-1 *multiplicative* recursion (Mode A / "1+2^d") and place it on the ladder;
- recursive lower-bound families for A038375 better than base-path (lambda 3.383 << true growth ~ n/2);
  the maximizers are circulant/Paley (THM-336/338) — is there a circulant-recursive construction?
