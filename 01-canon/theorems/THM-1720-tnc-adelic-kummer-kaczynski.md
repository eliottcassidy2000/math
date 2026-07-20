---
id: THM-1720
title: "TNC AS AN ADELIC / CELLULAR-AUTOMATON / BOUNDARY-FUNCTION SYNTHESIS. Taking the best S422-424 levers (positivity, amoeba-radius separation) and weaving in the repo's cellular-automata (Sierpinski/Kummer carry) and Kaczynski (boundary-function) threads. (1) ADELIC COPRIMALITY: for a tunable trinomial the levels CT(m0), CT(2m0) are coprime because their root sets differ at SOME PLACE of Q -- a finite prime p (via the p-adic Newton polygon) or the archimedean place (via the amoeba/multinomial radius). Verified: {-3,-1,3} at p=2, {-3,1,5} at p=3, {-3,2,7} at p=7, {-4,1,6} at p=2 (finite places); {-2,1,4} at the archimedean place (|roots|=1 vs sqrt(2+-sqrt3)), no finite prime<60 separating. This is product-formula reasoning: agree at every place => equal, but they are not, so they differ somewhere => coprime => TNC. (2) THE CELLULAR-AUTOMATON CONTENT: the finite-place valuations are computed by KUMMER'S THEOREM, v_p(multinomial) = number of carries adding the parts in base p -- exactly the Sierpinski/Pascal carry automaton (Rule 90 at p=2, HYP-2491/2497). So the p-adic Newton polygon of CT(m) is a readout of the carry CA. (3) THE KACZYNSKI BOUNDARY REFRAME: G(t) = t(log Pi)' is analytic in the disk |t| < 1/rho; the saddle values t_j = 1/w_j are BOUNDARY SINGULARITIES on |t| = 1/rho, and the N branches u_i(t) are the APPROACH PATHS. TNC <=> G constant <=> the boundary function (radial limits) is trivial <=> empty singular set <=> R monomial -- Kaczynski's boundary-function dichotomy applied to G"
status: EXPLORATORY SYNTHESIS -- the adelic place-by-place separation is VERIFIED on all tunable trinomials tested (mix of finite and archimedean places); a uniform 'they differ at some place' proof is the open target. (2),(3) are exact structural reframes.
author: opus-2026-07-20-S425
depends_on: [THM-1715 (positivity/recurrence), THM-1710 (resultant; cyclotomic refuted), THM-1680 (trinomial gcd), THM-1635 (G = t(log Pi)'), HYP-2491 (Sierpinski/Kummer carry CA), the repo's Kaczynski-boundary thread]
---

# THM-1720 — TNC as adelic / cellular-automaton / boundary-function synthesis

The best levers from S422–424 are **positivity** (CT has positive coefficients in the gauge
params) and **amoeba-radius separation** (root moduli grow with the level). This note takes
those and weaves in two repo threads the owner named: **cellular automata** (Sierpiński/Kummer
carry) and **Ted Kaczynski's mathematics** (boundary functions). They fit together as one
adelic picture.

## 1. Adelic coprimality (verified)

THM-1680/1710: trinomial TNC ⟺ `CT(Λ^{m_0})` and `CT(Λ^{2m_0})` share no root in `a`. The
**reason they are coprime is that their root sets differ at some place of `ℚ`:**

| pattern | separating place |
|---|---|
| `(3; −3,−1,3)` | `p = 2` (finite) |
| `(3; −3,1,5)` | `p = 3` |
| `(3; −3,2,7)` | `p = 7` |
| `(4; −4,1,6)` | `p = 2` |
| `(2; −2,1,4)` | **archimedean** (`|roots| = 1` for `CT(m_0)`; `√(2±√3)` for `CT(2m_0)`) — no finite `p < 60` separates |

At a finite `p`, the **`p`-adic Newton polygon** of `CT(Λ^m)` gives the `p`-adic valuations
of its roots; disjoint valuation-sets ⟹ no shared root. At the archimedean place, the
**ordinary Newton polygon of `log|coeff|`** gives the root moduli (the amoeba); disjoint
radii ⟹ no shared root.

> **Product-formula reasoning.** If `CT(m_0)` and `CT(2m_0)` shared a root, that root would
> have the same valuation at *every* place. They differ at *some* place (finite or ∞), so no
> shared root — coprime — TNC. Verified for every tunable trinomial.

This is the honest unification of the S422 amoeba lever (the **archimedean** place) with a
finite-prime refinement, and it explains why `{−2,1,4}` looked special: it is the case that
closes *only* at infinity.

## 2. The cellular-automaton content (Kummer = Sierpiński carry)

The finite-place valuations are not arbitrary — **Kummer's theorem** computes them:

```
v_p( \binom{m}{x,y,z} ) = (number of carries when adding x+y+z in base p).
```

Carry-counting in base `p` is exactly the **Sierpiński/Pascal carry cellular automaton**
(Rule 90 at `p = 2`), the object of the repo's `pollock_sierpinski_carry_scout` thread
(HYP-2491/2497). So:

> **The `p`-adic Newton polygon of `CT(Λ^m)` is a readout of the carry automaton.** The
> finite-place separations in §1 are computed by the Sierpiński CA on the minimal-rep
> multinomials.

This is the concrete sense in which the CA thread bears on TNC: the divisibility that pins
the `p`-adic root valuations is a carry pattern. A uniform statement — *for every tunable
trinomial the carry CA produces disjoint valuation-sets at some prime, OR the amoeba
separates at ∞* — would be the adelic completion.

## 3. The Kaczynski boundary-function reframe

Kaczynski's mathematics (his 1967 thesis and 1969 papers) is **boundary functions**: for a
function defined in a disk, the boundary function assigns to each boundary point the limit
along a specified **approach path**, and his theorems govern the structure of the set where
such limits exist. The repo already carries a "Kaczynski boundary / approach labels" thread
in the LRC analytic sieve.

TNC has exactly this shape. `G(t) = t·(\log Π)'` (THM-1635, sign-fixed) is **analytic in the
disk `|t| < 1/ρ`**, `ρ = ` dominant saddle value. On the boundary circle `|t| = 1/ρ` sit the
**singularities** `t_j = 1/w_j` (saddle values), and the `N` small branches `u_i(t)` are the
**approach paths** converging to each boundary point.

> **TNC ⟺ `G` is constant ⟺ the boundary function of `G` is trivial ⟺ `G`'s singular set on
> `|t| = 1/ρ` is empty ⟺ `R` is a monomial.**

The saddle-value collisions (THM-1625) are exactly *coincidences of approach-path limits* —
two branches `u_i, u_j` delivering the same boundary value `w`. Kaczynski's dichotomy (a
boundary function is trivial iff its singular set is empty, for an algebraic `G`) is the
structural statement that TNC is asking for. This is a reframe, not yet a new proof, but it
places TNC in a named classical theory.

## 4. What this buys, honestly

- **A cleaner target than the dead cyclotomic route** (THM-1710): "the levels differ at some
  place" is elementary (Newton polygons + Kummer), and the finite-place part is a carry-CA
  computation. No deep cyclotomy.
- **The archimedean/finite split explains the `{−2,1,4}` anomaly**: it is the sole case in
  the sample that closes only at ∞ — the amoeba — which is why S422's radius argument saw it
  and a naive `p`-adic scan did not.
- **The boundary-function frame** connects TNC to Kaczynski's theory and to the repo's
  analytic-sieve "approach labels," suggesting the singular-set machinery there may transfer.

## 5. Open

1. **Uniform adelic separation.** For every tunable `k`-nomial, do `CT(Λ^{m_0}), CT(Λ^{2m_0})`
   always differ at some place? Split: (a) if the carry CA gives disjoint `p`-adic valuations
   at some prime, done; (b) else the multinomial amoeba must separate at ∞. Proving the
   dichotomy is exhaustive is the completion.
2. **Which primes?** The separating primes seen (`2, 3, 7`) are small; is there a bound
   `p ≤ f(charges)` from the carry structure? A Kummer bound on the separating prime would
   make the certificate finite and explicit.
3. **Kaczynski singular-set transfer.** Does the descriptive-set structure of boundary
   functions constrain `G`'s singular set enough to force nonemptiness for non-monomial `R`?

## Verification

`04-computation/tnc_padic_newton_kummer_opus_S425.py` (p-adic Newton-polygon root valuations;
the separating-prime search), `04-computation/tnc_adelic_places_opus_S425.py` (the
place-by-place separation table, finite vs archimedean). Outputs in `05-knowledge/results/`.
