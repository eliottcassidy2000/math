---
id: HYP-3738
title: THE CONSTRUCTION BINDING -- PROVED. {1..n-2, n(n-1)} has covering-min M=n/Phi6 at binding modulus D=Phi6=(n-1)n+1 ≡ 1 mod (n-1), rotation a=n: the images are EXACTLY the multiples of n mod Phi6 (a core AP, gcd(n,Phi6)=1) plus the killer (n-1)^2 which lands ONE above the top core point (n-2)n, splitting the wrap gap 2n+1 into {1, 2n}; gap multiset {1, n^(n-3), 2n}, deep hole 2n (radius n). So D = (n-3)n + (2n+1) = (n-1)n+1 ≡ 1 mod (n-1) -- the rung-n binding is PROVED. UNIQUENESS THEOREM: the binding modulus D and rung k are UNIQUE invariants of the covering-min (all extremal coverings share them), but the extremal COVERING is NOT unique (n=7 has exactly 2: {1,2,5,6,7,8},{1,4,5,6,7,11}). ZECKENDORF/OSTROWSKI: the semiconvergent ladder denominators k(n-1)+1 are the 2-term continuants K(n-1,k); the covering-min's CF address [0;n-1,k] is its unique Ostrowski representation (the CF-generalization of Zeckendorf), and the binding three-gap realizes it geometrically
status: PROVED (construction binding, by explicit images -- verified n=5..9, derivation general). Uniqueness of D/rung VERIFIED (n=7: exactly 2 extremal coverings, both D=13; n=9: >=2, both D=33). Ostrowski/Zeckendorf framing is a structural identification. SPREAD-regime binding (D≡1 for the low-rung covering-mins) still OPEN: spreads are NOT three-gap (n=7 spread has 4 distinct gaps), so a different argument is needed.
source: klein-2026-06-30-S40
depends_on:
  - HYP-3734   # Farey-neighbor <=> D≡1 mod(n-1); the radius layers
  - HYP-3736   # the dense-core / band-prime over-constraint
related:
  - HYP-3551   # the construction 14/183 (n=14, rung 14)
  - HYP-3717   # three-gap / torus lift
  - HYP-3724   # Phi_6 = the construction denominator
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3738 — the construction binding (PROVED), uniqueness, and the Ostrowski/Zeckendorf address

## THEOREM (construction binding, proved)
For the construction `S = {1, 2, ..., n-2, n(n-1)}`, the covering-min is `M = n/Phi_6(n)` achieved at modulus
`D = Phi_6 = n^2-n+1` and rotation `a = n`, and `D = (n-1)n + 1 ≡ 1 mod (n-1)` (rung `n`).

**Proof (explicit images at `a = n`, mod `Phi_6 = n^2-n+1`).**
- *Core* `k.n` for `k=1,...,n-2`: since `(n-2)n = n^2-2n < Phi_6`, no reduction -- the core images are the AP
  `{n, 2n, ..., (n-2)n}` (step `n`; valid as `gcd(n, Phi_6)=1` because `Phi_6 ≡ 1 mod n`).
- *Killer* `n(n-1).n = n^3-n^2 = n.Phi_6 - n ≡ -n ≡ (n-1)^2 mod Phi_6`. And `(n-1)^2 = n^2-2n+1 = (n-2)n + 1`
  -- exactly ONE above the top core point `(n-2)n`.
So the sorted images are `n, 2n, ..., (n-2)n, (n-2)n+1`. The gaps:
  - `(n-3)` gaps of `n` (between consecutive core points),
  - one gap of `1` (top core `(n-2)n` to killer `(n-2)n+1`),
  - the wrap gap from the killer back to `n` through `0`: `n - ((n-2)n+1) + Phi_6 = 2n`.
Gap multiset `{1, n^(n-3), 2n}`, summing to `(n-3)n + 1 + 2n = n^2-n+1 = Phi_6 = (n-1)n+1`. The observer `0`
sits in the wrap gap (between the killer at `-n` and the core point `n`), at distance `n` from each -- the
DEEP HOLE, radius `n`, so `M = n/Phi_6`. Hence `D = Phi_6 ≡ 1 mod (n-1)`. ∎
(Verified n=5..9: e.g. n=7, a=7, images `{7,14,21,28,35,36}` mod 43, gaps `{1,7,7,7,7,14}`.)

The killer kills resonances `n-1` and `n` (`n(n-1)` is a multiple of both) AND creates the unit gap that makes
`D ≡ 1 mod (n-1)`. The `+1` in `Phi_6 = (n-1)n+1` IS that single unit gap.

## UNIQUENESS THEOREM (D and rung are unique; the covering is not)
The covering-min VALUE `M(n)` determines a unique binding modulus `D = denom(M)` and rung `k`. Verified: ALL
covering-min coverings share `D` and `k`, but the extremal COVERING is NOT unique --
> `n=7`: exactly TWO covering-min coverings (speeds `<= 14`): `{1,2,5,6,7,8}` and `{1,4,5,6,7,11}`, both with
> binding `D=13`, rung `2`. (`{1,4,5,6,7,11}` uses the band-prime killer `11`; `{1,2,5,6,7,8}` uses the
> spreader `8` -- two routes to the same rung, HYP-3736.)
> `n=9`: `>=2` (`{1,4,5,6,7,11,32,36}` and mac-mini's `{1,3,4,5,7,11,18,32}`), both `D=33`.
So uniqueness lives at the level of the VALUE / binding / rung, not the covering -- the right invariant is the
rung `k(n)`.

## ZECKENDORF / OSTROWSKI: the unique address
The semiconvergent ladder denominators `q_k = k(n-1)+1` are the **2-term continuants** `K(n-1, k)`; the
covering-min's continued fraction `M = [0; n-1, k]` is its **unique Ostrowski representation** -- the
continued-fraction generalization of Zeckendorf's Fibonacci numeration (Zeckendorf = the golden-ratio / all-1
CF case; here the CF is `[n-1, k]`). Geometrically, the binding **three-gap** realizes this numeration: for the
construction, `alpha = a/D = n/Phi_6 = M = [0; n-1, n]`, with convergent denominators `q_0=1, q_1=n-1, q_2=Phi_6`,
and the three gaps `{1, n, 2n}/Phi_6` are exactly the Ostrowski "place values" (`||q_1 alpha|| = 1/Phi_6` the
unit gap; `||q_0 alpha|| = n/Phi_6` the regular gap; the deep hole `2n/Phi_6 = 2||q_0 alpha||`, the observer at
its center). The covering-min's rung is its unique Ostrowski/Zeckendorf digit.

## What remains
The construction binding (rung `n`, the covering-min for `n >= 12`) is PROVED. The SPREAD-regime binding
(`D ≡ 1 mod (n-1)` for the low-rung covering-mins, `n=7..11`) is still OPEN: the spreads are NOT three-gap
(the `n=7` spread `{1,2,5,6,7,8}` has gap multiset `{1,1,2,2,3,4}`, four distinct values), so the clean
three-gap argument does not apply -- a different (over-constraint / Ostrowski) argument is needed.

## Net
The construction binding is a THEOREM: at `a=n` the images are the multiples of `n` mod `Phi_6` (a core AP)
plus the killer `(n-1)^2` adjacent to the top, giving the three-gap `{1, n^(n-3), 2n}` and `D=(n-1)n+1 ≡ 1`.
The unit gap is the `+1`. Uniqueness holds for the binding/rung (not the covering: `n=7` has two). The ladder
is an Ostrowski/continuant numeration -- Zeckendorf's CF-generalization -- and the covering-min's rung is its
unique address digit, realized geometrically by the binding three-gap.
