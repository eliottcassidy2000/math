---
id: HYP-3736
title: The OVER-CONSTRAINT that sets k_min -- the radius-1 band-prime KILLER-OR-TRANSVERSAL mechanism + the DENSE-CORE TRANSVERSAL LEMMA. PROVED LEMMA: {1..m} is a +-transversal mod EVERY prime p<=2m+1 (each unit pair {u,p-u} has min <= (p-1)/2 <= m). Consequence: a covering-min covering must, at each radius-1 band prime p in (n, 2/M), EITHER include a multiple of p (a KILLER -- a large speed) OR be a +-transversal mod p; the dense core {1..n-2} supplies FREE transversals for all band primes p<=2n-3, so the construction (dense core + 1 big killer for resonances n-1,n) reaches rung n; spreads trade core speeds for spreaders+band-prime-killers to reach a LOWER rung -- the budget tradeoff (n-1 = core + killers + spreaders) sets k_min. The large speeds in the covering-min ARE the band-prime killers (why small-Smax searches miss it)
status: VERIFIED/PROVED in part. LEMMA proved + verified (m=3..11). Killer-or-transversal mechanism verified on all covering-min spreads n=7,8,9,11 (each band prime is killed-by-inclusion or core-transversal). Construction band-prime handling verified n=7..12 (dense core = transversal mod all its band primes). The exact k_min remains the arithmetic of the budget tradeoff; D≡1 mod(n-1) (covering-min on the ladder) still OPEN (verified n=7..14, HYP-3734).
source: klein-2026-06-30-S39
depends_on:
  - HYP-3734   # the Farey-neighbor reduction + the radius-layer over-constraint
related:
  - HYP-3732   # the Farey-rung invariant
  - THM-523    # resonance-killing (the radius-0 layer)
  - HYP-3731   # the danger-circulant IP
results:
  - 04-computation/farey_rung_spread_family_klein.py
  - 05-knowledge/results/farey_rung_spread_family_klein.out
---

# HYP-3736 — the over-constraint that sets k_min: band-prime killer/transversal + the dense-core lemma

## The dense-core transversal LEMMA (proved)
> **Lemma.** For any prime `p <= 2m+1`, the set `{1,2,...,m}` contains a representative of every unit pair
> `{u, -u}` of `(Z/p)*` (it is a `+-`transversal). *Proof.* Each pair `{u, p-u}` has smaller element
> `min(u,p-u) <= (p-1)/2 <= m`, which lies in `{1,...,m}`. ∎

Equivalently: `{1,...,m}` is a radius-1 covering of `Z/p` for every prime `p <= 2m+1` (every rotation `a` has
some `s<=m` with `s.a in {1,-1} mod p`).

## The over-constraint at the radius-1 band (leveraging it)
For a covering with gap `M`, every modulus `D'` demands covering radius `floor(M.D')`. The **radius-1 band** is
`D' in (n-1, 2/M)`. At a band PRIME `p` (with no multiple of `p` forced, since `p>n`), radius-1 means a
`+-`transversal mod `p` -- which needs `>= (p-1)/2` speeds. So at each band prime the covering must
**EITHER include a multiple of `p` (a KILLER) OR be a `+-`transversal mod `p`.** Verified on every covering-min
spread:
```
n=7  M=2/13 : band primes {11}        -> 11: core-transversal
n=8  M=2/15 : band primes {11,13}     -> 11: KILLER(11),  13: transversal
n=9  M=4/33 : band primes {11,13}     -> 11: KILLER(11),  13: transversal
n=11 M=3/31 : band primes {13,17,19}  -> 13,17,19: ALL KILLERS (13,17,19 in S)
```
**The large speeds in the covering-min ARE the band-prime killers** (`11` at `n=8,9`; `13,17,19` at `n=11`).
This is exactly why a small-`Smax` search misses the covering-min: the killers are large (`>= p > n`).

## Why the construction reaches rung n, and how spreads beat it
The construction `{1,...,n-2, n(n-1)}` has dense core `{1,...,n-2}` = `+-`transversal mod EVERY prime
`p <= 2(n-2)+1 = 2n-3` (lemma, `m=n-2`). Its radius-1 band primes are all `<= 2/M < 2n-1`, hence `<= 2n-3`, so
the **dense core handles every band prime for FREE**; the one big speed `n(n-1)` kills resonances `n-1, n`.
That is the whole covering -- rung `n`, no separate band-prime killers needed (verified `n=7..12`).

To reach a **lower rung** a spread must place its speeds into the thinner annulus of the lower-rung binding
`D = k(n-1)+1` (radius `k < n`), which forces giving up dense-core speeds for spreaders. Removing core speeds
breaks the free transversal coverage, so the spread must re-buy the now-uncovered band primes with KILLERS
(large speeds) or partial transversals. The **budget identity**
> `n-1 = (resonance-killers 2..n) + (band-prime killers/transversals) + (spreaders)`
is the tradeoff that sets `k_min`: a lower rung needs more spreaders (thinner annulus) AND the band primes
still must be paid for; when the budget can't cover both, that rung is infeasible. As `n` grows the band
`(n, 2(n-1))` holds more primes (`~ n/ln n`), the tradeoff tightens, `k_min` rises (`2,2,4,4,3,...`), and by
`n=12` only the construction (core does all the band work) survives.

## Status of the two targets
- **k_min characterization (this HYP):** the mechanism is identified and verified -- band-prime
  killer/transversal, the dense-core lemma giving free transversals up to `p<=2n-3`, the budget tradeoff. The
  exact `k_min(n)` is the optimum of that tradeoff (genuinely arithmetic; `2,2,4,4,3` for `n=7..11`).
- **`D ≡ 1 mod (n-1)` (covering-min on the ladder):** still OPEN. The over-constraint supports it -- the
  binding `D=k(n-1)+1` and the dense-core structure both carry the factor `n-1` -- but no proof yet. Verified
  `n=7..14` (HYP-3734).

## Net
The over-constraint that sets `k_min` is the radius-1 band: every band prime `p in (n,2/M)` must be paid for
by a KILLER (a multiple of `p`, i.e. a large speed) or a `+-`transversal mod `p`. The proved dense-core lemma
(`{1..m}` is a transversal mod every prime `<= 2m+1`) shows the construction gets all band primes free from its
core (rung `n`), while spreads must spend large speeds on band-prime killers to thin the annulus to a lower
rung. The covering-min's large speeds are precisely these band-prime killers -- the concrete reason it sits at
an irregular rung and eludes bounded search.
