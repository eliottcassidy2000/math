---
id: HYP-3741
title: THE WITNESS HIERARCHY -- the k-witness is the radius-1 analog of THM-523's q-witness, unifying into one family. WITNESS THM (PROVED): for a speed set S, prime p, radius r>=0, if some rotation a has NO runner s.a in {-r,..,r} mod p (all runners avoid the r-neighborhood of 0 at t=a/p), then M(S) >= (r+1)/p, witnessed by t=a/p. r=0 = the q-witness (no multiple of p -> gap>=1/p, THM-523); r=1 = the k-witness (miss +-pair {k,p-k} mod p -> t=k^{-1}/p, gap>=2/p). This is the CONSTRUCTIVE DUAL of the radius-layer over-constraint (HYP-3734/3736). TRANSVERSAL M-OPTIMALITY (point 1): the consecutive transversal {1..(p-1)/2} is the UNIQUE minimal radius-1 cover mod p -- removing any speed k unreps pair {k,p-k} and the k-witness fires (verified p=13,17,19). THREE-GAP RIGIDITY (point 3): the construction is a STRICT local M-min -- perturbing the killer breaks resonance n-1 (no longer a covering), the r=0 q-witness gives M=1/(n-1); dropping a core speed gives M=1/(n-2) (verified n=7,9,12)
status: PROVED (witness hierarchy -- direct from the definition of gap_p; verified). k-witness verified p=13,17,19 all small k. Three-gap rigidity verified n=7,9,12 (all perturbations RAISE M via a resonance break). HONEST: the witnesses prove the over-constraint NECESSITY (each radius layer must be handled) + the construction rigidity (q-witness tight at 1/(n-1)); the FULL covering-min lower bound combines all witnesses (the budget) + the binding (still open for spreads n=7..11).
source: klein-2026-06-30-S42
depends_on:
  - HYP-3736   # the killer-or-transversal budget (radius layers)
  - THM-523    # the q-witness (= r=0 of the hierarchy)
related:
  - HYP-3738   # the construction binding (three-gap); rigidity is its local-min statement
  - HYP-3742   # M-uniqueness; spreads essential
  - HYP-3734   # the radius-layer over-constraint (witnesses are its dual)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3741 — the witness hierarchy: k-witness + three-gap rigidity

## The WITNESS HIERARCHY (proved)
> **Theorem.** Let `S` be a speed set, `p` a prime, `r >= 0` an integer. If there is a rotation
> `a in (Z/p)*` such that no `s in S` has `s.a ≡ d (mod p)` for any `d in {-r, ..., r}` (every runner is
> `> r` from `0` at `t = a/p`), then `M(S) >= (r+1)/p`, witnessed by `t = a/p`.

*Proof.* `min_s ||s.(a/p)|| = (min_s dist(s.a mod p, 0))/p >= (r+1)/p`, and `M(S) = max_t min_s ||st|| >= ` this
value. ∎ (Direct from the definition of the gap.)

Two levels recover the known and the new witness:
- **`r = 0`: the q-witness (THM-523).** No runner at `0` for some `a` iff no `s ≡ 0 mod p` (no multiple of
  `p`). Then `gap >= 1/p`. With `p = q <= n` a missed resonance, this is exactly THM-523's `q`-witness
  (`M >= 1/q >= 1/n`).
- **`r = 1`: the k-witness (new).** Some `a` with no runner in `{-1,0,1}` iff `S` misses the `+-`pair
  `{a^{-1}, -a^{-1}} mod p`. Taking `k = a^{-1}`, the witness is `t = k^{-1}/p` and `gap >= 2/p`.

This is the **constructive dual** of the radius-layer over-constraint (HYP-3734/3736): the over-constraint says
"to have gap `< (r+1)/p` you must cover the `r`-neighborhood mod `p`"; the witness hierarchy *proves* it by
exhibiting the explicit lonely `t` when you do not. The witnesses live at SMALL prime moduli (`p = 13, 17, 19`).

## Transversal M-optimality (point 1, via the k-witness)
The **consecutive transversal `{1, ..., (p-1)/2}` is the unique minimal radius-1 cover of `Z/p`**: it reps every
pair `{j, p-j}` by `j <= (p-1)/2` (the dense-core lemma, HYP-3736). Remove any speed `k`: pair `{k, p-k}` is
unrepped (its other rep `p-k > (p-1)/2`), and the **k-witness fires** at `t = k^{-1}/p`, `gap >= 2/p`. Verified
`p = 13, 17, 19`, every small `k`: removing `k` -> witness rotation `k^{-1}`, min runner-distance `2`. So no
element of the consecutive transversal can be dropped without a witness -- it is rigid and optimal.

## Three-gap rigidity (point 3, via the q-witness)
The construction `{1, ..., n-2, n(n-1)}` is a **strict local M-minimum**; every perturbation RAISES `M`, and the
mechanism is the `r=0` q-witness breaking a resonance (verified `n = 7, 9, 12`):
```
killer n(n-1) -> n(n-1)±1 or +n : no longer a multiple of n-1 -> resonance n-1 MISSED
                                   -> q-witness t=1/(n-1) -> M = 1/(n-1)   (was n/Phi6)
drop a core speed (n-2 -> n-1)  : resonance broken -> M = 1/(n-2)
```
So the unit gap (the `+1` in `Phi_6 = (n-1)n+1`, created by the killer adjacency, HYP-3738) is RIGID: the
killer must be exactly a multiple of both `n-1` and `n` (minimally `n(n-1)`) or the q-witness fires. Deleting
the AP point / perturbing the killer breaks the unit-gap and the deep hole at `1/(n-1)` reappears.

## Back-and-forth: the two ends of one hierarchy
- **q-witness (`r=0`, resonances)** -> the construction's killer rigidity (three-gap, point 3). Tight: the
  perturbed `M` equals the witness value `1/(n-1)`.
- **k-witness (`r=1`, band primes)** -> the consecutive transversal optimality (point 1). The covering must
  cover every `+-`pair mod each band prime or the witness fires.
Together they are the constructive backbone of the over-constraint: a covering with small gap must handle
EVERY radius layer at EVERY relevant prime; the construction does so minimally (core transversals + one killer),
and is rigid.

## Honest scope
The witness hierarchy PROVES the over-constraint necessity (each layer must be handled) and the construction
rigidity (q-witness tight). It does NOT by itself give the full covering-min lower bound `M(S) >= M(n)` -- a
single `r`-witness gives `>= (r+1)/p`, often weaker than `M(n)` when the true binding is elsewhere (e.g. the
`n=9` construction misses pair `{8,9}` mod 17, witness `gap >= 2/17`, but its true `M = 9/73` binds at 73). The
full lower bound combines all witnesses (the budget, HYP-3736) plus the binding (proved for the construction,
HYP-3738; open for spreads `n=7..11`).

## Net
The k-witness is the radius-1 member of a WITNESS HIERARCHY whose radius-0 member is THM-523's q-witness:
runners avoiding the `r`-neighborhood mod `p` force `gap >= (r+1)/p` at an explicit `t`. This makes the
radius-layer over-constraint constructive, proves the consecutive transversal is the unique minimal radius-1
cover (transversal M-optimality), and explains three-gap rigidity (perturbing the construction trips the
q-witness, `M -> 1/(n-1)`). The covering-min is the structure that trips NO witness at any radius.
