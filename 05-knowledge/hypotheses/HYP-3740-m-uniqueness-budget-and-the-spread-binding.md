---
id: HYP-3740
title: M-UNIQUENESS via the killer-or-transversal budget + the spread-binding open piece. NEGATIVE RESULT: the clean rung-k construction (dense core {1..n-2} + CRT killer) FAILS for k<n -- the forced-large killer creates a deeper hole at a FINER modulus (n=7,k=2: killer 168=13^2-1 gives M=28/169 at binding 13^2, not 2/13), so CLEAN three-gap constructions exist ONLY at rung n; low rungs need spreads (messy gaps). M-UNIQUENESS: (proved for n>=12) the radius-1 band budget forces the construction (mac-mini HYP-3737) and its binding is D=Phi6≡1 mod(n-1) (klein HYP-3738), so M(n)=n/Phi6 is the UNIQUE covering-min, binding D=Phi6, for n>=12 incl. n=14=14/183; (all n) the value M(n) fixes a unique binding D=denom(M) and rung shared by every extremal covering (the covering itself is NOT unique: n=7 has 2). The SPREAD-regime binding (D≡1 for n=7..11) stays OPEN -- reframed: the covering-min is conjecturally a best one-sided approximation to 1/(n-1)
status: NEGATIVE RESULT verified (clean rung-k construction fails for k<n: n=7,8,9). M-UNIQUENESS proved for n>=12 (combining mac-mini's band-forcing + klein's construction-binding); invariant-uniqueness holds for all n (value determines binding/rung). Spread binding OPEN (clean-construction route now closed).
source: klein-2026-06-30-S41
depends_on:
  - HYP-3738   # the construction binding (PROVED: D=Phi6, three-gap)
  - HYP-3736   # the killer-or-transversal budget
related:
  - HYP-3737   # mac-mini: radius-1 band forces construction for n>=12
  - HYP-3734   # Farey-neighbor <=> D≡1 mod(n-1)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3740 — M-uniqueness via the budget; the spread-binding open piece

## NEGATIVE RESULT: clean three-gap constructions exist only at rung n
To prove the spread binding I tried to GENERALIZE the construction: dense core `{1,...,n-2}` + a "rung-k
killer" `kappa` solving `kappa.k ≡ (n-2)k+1 mod D` (lands one above the top core point at rotation `a=k`) and
`kappa ≡ 0 mod n(n-1)` (kills resonances `n-1, n`). At `a=k` this gives the clean three-gap `{1, k^(n-3), 2k}`
on `Z/D`, `D=k(n-1)+1` -- exactly the construction structure, which WOULD prove `D ≡ 1 mod (n-1)`.

**It fails for `k < n`.** The CRT forces `kappa` LARGE, and the large killer creates a deeper hole at a FINER
modulus, moving the binding off `D`:
```
n=7 k=2: kappa=168 (=13^2-1)  M=28/169  (binding 13^2, not 2/13)   FAIL
n=7 k=3: kappa=588            M=98/589                              FAIL
n=8 k=2: kappa=224            M=32/225  (binding 15^2)              FAIL
n=9 k=2: kappa=288            M=36/289  (binding 17^2)              FAIL
n=9 k=3: kappa=1224           M=153/1225                           FAIL
```
The binding lands at `D^2` because `kappa ≡ -1 mod D^2`-type, so `kappa.a ≡ -a` makes a slope-`-1` runner that
digs a deep hole at the fine modulus. Only `k=n` (`kappa = n(n-1)`, the MINIMAL killer) avoids this -- that is
why the construction is special and the three-gap binding (HYP-3738) lives at rung `n` alone. **Low rungs have
no clean construction; they are achieved only by spreads** (which kill band primes by including the small
prime `p` itself, not a large multiple -- e.g. `n=11` uses speeds `13,17,19`), whose binding gaps are messy
(the `n=7` spread has gaps `{1,1,2,2,3,4}`, not `{1,2,2,2,2,4}`). So the spread binding cannot be reduced to
the construction proof.

## M-UNIQUENESS
**(i) n >= 12 -- PROVED.** mac-mini's radius-1 band over-constraint forces the construction for `n >= 12`
(HYP-3737), and the construction binding is `D = Phi_6 = (n-1)n+1 ≡ 1 mod (n-1)` (klein HYP-3738, proved by the
explicit `{1, n^(n-3), 2n}` three-gap). Together: **for `n >= 12` the covering-min is uniquely
`M(n) = n/Phi_6(n)`, with unique binding `D = Phi_6`** -- in particular `n=14: M = 14/183`, `D = 183 = 14.13+1`,
rung 14.

**(ii) all n -- invariant uniqueness.** The covering-min VALUE `M(n) = j/D` (lowest terms) fixes a unique
binding modulus `D = denom(M(n))` and rung `k = (D-1)/(n-1)`; every extremal covering shares them. The
extremal COVERING is NOT unique -- `n=7` has exactly two (`{1,2,5,6,7,8}` and `{1,4,5,6,7,11}`), both `D=13`.
So the budget's up-set structure (rung feasibility is monotone, HYP-3736) gives a unique `k_min`, and `M(n)`
is the unique value at that rung -- the right invariant is the rung, realized by possibly several coverings.

## The spread-binding open piece (reframed)
For `n = 7..11` the covering-min is a SPREAD at rung `k_min < n`, with `D = k_min(n-1)+1 ≡ 1 mod (n-1)`
VERIFIED but unproved. The clean-construction route is now CLOSED (negative result above). The natural
remaining frame: the covering-min is conjecturally a **best one-sided rational approximation to `1/(n-1)` from
below** within `(1/n, 1/(n-1))` -- and those are exactly the semiconvergents / ladder fractions
`k/(k(n-1)+1)`. Proving the covering-min is a best approximation (not merely a fraction in the interval) would
close `D ≡ 1 mod (n-1)`. The over-constraint/budget gives the rung lower bound (HYP-3736) but not yet the
on-the-ladder property for spreads.

## Net
The clean rung-k construction fails for `k < n` (the forced-large killer digs a deeper hole at `D^2`), so the
three-gap binding proof is rung-`n`-only and the spreads need a different argument. M-uniqueness is PROVED for
`n >= 12` (band-forcing + construction-binding -> `M = n/Phi_6` unique, `D = Phi_6 ≡ 1`, incl. `n=14`); for all
`n` the binding `D` and rung are unique invariants (the covering is not). The spread binding (`n=7..11`) stays
open, now framed as "the covering-min is a best one-sided approximation to `1/(n-1)`."
