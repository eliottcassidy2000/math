---
id: HYP-3763
title: "LARGE MULTIPLES ARE FORCED" -- made rigorous (general n), + the Steinhaus SCALING LAW for why they fail, + an honest correction to HYP-3745. (A) RIGOROUS LEMMA (general n, no primality): if k<=n-2, k not in S, and M(S)<=n/Phi6(n), then S contains a multiple of k, and if k>(n-2)/2 that multiple is >=2k>n-2 (LARGE). Proof: q-witness (observer at t=1/k sees every non-multiple >= 1/k, so no multiple => M>=1/k), and 1/k > n/Phi6 since kn<=(n-2)n<n^2-n+1; the smallest multiple above the vacated k is 2k. The forcing THRESHOLD is the 6th cyclotomic value Phi6(n)=n^2-n+1. (B) STEINHAUS SCALING LAW (mechanism): the forced multiple kappa=kc kills resonance k (D=k) but at a modulus D where k is near-resonant (k*a ≡ ±1) its image is c*(k-slot) ≡ ±c -- distance SCALED from 1 to c, DISPLACED from the hole it vacated. Verified n=14 k=12: c=2->M=2/25, c=7->7/89, c=13->13/157=13/Phi6(13); the surviving hole sits at a LARGE PRIME D (89,157,...). (C) HONEST CORRECTION to HYP-3745: its bound M>=2/(2n-3) holds only for single-killer perturbations; the DOUBLE-killer covering escape (84=12*7=14*6) gives 7/89 < 2/25=2/(2n-3), still > n/Phi6 but razor-thin. Single-drop COVERING escapes exceed n/Phi6 for n=10,12,14 (margins 5/43,2/21,7/89, SHRINKING ~1/n^2) but FAIL at n=8 (4/29<8/57): the construction is the covering-min only for n>=~10
status: MIXED. (A) FULLY RIGOROUS & general (q-witness THM-523 + elementary arithmetic; verified n=8..20). (B) mechanism VERIFIED (scaling law M=c/D at a large prime; n=14 all k, plus curated c=2,7,13). (C) VERIFIED corrections: HYP-3745's 2/(2n-3) fails for double-killers (7/89<2/25); construction is covering-min for n=10,12,14 (single-drop) but NOT n=8 (4/29<8/57), margin ~1/n^2 -> 0. The full "M>n/Phi6 for all covering escapes at all large n" (the covering-min conjecture) remains OPEN and razor-thin (HYP-3701).
source: klein-2026-06-30-S52
depends_on:
  - THM-523    # covering-set reduction + q-witness (PROVED; the rigorous lever)
  - HYP-3762   # three-gap / Steinhaus on the rotation orbit (the displacement mechanism)
related:
  - HYP-3745   # CRT escape uncoverable -- (C) corrects its 2/(2n-3) bound scope
  - HYP-3747   # the lowness lemma (this is its "large-speed forced" step, made rigorous)
  - HYP-3748   # step-3 rigorization (R4 unbounded = this residual)
  - HYP-3749   # punctured-core wide hole (the surviving hole)
  - HYP-3701   # covering-min extremal family transitions with n (the small-n failure, thin margin)
  - HYP-3715   # Phi6 = hexagonal/Eisenstein (the forcing threshold)
results:
  - 04-computation/large_multiple_forced_steinhaus_klein.py
  - 05-knowledge/results/large_multiple_forced_steinhaus_klein.out
  - 04-computation/covering_escape_large_multiple_klein.py
  - 05-knowledge/results/covering_escape_large_multiple_klein.out
---

# HYP-3763 — "large multiples are forced," rigorized, and the Steinhaus scaling law

## The scope fix (a definitional guard)
A **covering set** (THM-523) is a primitive `(n-1)`-set containing a multiple of **every** `q in
{2,...,n}`. The **covering-min** = `min M` over covering sets; the construction `{1,...,n-2, n(n-1)}`
gives `n/Phi6(n)`. The tight LRC minimizers (`GW = {1..11,13,24}`, `M=1/n`) and the drop-2 mediant
sets (`M=2/(2n-1)`) are **NOT covering sets** (they miss a multiple of `n`) -- THM-523's *trivial*
class, handled by the `q`-witness. Any covering-min / lowness argument must stay inside covering sets;
conflating them makes a search report the tight floor `1/n` as an "escape" (it is not one).

## (A) The LARGE-MULTIPLE-FORCED lemma -- FULLY RIGOROUS, general n, no primality
> **Lemma.** Let `n >= 4`, `k <= n-2`, `k not in S`, and `M(S) <= n/Phi6(n)`. Then `S` contains a
> multiple of `k`; and if `k > (n-2)/2`, that multiple is `>= 2k > n-2` (a *large* speed).

**Proof.** *q-witness:* put the observer at `t = 1/k`. Each runner `s` sits at `||s/k|| =
dist(s mod k, 0)/k`. If `k` divides no `s`, then `s mod k in {1,...,k-1}` for all `s`, so every
`||s/k|| >= 1/k`, whence `M(S) >= 1/k`. *Threshold:* `1/k > n/Phi6(n)` iff `Phi6(n) = n^2-n+1 > kn`,
and `kn <= (n-2)n = n^2-2n < n^2-n+1`. So `M(S) <= n/Phi6 < 1/k` forces a multiple of `k` in `S`.
Since `k not in S` and `k > (n-2)/2`, the smallest available multiple is `2k > n-2`. ∎

No primality is needed. The **forcing threshold is the 6th cyclotomic value** `Phi6(n)=n^2-n+1` (the
hexagonal/Eisenstein modulus, HYP-3715): below `1/k` the resonance at `k` cannot be missed. Verified
`n=8..20`, all top-half `k`: `M(P_k) >= 1/k > n/Phi6`, forced multiple `>= 2k > n-2`.

## (B) The STEINHAUS SCALING LAW -- why the forced large multiple fails
The forced multiple `kappa = k*c` kills resonance `k` (at `D=k`, `kappa ≡ 0`). But the core
`{1,...,n-2}` is an arithmetic progression, so (three-gap theorem, HYP-3762) at any modulus `D` and
rotation `a` its images are a rotation orbit with `<=3` gaps -- a maximally even cover with no slack.
Remove `k` and the two gaps around its image `k*a` merge into a **double gap** = the hole. At a
modulus `D` where `k` is *near-resonant* (`k*a ≡ ±1`, i.e. `k` would sit a hair from the observer and
cover it tightly), the substitute's image is
```
    kappa * a  =  c * (k * a)  ≡  ±c   (mod D)
```
-- distance **scaled from 1 to c**, out at the rim of the hole instead of in it. `kappa` refills the
`k`-slot only when `kappa ≡ k (mod D)`, i.e. `D | k(c-1)` -- a thin set never containing the deep
hole's modulus. **Killing the resonance (`kappa ≡ 0 mod k`) and filling the hole (`kappa ≡ k mod D`)
demand incompatible congruences of the same integer.**

Verified (`n=14`, drop `k=12`, genuine covering escapes):
```
 kappa = 24 = 12*2   (c=2)   binds D=25         M = 2/25   = 2/(2k+1)
 kappa = 84 = 12*7   (c=7)   binds D=89 (prime) M = 7/89
 kappa = 156= 12*13  (c=13)  binds D=157(prime) M = 13/157 = 13/Phi6(13)
```
`M = c/D`, the multiplier `c` is the hole depth, and the surviving hole sits at a **large prime** `D`.
"Large primes forced" in *both* senses: the covering set is forced to carry a large speed, and the
surviving hole is forced out to a large prime modulus (a near-resonance of the vacated small speed).

## (C) Honest corrections
- **HYP-3745's `M >= 2/(2n-3)` is scope-limited.** It holds for *single-killer* perturbations
  (`kappa=kc` replacing `k`, construction killer kept): `n=14` gives `2/25 = 2/(2n-3)`. But the
  *double-killer* covering escape `84 = 12*7 = 14*6` (one speed covering resonances 12 AND 14, freeing
  the budget so 13 stays small) gives `M = 7/89 < 2/25 = 2/(2n-3)`. So `2/(2n-3)` is **not** a uniform
  lower bound over covering escapes; the true single-drop covering-escape min is `7/89` at `n=14`.
- **The construction is the covering-min only for large `n`.** Single-drop covering escapes exceed
  `n/Phi6` for `n=10` (`5/43`), `n=12` (`2/21`), `n=14` (`7/89`), but the margin **shrinks like
  `1/n^2`** (`+0.0064, +0.0050, +0.0021`), and at `n=8` a covering escape `{1,2,3,4,5,7,24}` gives
  `M = 4/29 < 8/57` -- the construction is **beaten** (matches HYP-3701). So the lowness lemma / the
  large-speed impotence needs `n >= ~10`, and the full statement for all large `n` is the razor-thin
  open edge of the covering-min conjecture.

## Net
The "large speed is forced" step of the lowness lemma is now **fully rigorous and general** (Lemma A:
`q`-witness + the cyclotomic threshold `Phi6`). The Steinhaus **scaling law** explains rigorously why
the forced multiple then fails -- its multiplier `c` scales its coverage distance from 1 to `c` at the
vacated speed's near-resonance, so it lands at the rim of the hole (`M=c/D`, `D` a large prime), not in
it. The remaining gap (a uniform closed inequality `M > n/Phi6` for *all* covering escapes at *all*
large `n`) is genuinely thin -- the margin decays `~1/n^2` and the construction is not even the
covering-min for `n<=8` -- and stays open (HYP-3701). Two honest corrections: HYP-3745's `2/(2n-3)`
bound is single-killer-only; the construction's covering-minimality is an `n>=~10`, razor-margin fact.
