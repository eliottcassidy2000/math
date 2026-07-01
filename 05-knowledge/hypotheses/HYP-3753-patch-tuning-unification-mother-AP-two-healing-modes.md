---
id: HYP-3753
title: THE PATCH-TUNING UNIFICATION -- GW (lonely-runner FLOOR) and the covering-min construction are BOTH single-element patches of the mother AP {1..n-1} (skip element k, insert element r). The patch-tuning map M(k,r) has GENERIC value = the resonance hole 1/k (skipping k breaks resonance k), healed by special patches. TWO canonical healing modes: (A) COVERING-PATCH = skip n-1, patch lcm(n-1,n)=n(n-1) -- UNIVERSAL: always M=n/Phi6 (the covering-min construction), covers all resonances <=n, for ALL n (verified 8,10,12,14). (B) GW-PATCH = skip n-2, patch 2(n-2) (DOUBLE) -- reaches the LRC FLOOR 1/n IFF n==2 mod 6 (gcd(n-2,2)>1 and gcd(n-2,3)>1, the Jacobsthal/gcd criterion HYP-2893); verified n=8,14 hit the floor, n=10,12 give near-floor 2/19,2/23. For n=14 (==2 mod 6, the LRC14 target) the FLOOR VARIETY of single-patches is UNIQUELY {skip 12, patch 24}=GW. So the two hardest extremal objects of the two problems (tightest lonely runner; smallest covering) are single edits of the SAME AP, via two healing modes: DOUBLE (floor) and LCM (covering)
status: VERIFIED (n=8,10,12,14). Covering-patch = n/Phi6 universal (all covering); GW-patch = floor 1/n iff n==2 mod 6 (verified 8,14 hit, 10,12 miss); n=14 floor variety = {GW} uniquely (scan r<=4n). Merges the owner's patch-tuning idea with the covering-min (HYP-3738), the lowness lemma (HYP-3747), and the tight-set classification (HYP-3750).
source: klein-2026-06-30-S49
depends_on:
  - HYP-3750   # tight-set classification (floor variety = these patches)
  - HYP-2893   # Goddyn-Wong = Jacobsthal criterion (the GW-patch floor condition)
related:
  - HYP-3738   # the construction binding (= the covering-patch)
  - HYP-3747   # the lowness lemma (the core = the mother AP)
  - HYP-3745   # CRT escape (patch families)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3753 — the patch-tuning unification

## The unification
The **mother AP** is `{1, ..., n-1}` (the `n-1` speeds of a covering / the tight LRC set). A **single-element
patch** skips one element `k` and inserts one element `r`: `S = {1,...,n-1} \ {k} ∪ {r}`. Both of the problem's
extremal objects are single patches of this one AP:
- **GW** (tightest lonely runner, `M = 1/n`) `= skip (n-2), patch 2(n-2)`. For `n=14`: `{1..11,13,24}`.
- **Covering-min construction** (`M = n/Phi_6`) `= skip (n-1), patch n(n-1) = lcm(n-1,n)`. For `n=14`:
  `{1..12,182}`.
Same AP, two single edits.

## The patch-tuning map M(k, r)
Skipping `k` breaks resonance `k`, so the **generic** patch gives the resonance hole `M = 1/k` (the `q`-witness
at `D=k`). Special `r` "heal" it (verified `n=14`, `r <= 200`): skip 11 generic `1/11` (170/187 of `r`), skip 12
generic `1/12` (172/187), skip 13 generic `1/13` (173/187). Only a handful of `r` heal each skip.

## Two healing modes
- **(A) COVERING-PATCH -- LCM healing (universal).** `skip (n-1), patch lcm(n-1,n) = n(n-1)`. Verified
  `n=8,10,12,14`: always `M = n/Phi_6` (`8/57, 10/91, 12/133, 14/183`) AND covers all resonances `2..n`. This
  is the covering-min construction, available for EVERY `n`. The single big speed `n(n-1)` kills resonances
  `n-1` and `n` at once (LCM) -- the healing that yields a COVERING.
- **(B) GW-PATCH -- DOUBLE healing (arithmetic-gated).** `skip (n-2), patch 2(n-2)`. Reaches the LRC **floor**
  `M = 1/n` IFF `n ≡ 2 (mod 6)`. Proof of the condition (HYP-2893, mother size `m=n-1`): tight iff
  `gcd(n-2, j) > 1` for `j in [m-v+1, 2m-2v+1] = [2,3]` (`v=n-2`), i.e. `2 | (n-2)` and `3 | (n-2)` <=> `n` even
  and `n ≡ 2 mod 3` <=> `n ≡ 2 mod 6`. Verified: `n=8` (`≡2`, `M=1/8` floor), `n=14` (`≡2`, `M=1/14` floor);
  `n=10,12` (`≢2`, near-floor `2/19, 2/23`). **`n=14` is `≡ 2 mod 6`** -- the LRC-14 target sits exactly on the
  floor-achieving residue class.

## The floor variety is a single point (n=14)
Scanning ALL single-patches of `{1..13}` for `M = 1/14` (`r <= 4n`): the ONLY one is `skip 12, patch 24 = GW`.
So among single patches, GW is the UNIQUE floor-achiever at `n=14` -- the tightest lonely runner is a rigid,
unique single edit of the mother AP (cf. the tight-set classification HYP-3750; the cross-type tight sets are
multi-patches / lifts, not single patches).

## Extensions / the moduli picture
- The pair `(k, r)` is a **patch moduli space** over the mother AP; the M-landscape has the generic
  resonance-hole plateau `1/k` and two special healing loci: the FLOOR variety (`M=1/n`, double-healing, gated
  `n≡2 mod 6`) and the COVERING variety (kills all resonances, LCM-healing, min `= n/Phi_6`, universal).
- **Multi-patches** (skip `j` elements, patch `j`) generate the broader tight family of HYP-3750 (the
  duplication+drop near-APs, GW-type + cross-type) -- the cross-type is a 2-patch, not a single patch.
- The two healing modes are DUAL: DOUBLE (`r=2v`, minimal residue shift, the lonely-runner floor) vs LCM
  (`r=lcm(n-1,n)`, kills two resonances, the covering-min). The mother AP + {double, lcm} = the two extremal
  problems.

## Net
GW and the covering-min construction are two single-element patches of the mother AP `{1..n-1}`: GW `= skip(n-2)
patch 2(n-2)` (DOUBLE healing -> the floor `1/n`, iff `n≡2 mod 6`), the construction `= skip(n-1) patch
lcm(n-1,n)` (LCM healing -> the covering `n/Phi_6`, universal). The patch-tuning map `M(k,r)` is the resonance
hole `1/k` generically, healed at these special fibers. For `n=14` the floor variety is uniquely GW. The
tightest lonely runner and the smallest covering are two heals of the same defect in the same AP -- one problem,
two healing modes.
