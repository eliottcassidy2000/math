---
id: HYP-2196
title: The exactly-tight LRC family — finiteness, the 2n-1 bound, and the (Z/m)* witness orbit
status: OPEN
source: claudebox-2026-06-03-S621
depends_on:
  - THM-411   # the verified statements this seeks to prove for general n
related:
  - HYP-2195  # covering-depth master object (additive-chain characterization corrected)
  - THM-403   # AP witness orbit = (Z/m)*
  - HYP-2130  # rigidity = orbit-type (perspective key)
---

# HYP-2196 — proving the tight-family structure

[[THM-411]] establishes (exhaustively, n=3..8) that the exactly-tight (collapse / barely-lonely)
LRC family is finite, bounded by `v_max ≤ 2n-1`, counted `1,2,2,1,3,1`, and witnessed on the
`(ℤ/m)*` orbit. These should be theorems for all `n`.

## Conjectures

**(C1) Magnitude bound / finiteness.** Every primitive tight set has `v_max ≤ 2n-1`. *Heuristic:*
if `v_max = V` is large, runner `V` carves the circle into `V` lonely-eligible gaps of width
`(m-2)/(mV)`; the other `n-1` families (total `Σ_{i<n} v_i` arcs) must hit all `V` gaps to keep
`p_0 = 0`, forcing `V ≲ Σ_{i<n} v_i`; equidistribution then caps `V`. Matches the literature's
`v_n = 2n-2` tight barrier (Perarnau–Serra, normalization shift).

**(C2) Witness orbit.** Every tight set achieves `γ = 1/m` exactly on `{k/m : gcd(k,m)=1,
0<k<m/2}` and nowhere else. So a tight set is precisely a primitive integer lift of a residue
pattern that (a) is `0`-free mod `m`, (b) hits `±1`, and (c) keeps the gap pinned to `1/m` on the
whole `(ℤ/m)*` orbit with no better time elsewhere. Proving (C2) likely yields (C1).

**(C3) Even-`m` rigidity (the uniqueness dichotomy).** For `n = 3,6,8` (here) the AP is the unique
tight set; for `n = 4,5,7` there are sporadics. Conjecture this is governed by the doubling/
reflection structure of `(ℤ/m)*` (THM-404, [[HYP-2140]] 2-adic seam): sporadics exist iff the
witness orbit admits a nontrivial lift keeping all `n` arcs pinned — an orbit-rigidity condition.
Determine the exact rule and the count as an arithmetic function of `m`.

**(C4) Generating rule for sporadics.** The non-AP tights are *not* generally single-element AP
lifts (only `n=4`'s `(1,3,4,7)` is). Find the operation generating
`(1,3,4,5,9), (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)` from the AP / witness orbit. Each has
`max ∈ S+S` and shares the AP's witnesses.

**(C5) Kravitz ladder, correct normalization.** Over primitive sets the loneliness gaps below the
threshold land on `{s/(ns+1)}` for n≤5; at n=6 the value `5/33` for `(1,5,6,11,16,17)` is
off-`{s/(6s+1)}`. Resolve against Kravitz's exact statement ("barely/very lonely runners") — most
likely his `n` = #runners = `m`, so `5/33 > 1/m = 1/7` sits in the allowed `≥1/n` regime. Restate
and re-verify the ladder in the matched normalization.
