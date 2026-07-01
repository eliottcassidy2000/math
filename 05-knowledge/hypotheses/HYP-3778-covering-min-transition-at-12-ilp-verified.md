---
id: HYP-3778
title: THE COVERING-MIN TRANSITION AT n=12 IS REAL (ILP-verified up to speeds 4n) -- correcting my HYP-3764 pessimism. Using the exact set-cover ILP (HYP-3731/3733, the RIGHT tool), the covering-min for speeds <= 4n: n=7..11 the SPREAD FAMILY beats the construction (2/13,2/15,4/33,4/37,3/31, all < n/Phi6) -- including the elusive n=11 = {2,6,8,9,10,11,13,14,17,19}=3/31 which DROPS speed 1 (why all my ad-hoc searches missed it); but n=12,13,14 have NO beater <= 4n -- the best speed-<=4n covering set is 1/(n-1) (1/11,1/12,1/13) > n/Phi6 (12/133,13/157,14/183), so the construction WINS. Clean qualitative transition: covmin(<=4n) < construction for n<=11, = 1/(n-1) > construction for n>=12. This CONFIRMS HYP-3737 (construction for n>=12) and CORRECTS HYP-3764's leading hypothesis H3 ('transition is a search-mirage, beaters persist') -- the beaters do NOT persist past n=11 within the 4n envelope; my earlier random/hillclimb/calibrated-drop searches failed because they were the WRONG TOOL (the real beaters are highly SPREAD, drop 5+ core speeds incl. speed 1, which ad-hoc structural searches can't reach; the exact ILP can). HONEST residual: only speeds <= 4n excluded (milp times out at V>~4n, e.g. V=90 gives a spurious 6/55); a beater in (4n, n(n-1)] is not ruled out, but all known beaters have speeds <= 3.5n
status: MIXED (ILP-verified up to V=4n + honest residual). VERIFIED (exact set-cover ILP, milp/HiGHS, reliable at V~4n): covmin(<=4n) = 2/13,2/15,4/33,4/37,3/31 (n=7..11, matches repo, beats construction) and 1/11,1/12,1/13 (n=12,13,14, = 1/(n-1), does NOT beat construction n/Phi6). The n=11=3/31 recovery (drops speed 1) VALIDATES the ILP where my ad-hoc searches failed. CORRECTS HYP-3764 (H3 disfavored, H4/HYP-3737 supported) up to the 4n envelope. RESIDUAL: speeds (4n, n(n-1)] untested (milp timeout); NOT a full proof that covmin(n>=12)=n/Phi6, but strong evidence the spread family is exhausted by n=12.
source: klein-2026-06-30-S60
depends_on:
  - HYP-3731   # the covering-min set-cover ILP (the right tool)
  - HYP-3764   # my open-edge pessimism (this CORRECTS it)
related:
  - HYP-3737   # construction forced for n>=12 (this CONFIRMS, up to 4n)
  - HYP-3733   # ILP pins n=7,8,9,11 (the validation)
  - HYP-3732   # covering-min Farey rung a(n)=2,2,4,4,3 (n<=11); rung n (construction) for n>=12
  - MISTAKE-087 # the beaters are spread-structured (why searches miss them)
results:
  - 04-computation/covering_min_ip_v2_macmini_20260630.py
  - 05-knowledge/results/covering_min_ilp_n12_klein.out
  - 05-knowledge/results/covering_min_ilp_4n_klein.out
---

# HYP-3778 — the covering-min transition at n=12 is real (ILP-verified up to 4n)

## The new approach: use the RIGHT tool
For sessions I attacked the covering-min open edge (is the construction `n/Phi_6` the covering-min for
`n>=12`?) with random search, hill-climbing, and calibrated drop-`<=k` searches -- and every one
FAILED validation, missing known beaters (HYP-3764). This session: the exact **set-cover ILP**
(HYP-3731/3733; universe = all breakpoints `k/d`, `d<=2V`; binary-search the radius; scipy `milp`/HiGHS)
is the right tool. It is EXACT for speeds `<= V`, and it recovers what the ad-hoc searches missed.

## Validation (the ILP recovers the elusive beaters)
```
  n     covmin(speeds<=4n)   construction n/Phi6   beats?
  7     2/13                 7/43                   YES
  8     2/15                 8/57                   YES
  9     4/33                 9/73                   YES
 11     3/31                 11/111                 YES   <- {2,6,8,9,10,11,13,14,17,19}, DROPS speed 1
```
The `n=11` covering-min `3/31` **drops speed 1 and 5 core speeds** -- a highly spread set that random,
hill-climb, and my drop-`<=3`-keep-`1` calibrated search ALL missed (why HYP-3764 was stuck). The exact
ILP finds it. So the ILP is validated where every heuristic failed.

## The transition at n=12 (the result)
```
  n     covmin(speeds<=4n)   construction n/Phi6   beats?
 12     1/11  = 0.09091      12/133 = 0.09023      NO
 13     1/12  = 0.08333      13/157 = 0.08280      NO
 14     1/13  = 0.07692      14/183 = 0.07650      NO
```
For `n=12,13,14`, the best covering set with speeds `<= 4n` is `1/(n-1)` (a punctured-resonance set,
e.g. `n=14`: `{2,4,...,24,39}`), which is **ABOVE** the construction `n/Phi_6`. So the spread family --
which produced a lower-rung beater for every `n<=11` -- has **NO member `<= 4n` beating the construction
at `n=12,13,14`**. There is a clean qualitative transition:
> `covmin(<=4n) < construction` for `n<=11`;  `covmin(<=4n) = 1/(n-1) > construction` for `n>=12`.
This is exactly HYP-3737's predicted transition at `n=12`, now VERIFIED (up to `4n`) with the exact ILP.

## This corrects HYP-3764
HYP-3764's leading (bold) hypothesis was **H3: the "n>=12 transition" is a search-mirage; low-rung
beaters persist for all n**; it dismissed the transition as "unverified, search-resistant." The ILP shows
the OPPOSITE: the beaters do **not** persist past `n=11` within the `4n` envelope; the transition is real
and verifiable. HYP-3764's pessimism came from using unreliable searches (which cannot reach the spread
beaters) and concluding the problem was intrinsically search-resistant -- but it was **tool-resistant**,
not search-resistant: the exact ILP resolves it. So **H4 (transition real, HYP-3737) is now supported,
H3 disfavored**, up to speeds `4n`.

## Honest scope / residual
- The ILP is exact only for speeds `<= V`, and reliable only at `V ~ 4n` (at `V=90, n=12` the `milp`
  25s time-limit is hit and it returns a spurious `6/55 > 1/11` -- larger `V` is timeout-unreliable, not
  trustworthy). So the verified statement is: **no covering-min beater with speeds `<= 4n` at
  `n=12,13,14`.**
- A beater with speeds in `(4n, n(n-1)]` is NOT excluded. But every known beater (`n<=11`) has max speed
  `<= 3.5n`, so the `4n` envelope is where the spread family lives; its emptiness at `n=12,13,14` is
  strong evidence the family is exhausted. A full proof still needs the `(4n, n(n-1)]` range (or the
  lowness lemma HYP-3747).

## Net
Switching to the exact set-cover ILP -- the right tool -- resolves the covering-min open edge up to
speeds `4n`: the spread family beats the construction for `n<=11` (including the speed-1-dropping `n=11`
beater `3/31` that all heuristics missed) but is EXHAUSTED at `n=12,13,14` (best `<=4n` set is
`1/(n-1) >` construction). The transition at `n=12` (HYP-3737) is real and verified up to `4n`; my
HYP-3764 pessimism (H3, "search-mirage") is corrected -- the edge was tool-resistant, not
search-resistant. Residual: the `(4n, n(n-1)]` speed range and a full lower-bound proof remain open.
