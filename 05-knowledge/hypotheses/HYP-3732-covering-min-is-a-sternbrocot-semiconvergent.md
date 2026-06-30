---
id: HYP-3732
title: The primitive covering-min is a STERN-BROCOT SEMICONVERGENT on the ray [0;n-1,k]=k/((n-1)k+1) -- the floor 1/n=[0;n-1,1] (k=1), the covering-min=[0;n-1,a(n)] (k=a(n)), the construction n/Phi_6=[0;n-1,n] (k=n), the interval top 1/(n-1)=[0;n-1,inf] ALL lie on ONE ray indexed by depth k; the covering-min is the SMALLEST achievable depth k>1. EXACT (ILP, V=4n suffices for n<=11): M_prim is a Farey neighbor of 1/(n-1) with depth a(n)=2,2,4,4,3 for n=7..11 (covering-min 2/13,2/15,4/33,4/37,3/31; achievability NON-monotonic in k -- a=4 achievable at n=9 but a=2,3 NOT). CORRECTION (MISTAKE-088): my first pass claimed n>=12 -> 1/(n-1) (clean) with a transition at n=12 and LRC14 hard core 1/13 -- that is a V-ARTIFACT (the ILP at V=72 < n(n-1) under-resourced n>=12 and missed the construction). Since n/Phi_6 < 1/(n-1) for ALL n (n^2-n < n^2-n+1) and the construction is a valid primitive covering set, M_prim(n) <= n/Phi_6 < 1/(n-1) ALWAYS, so 1/(n-1) is NEVER the covering-min. For n>=12 the exact covering-min is UNKNOWN (<= n/Phi_6=depth-n; the small-depth spread family may or may not still beat it; needs V ~ n(n-1)). The depth a(n) is the irregular core (exact for n<=11); the ray frame is the unification of S46/S47/this
status: PARTIALLY VERIFIED + SELF-CORRECTED. SOLID: the Stern-Brocot ray frame; a(n)=2,2,4,4,3 and M_prim for n=7..11 (ILP V=4n, n=9 cross-checked vs exhaustive). RETRACTED (MISTAKE-088, V-artifact): n>=12 -> 1/(n-1), the n=12 transition, LRC14=1/13, the HYP-2566 pinning. For n>=12 the covering-min is <= n/Phi_6, exact value OPEN (ILP needs V~n(n-1); klein-S36 flagged the under-resourcing).
source: mac-mini-2026-06-30-S52
related:
  - HYP-3731  # klein-S36 + mac-mini-S51 (ceded): the covering-min ILP (the tool); klein flagged the n>=13 under-resourcing
  - HYP-3733  # mac-mini-S51 (renamed from 3731): the ILP method + OCF bridge + 4 reframings
  - HYP-3725  # the construction n/Phi_6 = [0;n-1,n] is the depth-k=n point on the SAME ray
  - HYP-3722  # the convergent [0;n-1,n] CF (the construction's ray position)
  - HYP-2566  # uniform looseness -- NOT pinned (the n>=12 pinning was the V-artifact); still open
results:
  - 04-computation/covering_min_sternbrocot_analysis_macmini_20260630.py
  - 04-computation/covering_min_ip_extended_macmini_20260630.py
  - 05-knowledge/results/covering_min_ip_extended_macmini_20260630.out
---

# HYP-3732 -- the covering-min is a Stern-Brocot semiconvergent [0;n-1,a(n)] (with a self-correction)

The owner asked to study the covering-min's irregularities, connect them to other no-closed-form sequences,
and define new sequences. The real structure is a **Stern-Brocot ray**; the irregularity is the achievable
*depth* on it. (One over-reach -- a clean `n>=12` regime -- turned out to be a search-bound artifact; recorded
honestly below and in MISTAKE-088.)

## The Stern-Brocot ray -- the unification (SOLID)
Every value of interest is a **semiconvergent** `[0; n-1, k] = k/((n-1)k+1)` on one ray indexed by depth `k`,
between the floor `1/n` and the top `1/(n-1)`:

| object | depth k | value |
|--------|---------|-------|
| the floor | `k=1` | `1/n` |
| **the covering-min** | `k=a(n)` | `[0;n-1,a(n)]` |
| the construction `n/Phi_6` (S47/HYP-3725) | `k=n` | `n/Phi_6(n)` |
| the interval top | `k -> inf` | `1/(n-1)` |

The floor, the covering-min, and the construction are the **same family at different depths** -- unifying S46
(the construction's CF `[0;n-1,n]`), S47 (construction != covering-min: depth `n` vs depth `a(n)`), and this
session. The covering-min is the **smallest achievable depth `k>1`**.

## The irregular core a(n) (EXACT for n=7..11)
`M_prim` is a **Farey neighbor of `1/(n-1)`** (`den-(n-1)num=1`), so `M_prim=[0;n-1,a(n)]`:

| n | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|----|----|
| M_prim | 2/13 | 2/15 | 4/33 | 4/37 | 3/31 |
| **a(n)** | 2 | 2 | **4** | **4** | **3** |
| margin | 1/91 | 1/120 | 1/99 | 3/370 | 2/341 |

`a(n)=2,2,4,4,3` is the **smallest Stern-Brocot depth a primitive covering set realizes**; achievability is
**non-monotone** in `k` (at `n=9`, `k=4` works but `k=2,3` do not, so `M_prim=4/33` not the mediant `2/17`).
This non-monotone achievability IS the irregularity. (ILP V=4n; klein-S36 independently confirms n=7,8,9,11;
n=9 cross-checked vs the 2M-set exhaustive.)

## The self-correction (MISTAKE-088): there is NO clean `n>=12` regime
My first pass reported `M_prim(n)=1/(n-1)` for `n>=12` (a transition at `n=12`, LRC14 hard core `1/13`, and a
pinning of HYP-2566). **That is a search-bound artifact.** The ILP used `V=72`, but the construction lives at
speed `n(n-1)` (`182` at `n=14`); `V=72` could not see it. And `n/Phi_6 < 1/(n-1)` for ALL `n` (since
`n(n-1) < n^2-n+1`), with the construction a valid primitive covering set, so `M_prim(n) <= n/Phi_6 < 1/(n-1)`
**always** -- `1/(n-1)` is never the covering-min. The ILP's `1/(n-1)` for `n=12,13,14,15` was simply the best
*low-speed* primitive set; the construction (and possibly a large-speed spread set) beats it.

**Correct status for `n>=12`:** `M_prim(n) <= n/Phi_6(n)` (the construction, depth `k=n`); the exact value is
**OPEN**. Whether the small-depth spread family (which wins for `n<=11`) keeps beating the construction at
`n>=12`, or whether the achievable depth grows toward `n`, is unknown -- the ILP needs `V ~ n(n-1)` (a much
larger universe) to decide. klein-S36 flagged exactly this under-resourcing.

## Connections to the repo's no-closed-form sequences (the owner's ask)
- **CF / Stern-Brocot family.** The covering-min joins the construction (`[0;n-1,n]`, HYP-3722/S46) and klein's
  recursive **Farey/Euclidean climb** (S32) and **Sylvester/Egyptian** tower (HYP-3724): all are `[0;n-1,k]`
  semiconvergents. The covering-min is the *shallow* achievable end; the construction the depth-`n` end. The
  whole LRC-covering story is "which rung of one ladder."
- **Extremal vs additive irregularity (a taxonomy).** The repo's no-closed-form sequences split: *additive*
  ones (W(n)=1,2,8,32,158,928 = a composition sum; tournament counts = Burnside sums) have formulas;
  *extremal* ones (the covering-min = min over covering sets; the width of `G_n` = max antichain) are
  irregular. The covering-min is extremal -> its irregularity is the non-monotone achievable depth `a(n)`.

## New sequences defined
1. **`a(n)`** -- the covering-min's Stern-Brocot depth (smallest achievable `k>1`): `2,2,4,4,3` (`n=7..11`); a
   "Diophantine achievability depth" (how deep the divisibility-covering lets the escape go). Open for `n>=12`.
2. **`M_prim(n)`** -- the primitive covering-min: `2/13,2/15,4/33,4/37,3/31` (`n<=11` exact), `<= n/Phi_6`
   (`n>=12`, open).
3. **the achievability map** `k -> (is [0;n-1,k] realizable by a primitive covering set?)` -- non-monotone,
   number-theoretic; the covering-min is its smallest `k>1`. This is the covering-min's true content.

## What it buys (and the honest caveat)
The ray turns the "irregular" covering-min into a structured object: the floor, the covering-min, and the
construction are one Stern-Brocot ray, and the irregularity is the achievable depth `a(n)` (exact, irregular,
for `n<=11`). The over-reach -- a clean `n>=12` formula -- was a V-artifact (MISTAKE-088); the `n>=12`
covering-min is open (`<= n/Phi_6`). The clean question now is the **achievability map** on the ray.
