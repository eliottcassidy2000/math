---
id: HYP-3724
title: Phi_6 IS THE SYLVESTER MAP -- the covering-modulus map Phi_6(x)=x^2-x+1 iterated from 2 gives SYLVESTER's sequence 2,3,7,43,1807,...; the apex primes 2,3,7,43 ARE Sylvester numbers (7=s_3 is the LRC(14) apex, the genus-1 boundary; 1807=13*139 is the first composite, where the prime-apex tower breaks). The recursive product structure: prod_{j<k} s_j = s_k-1, so the killer at apex s_k = s_k(s_k-1)=prod_{j<=k}s_j=s_{k+1}-1 (the product of ALL apexes so far) and the covering modulus = s_{k+1} (the next Sylvester); the apex covering-min M(s_k)=s_k/s_{k+1} (consecutive Sylvester ratio). EGYPTIAN: sum 1/s_k = 1 (the greedy/fastest unit-fraction for 1), remainder after N = 1/(s_{N+1}-1)=1/prod_{j<=N}s_j (doubly exponential). With Phi_6=2T+1 (HYP-3723): Sylvester = iterated 'twice-triangular-plus-one'; the LRC covering = the Sylvester/Egyptian dynamical system on the triangle
status: VERIFIED (Sylvester via Phi_6: 2,3,7,43,1807; sum 1/s_k -> 1 with remainder 1/(s_{N+1}-1); prod_{j<k}s_j=s_k-1; killer=s_k(s_k-1)=s_{k+1}-1; M(n)=n/Phi_6(n)=1/[n-1;n]). Synthesis (recursive/abstract). The clean Sylvester recursion is at the APEX (n=apex prime); LRC(14) uses N=14 (not itself Sylvester) but its apex 7=s_3 is on the tower.
source: klein-2026-06-29-S31
depends_on:
  - HYP-3723   # Phi_6 = 2T+1; the functional decomposition; the covering modulus = twice the staircase
  - HYP-3715   # M = n/Phi_6(n); killer = lcm(n-1,n)
related:
  - HYP-1865   # prior LRC Egyptian-split ledger (n=16)
  - HYP-2442   # Erdos-Moser tower
  - HYP-3705   # Phi_6 = Eisenstein norm (the Q(sqrt-3) column)
results:
  - 04-computation/torus_lift_threegap_klein.py
---

# HYP-3724 — Phi_6 is the Sylvester map; the apex tower is the Egyptian-fraction tower

## Phi_6 = the Sylvester map
The covering-modulus map `Phi_6(x) = x^2 - x + 1`, iterated from `2`, IS Sylvester's sequence:
> `s_1 = 2`, `s_{k+1} = Phi_6(s_k) = s_k^2 - s_k + 1`  ->  `2, 3, 7, 43, 1807, 3263443, ...`
The **apex primes `2, 3, 7, 43` ARE the Sylvester numbers**; `7 = s_3` is the LRC(14) apex (the genus-1
boundary, where `genus(X_0(2p))` jumps to 1). `s_5 = 1807 = 13 . 139` is the first COMPOSITE Sylvester
number -- where the *prime*-apex tower breaks (so the LRC genus-tower of prime apexes is `2,3,7,43`, length 4).

## The recursive product structure (the killer is the product of all previous apexes)
Sylvester's defining property `s_{k+1} = (prod_{j<=k} s_j) + 1` gives, with `prod_{j<k} s_j = s_k - 1`:
> killer at apex `s_k` `= s_k(s_k - 1) = prod_{j<=k} s_j = s_{k+1} - 1`,   modulus `= s_{k+1}`,
> apex covering-min `M(s_k) = s_k / s_{k+1}` (a consecutive Sylvester ratio).
So the covering structure at apex `s_k` is built from the PRODUCT of all apexes up to `s_k` (the killer) plus
one (the modulus = the next apex). The tower multiplies in one apex per level. (Verified `s_k=2,3,7,43`.)

## The Egyptian fraction (the conserved quantity)
> `sum_k 1/s_k = 1`,   `1 - sum_{k<=N} 1/s_k = 1/(s_{N+1} - 1) = 1/prod_{j<=N} s_j`.
The apex reciprocals `1/2 + 1/3 + 1/7 + 1/43 + 1/1807 + ... = 1` -- the GREEDY (Sylvester-Fibonacci) unit
fraction for `1`, the *fastest* Egyptian fraction (doubly-exponential remainder). So the apex tower exactly
TILES the unit by unit fractions; the LRC covering is the greedy Egyptian covering, and `M > 1/n` is the
greedy over-shoot (`M(s_k) = s_k/s_{k+1} > 1/s_k`).

## The covering-min reciprocal
`M(n) = n/Phi_6(n) = 1/((n-1) + 1/n) = 1/[n-1; n]` -- its reciprocal is the length-2 continued fraction
`(n-1) + 1/n`. At an apex `n = s_k`, `M = s_k/s_{k+1}`.

## The recursive/abstract picture (everything iterates Phi_6 = 2T+1)
Chaining with HYP-3723 (`Phi_6 = 2 T_{x-1} + 1`):
> triangle `T` -> `Phi_6 = 2T+1` -> iterate -> SYLVESTER `2,3,7,43,...` (the apex tower) -> EGYPTIAN
> `sum 1/s_k = 1` (the complete covering).
The LRC covering structure is the orbit of the Sylvester map `Phi_6` (= the triangle, doubled & shifted),
with the Egyptian fraction `sum 1/s = 1` as the conserved "total covering budget," tiled one apex at a time.
The operation algebra `{even/odd, +1/-1, x2 / /2}` (HYP-3723) generates each step; the iteration generates
the tower. Recursively and abstractly: **the lonely-runner covering is a Sylvester/Egyptian dynamical
system on the triangular numbers.**

## Honest scope
- VERIFIED: `Phi_6` = Sylvester map; apex primes = Sylvester; `sum 1/s_k = 1` with the Sylvester remainder;
  killer `= prod` of previous apexes `= s_{k+1}-1`; `M(n) = 1/[n-1;n]`.
- The clean Sylvester recursion is at the APEX (`n = apex prime`); `LRC(14)` uses `N = 14` (not itself
  Sylvester), but its apex `7 = s_3` is on the tower and `Phi_6(14) = 183` uses the same map. The prime-apex
  genus tower is finite: `2, 3, 7, 43` (then `1807` composite).
- This is structural/synthetic (a unifying recursive frame), not a new bound; it explains WHY `Phi_6` and the
  apex `7` are the natural objects (the Sylvester/Egyptian tower) and ties the covering layer to the triangle.
