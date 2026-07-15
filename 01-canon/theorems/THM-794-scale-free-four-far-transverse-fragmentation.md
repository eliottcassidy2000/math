---
id: THM-794
title: Scale-free four-far transverse fragmentation and capped closure
status: PROVED (elementary interval geometry and capped tail) + VERIFIED-EXACT (176-prime capped base and component canaries attached)
source: codex-2026-07-14-S2
renumber_note: claimed locally as THM-779, then THM-783, then THM-792;
  renumbered to THM-794 after the shared THM-792/793 claims landed.
depends_on:
  - THM-755   # the peel-relative capped-envelope coordinate
  - THM-795   # invariant-level good-set-state transverse-tooth cap
related: [HYP-6780, HYP-6815, HYP-6830, MISTAKE-145]
verification: 04-computation/lrc14_affine_slope_suspension_codex_S2.py
  (+ 05-knowledge/results/lrc14_affine_slope_suspension_codex_S2.out)
---

# THM-794 — Scale-free four-far transverse fragmentation and capped closure

## Statement

For every prime `N>110`, define

```text
B   = {1,2,...,9,15,110},
P_N = B union {N},
V_N = P_N union {1092N}.
```

Then:

1. `V_N` is a primitive covering thirteen-speed row in the exactly-`f=4`
   stratum, with nine-core `{1,...,9}` and far speeds
   `{15,110,N,1092N}`.
2. No nontrivial integer divides seven members. The largest divisor packet has
   size five in `P_N` and six in `V_N`. Under HYP-6830's convention that
   `c*(S)` is the largest scale dividing at least seven speeds, or one if none
   exists, `c*(P_N)=c*(V_N)=1`.
3. If `r_+(P)` denotes the number of positive-length interval components of

   ```text
   G'_P={t in R/Z : ||pt||>=1/14 for every p in P},
   ```

   then

   ```text
   r_+(P_N) >= N/1540 - 8/7.
   ```

   In particular `r_+(P_N)` is unbounded although the maximal high-support
   divisor scale is identically absent.
4. Nevertheless every `V_N` is closed by its top peel: the THM-755 capped
   envelope fires at `v=1092N`, and hence `M(V_N)>1/14`. This holds for the
   whole prime family, not only the four component-count canaries.

Thus good-set fragmentation is not controlled by the maximal high-support
divisor scale, even inside the first open four-dimensional LRC14 chart.

## Proof

Primitivity is immediate from `1 in V_N`. Moduli `2,...,9` divide their named
speeds, `10` and `11` divide `110`, and `12,13,14` divide
`1092N`, since `1092=lcm(12,13,14)`. Hence `V_N` is covering. Exactly the
four displayed far speeds exceed `14`.

For divisor packets, `2` divides exactly
`{2,4,6,8,110}` in `B`, while `3` divides four members. For every `d>=4`, at
most `floor(9/d)<=2` members of `{1,...,9}` are divisible by `d`, and the two
extra speeds can add at most two, so no such packet is larger than four. The
prime `N>110` shares no nontrivial divisor with `B`. Therefore the maximum in
`P_N` is five. Adding `1092N` adds at most one member to every old packet;
packets carrying the prime `N` contain only `N` and `1092N`. The full-row
maximum is consequently six, attained by the even packet.

It remains to prove transverse fragmentation. The fixed core good set contains
the interval

```text
J=[1/14,111/1540],            |J|=1/1540.
```

Indeed, on `J` the phases of `1,...,6` remain in safe affine branches below
`1/2`; the phases of `7,8,9` lie between `1/2` and `999/1540`; the fractional
part of `15t` lies between `1/14` and `125/1540`; and the fractional part of
`110t` lies between `6/7` and `13/14`. Every circular distance is therefore at
least `1/14`.

The `N`-runner is dangerous on the disjoint open teeth

```text
T_k=((k-1/14)/N,(k+1/14)/N),       k in Z.
```

If `N<=1760`, the asserted lower bound `N/1540-8/7` is nonpositive and hence
trivial. Assume for the remainder of this count that `N>1760`.

A tooth has closure strictly inside `J` whenever its center `k/N` is farther
than `1/(14N)` from both endpoints. After multiplying by `N`, the admissible
center interval has length

```text
N|J|-2/14 = N/1540-1/7.
```

An open real interval of length `L` contains at least `L-1` integers. Hence at
least `N/1540-8/7` full, pairwise disjoint danger teeth lie inside `J`. Each is
an open gap in `G'_{P_N}`. In the circle, only the two end pieces of `J` can
possibly reconnect outside `J`; the internal teeth still give at least as many
global positive-length good-set components as full gaps. This proves the stated
lower bound and the unboundedness. ∎

For the capped closure, first handle the finite base by exact rational interval
arithmetic. The companion script checks all `176` primes

```text
113 <= N < 1265
```

and verifies `PI_LO*(1092N)*|G'_{P_N}|>r_+(P_N)` in every case, where
`PI_LO=333/106<pi`. The smallest exact ratio `v|G'|/r_+` in this base is

```text
8523034/28875, attained at N=113.
```

The infinite tail is elementary. At most `N|J|+2=N/1540+2` of the
`N`-danger teeth can meet `J`, and each has length `1/(7N)`. Therefore

```text
|G'_{P_N}| >= 1/1540-(N/1540+2)/(7N)
             = 3/5390-2/(7N).
```

The union of the danger arcs of runners `p in P_N` has at most `sum(P_N)`
components. Hence both the full topological count `r_top` and the
positive-length count satisfy

```text
r_+(P_N) <= r_top(P_N) <= sum(P_N) = N+170.
```

Consequently

```text
PI_LO*(1092N)*|G'_{P_N}|-r_+(P_N)
 >= (PI_LO*1092*6/(7*1540)-1)N
    -(PI_LO*1092*2/7+170),
```

whose unique zero is

```text
11734415/9278 < 1265.
```

Thus the strict capped-envelope inequality holds for every integer `N>=1265`.
Together with the exact prime base, it holds for every prime `N>110`, proving
the fourth claim. ∎

## What the theorem does and does not say

This is not a counterexample to LRC14: the top peel closes the entire family.
The count `r_+` is the THM-755/Fourier convention: isolated equality points
have zero measure and create no jump interval. At the four component-count
canaries `N=211,503,1009,2003`, the positive-length counts are respectively
`66,104,174,310`, while the full closed-set topological counts are
`68,108,176,312`; the differences `2,4,2,2` are isolated equality points. The
exact ratios `v|G'_P|/r_+(P)` are

```text
63068083/152460, 727151/1155,
76395787/100485, 11045593/13020
```

The theorem refutes the proposed implication `r_+(P)<=B(c*(P))`; it supports a
peel-relative coordinate such as `r_+(P)/(v|G'_P|)`, not a raw component bound.
Its universal capped closure is specific to this family and is not a finiteness
or descent theorem for the general exactly-`f=4` chart.

## Tournament analysis and assumption challenge

Runner or divisor-packet vertices do not carry the proof. The operative
vertices are base safe components and inserted teeth, with incidence meaning
“this tooth cuts this component.” That object is bipartite rather than a
natural tournament: orienting runner pairs loses the number, width, and common
location of the cuts. Alternate vertices considered were runners, gaps,
section boundaries, endpoint events, residues, Fourier modes, divisor packets,
and peel obligations.

The companion HYP-6815 audit uses representations as tournament vertices. Its
pairwise observable is proof-critical row-pair separation; predicate-first and
compression-first gauges have a declared tie Hamiltonian path and differ by
eleven edge flips. That tournament is a diagnostic of information loss. The
proof carrier here is the owner-marked component/tooth incidence plus the named
peel, not the tournament orientation.
