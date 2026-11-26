# Mixed spine reversals and activity transfer on a closed branch

2026-11-26. Continues [the exterior-closure bound](EXTERIOR-CLOSURE-AND-MODULE-RUNS.md). The common step is to determine exactly what an edit passes along a boundary before bounding the number of states on the other side.

## 1. An unbounded family of unit reversals, but not a perfect-number crossing

Use the F1/F2 seam $H=hF-eE+RR'$: the left rooted profile $(F,E,R)$ counts all paths, root-starting paths and ordered two-path covers with first path ending at the root and optional second path; the right dual profile is $(h,e,R')$. Spine sites on the left are numbered in their increasing arc order, with a 1-bit pointing into the root. Address 9 has bit word 1000 (site 3 is the 1). Reverse physical arc $1\to3$. The exact profile change is $(0,2,0)$, so this edit lowers the seam count by $2e$ for any actual right core. In particular it lowers it by two on **every odd binary right address**, with unbounded central counts. For $c=2^a+1$, $a\ge1$, the center is $18\cdot2^a+17$, using $F(2^a)=2^a+1$.

Here is a formula valid for mixed bits and arbitrary distance. For $S$ containing both flipped sites, let $T$ be its complement in the left spine. Let $d(S),a(S),k(S)$ count paths on $S$ ending at a 1-bit, starting at a 0-bit, and unrestricted, respectively. Each takes value one on the empty set, representing an absent run. Subtract the increasing path exactly when $S$ has no interior site. Add all paths using the reversed arc: put any subset of the interior sites before the arc and the complementary subset after it, and retain exterior sites below and above in the prefix and suffix respectively. Each part is increasing. This gives $(\Delta d,\Delta a,\Delta k)$ using the actual first and last sites. No two disjoint runs can both contain the flipped pair. Summing at the seam root gives

    ΔF = Σ [Δd(S) a_old(T) + d_old(T) Δa(S)],
    ΔR = Σ [Δd(S) + Δk(S) d_old(T)],
    ΔE = Δa(whole spine).

For address 9 and masks 10,11,14,15, the individual contributions to $(\Delta F,\Delta R)$ are $(−1,−1),(−1,−1),(0,0),(2,2)$; the full-spine $\Delta a$ is two. Optional high zero padding leaves all three quantities invariant, since those sites are ordered sinks.

For a unit edit via this construction we must have $e=1$, so the center is $9(h+R')-1$. The midpoint is always **7 mod 9**, even for nonbinary right cores. An odd-exponent even perfect number has midpoint $U(2U-1)$ with $U=2^{p-1}\equiv1\pmod3$, hence is always **1 mod 9**. So this entire infinite lane misses those midpoints, including exponent 61. Doubling the supply of one-sided unit exits need not help at the prescribed integer.

The replay checks the mixed-bit formula against 2,560 literal orientations and verifies eight spine-plus-spoke forks at six centers:

| Left, right addresses | Spine edit | Spoke bit | Center → spine, spoke |
|---|---|---|---|
| 7,1 | 0→2 | 0 or 1 | 11 → 13,9 |
| 9,1 | 0→3 | 0 | 17 → 19,15 |
| 9,1 | 1→3 | 1 | 17 → 15,19 |
| 5,3 | 0→2 | 0 | 29 → 31,27 |
| 25,1 | 2→4 | 1 | 47 → 49,45 |
| 73,1 | 4→6 | 2 | 141 → 143,139 |
| 289,1 | 6→8 | 3 | 573 → 575,571 |

A genuinely different opposite exit would be an edit inside the right core satisfying $9\Delta(h+R')-\Delta e=2$ starting at $e=1$. Thus the edited terminal count must be 8 mod 9. An exhaustive literal orientation census found no such edit on right cores of order at most seven (64,080 order-seven orientations had one terminal path at the marked root, generating 1,345,680 directed edit tests). The fixed root uses the last label, with all labeled orientations tested; it costs no generality. This is finite, not a theorem about larger right cores.

## 2. Exact activity transfer shrinks the 128-link Moore search

Keep the full-chart difference setup and the actual three frozen stars of the previous checkpoint. Each closed branch has the same active size $p$. An outside low branch linked by a changed matching to a closed branch also has exactly $p$ active sheets: each of its active rows moves along the link, so its old image is active, and every active sheet in the closed branch has an active old exterior neighbor. The matching is a bijection. This is stronger than the generic cap 14 or 28 on that outside branch.

For $k=128$ changed links and $h$ high, $s$ closed branches, the high-incidence budget into the closed set is at most $z=4(h+s-44)$. An exterior low branch outside the closed set has all its changed neighbors high. If there are $l$ closed lows, write $d$ for the total number of their changed exits not internal to those $l$ vertices; these go to highs and consume $d$ units of the budget. The finite catalogue of degree-four neighbor quartets implies

    d ≥ D(l) = 0,4,6,6,4,4,2,2   for l=0,…,7.

For $l\le4$ this is the complete-graph degree bound. For $l=5$ there are at most two admissible quartets in a five-set, so at most two internal degree-four vertices, and $d$ is even. For $l=6,7$, zero deficit would mean a 4-regular low graph. Regenerate the 152 quartets and enumerate unions of two different quartets, padded to six or seven vertices. Such a graph has at least two different neighbor sets (take adjacent vertices). Then discard any set where some vertex has no contained quartet avoiding itself. Of 2,373 six-sets and 67,272 seven-sets, zero and 43 survive respectively. On each remaining seven-set the complement of a 4-regular graph is a 7-cycle or a disjoint 3- and 4-cycle. Enumerating these labeled patterns, at most five vertices have admissible neighbor quartets. Thus $d\ne0$ and even parity gives the bound.

At least $\max(0,D(l)-l(s-l))$ of the deficit exits from closed lows go to **outside** highs, and cannot also count as exits of a closed high. All exits from closed highs that do not go to outside lows consume one incidence each, using both incidences for a link internal to the closed highs and the high incidence at the closed-low end for a link from a closed low to a closed high. Thus at least $5(s-l)-z+\max(0,D(l)-l(s-l))$ exits from high closed branches go to outside lows. Dividing this by four and rounding up (or using zero if negative) forces a number $J$ of outside lows to have exact active size $p$. Put $L=54-h-l$ and $B=\max(0,L-\lfloor(55-h)/2\rfloor)$, so at least $B$ outside lows have cap 14 from the hub argument. In the prior integer second moment, use $J$ caps $p$, $\max(0,B-J)$ caps 14, and the rest up to $L$ cap 28. The remaining $57-s-L$ branches have cap 56. This conservatively credits as much overlap as possible between the exact-transfer branches and the cap-14 branches. Enforce $p\le14$ independently whenever $B>0$.

| $h$ | $s$ | $l$ | Necessary $p$ |
|--:|--:|--:|:--|
| 38 | 7 | 7 | 8 |
| 39 | 6 | 5 or 6 | 8,10 |
| 39 | 7 | 6 or 7 | 8 |
| 40 | 4 | 0 | 12,14 |
| 40 | 5 | 4 or 5 | 6,8,10,12,14 |
| 40 | 6 | 5 | 8,10 |
| 40 | 6 | 6 | 6,8,10,12 |
| 40 | 7 | 6 or 7 | 8 |

The three tight rows $(37,7),(38,6),(39,5)$ are now impossible, and the tight $(40,4)$ row allows only 12 or 14. No row establishes that a chart exists. In the other rows most closed branches must be low; their remaining quartets must agree on actual antipodal sheet identities, not just on support size.

Replay: `node transferred_boundaries.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt`. Optional `--order7` adds the longer right-core census.
