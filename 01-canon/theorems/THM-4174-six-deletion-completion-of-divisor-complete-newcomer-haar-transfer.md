---
id: THM-4174
title: "Six-deletion completion of divisor-complete newcomer Haar transfer"
status: >
  PROVED RELATIVE TO THM-4150/4156/4158/4170 + VERIFIED-EXACT STANDARD-LIBRARY
  WALL-LATTICE/SUBMASK CENSUS + SECOND-PATH NUMPY GROUP-LATTICE/MILP AND
  EXACT-TRANSVERSAL AUDITS; LRC(14) OPEN. Every positive newcomer q outside
  the fixed THM-4156 pool closes all 888,030 anchored seven-choice bodies at
  every common content with every pair of distinct positive odd tails. The
  exact deletion-arity failure filtration has sizes 61,9,1,0 at d=3,4,5,6.
source: codex-lrc-triple-deletion-synthesis-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4158-three-band-wrapped-carrier-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
related:
  - THM-4160-anchored-haar-deletion-cover-and-content-tower
  - THM-4166-two-deletion-repair-graph-haar-odd-tail-transfer
script: 04-computation/lrc14_multideletion_all_newcomers_exact_thm4174.py
output: 05-knowledge/results/lrc14_multideletion_all_newcomers_exact_thm4174.out
independent_four_script: 04-computation/lrc14_four_deletion_residual_audit_thm4174.py
independent_four_output: 05-knowledge/results/lrc14_four_deletion_residual_audit_thm4174.out
independent_five_script: 04-computation/lrc14_five_deletion_residual_audit_thm4174.py
independent_five_output: 05-knowledge/results/lrc14_five_deletion_residual_audit_thm4174.out
independent_six_script: 04-computation/lrc14_q50_six_deletion_arity_audit_thm4174.py
independent_six_output: 05-knowledge/results/lrc14_q50_six_deletion_arity_audit_thm4174.out
script_sha256: b8ad9ab60715ddde853a45459ccbc6d01f1d6e2c6a31f48eaa497765208722de
output_sha256: 8679f01a234987cc7784419d1013cac9c26a16f34ac2e00d4d94e6ee09327586
independent_four_script_sha256: 6185e541a234d7b0a6f672b44ed0fdbd0756db98368f1715e70ffcd850a5948f
independent_four_output_sha256: a7dacffa7ab64a497df3dab511249134d20bf5607561b7f968fe77b40c1beaec
independent_five_script_sha256: d018b7a84c3683cb3b1b320d0e17dbd6fba8dfe370752e7fee6d149e5e0b6b74
independent_five_output_sha256: 7246a6644c35b9a1afca039628ef7810dd51523535098344126c166d466a54e4
independent_six_script_sha256: 4536216aa8435fac609f9f5ed22541326fc0931a366866756b0e669e00a9b2cb
independent_six_output_sha256: ee13bbc4bb7fbdc51b8294635280d7bee6d501fe080f2815ac637d087edbd62b
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standard-library Fraction/submask implementation reconstructs all
  7,134 walls, performs 2,093,130 exact threshold comparisons, verifies zero
  equalities, and decides every residual transversal through depth seven.
  Normal, optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT WITH DECLARED SHARED GEOMETRY. Separate NumPy group-lattice arrays,
  SciPy MILP witnesses, and separately instantiated exact-recursion checks
  reproduce the d=4 and d=5 residual ledgers and the q=50 d=6 tau=8 control.
  The group-lattice/MILP route is algorithmically distinct; the recursions
  share branch order and all paths remain Python implementations sharing the
  exact wall-lattice identity. No independent C++ geometry is claimed.
---

# THM-4174 -- six-deletion completion of every newcomer

**PROVED RELATIVE TO THM-4150/4156/4158/4170 + VERIFIED-EXACT; LRC(14) REMAINS
OPEN.**

## 1. Statement

Retain the THM-4156 anchors and pool

```text
A={120,126,143},

P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},

O=P\A.                                                   (1)
```

For every positive integer `q notin P` and every `0<=d<=27`, define the
`d`-deletion repair hypergraph on `O` by

```text
E_d(q)={R in binom(O,d):
        mu(G_((P union {q})\R))>=4/63},                 (2)

G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                (3)
```

Let `B_3` be the exact 61-label failure set proved in THM-4170:

```text
B_3={3,6,22,24,25,46,48,50,55,64,70,72,75,83,93,96,
     100,103,105,110,122,127,128,140,147,153,158,166,
     172,173,183,186,192,206,210,220,256,260,270,282,
     294,306,313,320,325,332,346,366,372,384,416,440,
     462,512,516,520,550,567,744,768,924}.              (4)
```

Put

```text
B_4={25,50,96,100,105,128,192,210,256},
B_5={50},
B_6=empty.                                               (5)
```

The inheritance pass is exact. The closest mechanism is
[THM-4170](THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer.md),
whose 61-label set is the entire input universe here. The canonical hostile
is `q=50`, where five deletions produce 2,063 repair edges but still admit a
five-vertex transversal. The corrected near miss is therefore “many repair
edges imply certification”; edge count is not the invariant. The least-used
sidecar is deletion arity itself, tracked through the nested exact failure
filtration `B_3 superset B_4 superset B_5 superset B_6`.

> **Deletion-filtration theorem.** For every positive `q notin P`,
>
> ```text
> tau(E_3(q))>7 iff q notin B_3,
> tau(E_4(q))>7 iff q notin B_4,
> tau(E_5(q))>7 iff q notin B_5,
> tau(E_6(q))>7.                                        (6)
> ```
>
> Consequently, for every `K in binom(O,7)`, every positive integer `c`,
> and every two distinct positive odd integers `a,b`, there is an
> `x in R/Z` such that, with
>
> ```text
> H_(q,K)=A union {q} union K,
> ```
>
> one has
>
> ```text
> min_(v in 2cH_(q,K) union {a,b})||vx||>=1/14.          (7)
> ```

Thus every newcomer `q notin P` supplies exactly

```text
binom(27,7)=888,030                                     (8)
```

primitive divisor-complete core bodies. Relative specifically to THM-4170's
triple-deletion certificate, THM-4174 supplies the `61*888,030 = 54,169,830`
bodies at its exact failure labels. This is not a global novelty count against
every other canonical mechanism. Their arity split is

```text
d=4: 52*888,030=46,177,560,
d=5:  8*888,030= 7,104,240,
d=6:  1*888,030=   888,030.                             (9)
```

## 2. Why a transversal closes every body

Fix `q`, `d`, and `K in binom(O,7)`. If `tau(E_d(q))>7`, then `K` is not a
transversal, so some repair edge `R in E_d(q)` is disjoint from `K`. Hence

```text
A union {q} union K subset (P union {q})\R,
G_((P union {q})\R) subset G_(A union {q} union K).     (10)
```

The right-hand safe set has Haar mass at least `4/63`. THM-4150 therefore
closes every distinct positive odd tail pair after doubling.

The predicate is monotone upward in deletion arity. If every seven-set misses
a `d`-repair and `d<d'<=20`, extend that repair inside the 20-element
complement `O\K`; deleting more labels only enlarges the safe set. Thus every
qualifier at arity `d` also qualifies at arity `d'`. This proves all directions
of `(6)` outside the displayed residual sets once the exact finite results
below are established.

## 3. Exact wall-lattice computation

The fixed pool has exactly

```text
7,134 rational walls, 7,133 open cells,
L=lcm{14p:p in P}=18,241,159,416,480.                   (11)
```

At the midpoint of each cell, record `F subset O`, the exact set of failed
optional labels, and discard the cell if an anchor fails. The complete cell
histogram by `|F|` is

```text
0:150, 1:328, 2:518, 3:678, 4:728, 5:666, 6:472,
7:242, 8:102, 9:38, 10:20, 11:6, 12:2, 13:2, 14:2,
15:2, anchor-failing:3,177.                             (12)
```

Write a wall as `t/L`, and put `qt=kL+r`, `0<=r<L`. The numerator of the
safe-comb prefix integral through that wall, on denominator `14qL`, is

```text
J_q(t)=12kL+
  0,                 if 14r<=L,
  14r-L,             if L<14r<13L,
  12L,               if 13L<=14r.                      (13)
```

For a deletion mask `R`, sum the consecutive prefix differences exactly over
the cells whose failure mask satisfies `F subset R`. If the resulting integer
numerator is `N(q,R)`, then

```text
mu(G_((P union {q})\R))>=4/63
 iff 9N(q,R)>=8qL.                                     (14)
```

The new audit performs `61*binom(27,4)+9*binom(27,5)+binom(27,6)
=2,093,130` comparisons. There are no threshold equalities.

Each resulting hypergraph is passed to an exhaustive depth-seven transversal
recursion. At every search state it chooses an uncovered edge and branches on
its vertices; a greedy disjoint-edge packing gives a valid lower-bound prune.
All returned covers are checked against every active edge. A second path uses
global failure-group arrays and exact signed-integer vectorization, then checks
the same hypergraphs both by MILP witnesses and by a separately instantiated
exact-recursion cross-check.

## 4. The exact residual filtration

At arity four, exactly the nine labels in `B_4` retain a cover of size at
most seven. Their exact minimum covers are:

|`q`|repair 4-edges|`tau`|one minimum cover|
|---:|---:|---:|:---|
|25|925|6|`{85,88,145,168,193,252}`|
|50|4|1|`{193}`|
|96|421|6|`{8,95,145,168,193,252}`|
|100|538|5|`{85,193,240,252,290}`|
|105|1,328|6|`{80,85,88,145,168,252}`|
|128|819|6|`{16,80,85,88,168,252}`|
|192|1,557|7|`{16,88,95,168,193,252,290}`|
|210|891|5|`{85,145,168,193,252}`|
|256|520|6|`{16,80,88,168,193,252}`|

The other 52 labels of `B_3` have no cover of size at most seven and are
therefore rescued at arity four.

At arity five, eight of the nine labels in `B_4` have no seven-cover. The
sole residual is

```text
q=50: 2,063 repair 5-edges, tau=5,
cover {95,168,193,240,290}.                             (15)
```

At arity six, the same hostile becomes positive:

```text
q=50: 46,261 repair 6-edges, tau=8,
minimum cover {16,85,88,95,168,193,240,290}.            (16)
```

The exhaustive recursion rejects every cover through size seven and verifies
the displayed size-eight cover. Thus `B_6` is empty. The closest positive
six-edge margin at `q=50` is

```text
46099/15200966180400                                   (17)
```

above `4/63`; the closest failure lies

```text
30091/1200076277400                                    (18)
```

below it.

The failure anatomy at `q=50` is especially transparent. Its triple repair
hypergraph is empty. Its four-deletion hypergraph consists of only

```text
(88,95,176,193), (88,145,176,193),
(88,145,193,290), (145,168,193,290),                    (19)
```

all hit by the single label `193`. Five deletions disperse the support enough
to raise `tau` to five, and six deletions cross the exact target `tau=8`.
Edge count alone is not the invariant; the decisive coordinate is transversal
dispersion across deletion arity.

## 5. Haar transfer, content, and uniqueness

For every `q notin P`, equations `(4)--(6)` give some `d<=6` with
`tau(E_d(q))>7`; Section 2 supplies `(10)`, and THM-4150 proves `(7)` at
content one. For positive `c`, the surjective circle endomorphism

```text
m_c(y)=cy
```

satisfies

```text
G_(cH)=m_c^(-1)(G_H),             mu(G_(cH))=mu(G_H),   (20)
```

so the same transfer proves every content.

Every core is primitive because it contains `A` and
`gcd(120,126,143)=1`. The anchors contain a multiple of every modulus
`2,...,14`, as in THM-4156. Since `q notin P`, each core has exactly one
outside-pool label. Equality of two cores first identifies `q` and then `K`.
Across contents, `gcd(2cH)=2c`, so equality of physical body sets first
identifies `c` and then the core. Body labels are even and tails are odd, so
no body--tail collision occurs. This proves `(7)--(9)`. **QED.**

## 6. Scope, overlaps, and a sharp wrapped-carrier separation

The restriction `q notin P` is type-critical. If `q in O\K`, then

```text
A union {q} union K=A union (K union {q})
```

is an old THM-4156 body, and every old eight-choice body has eight such
presentations. If `q in K` or `q in A`, the union has only ten labels. Thus
including `q in P` would either duplicate the old family or leave the
eleven-body target; it is not an unclassified newcomer gap.

For `N>=290`, the scale-one newcomer classes with `1<=q<=N` have exact count

```text
(N-30)binom(27,7),                                      (21)
```

and are disjoint from the `binom(27,8)=2,220,075` old THM-4156 cores.

There is also a sharp content-stable separation from the full THM-4158
wrapped-carrier theorem. If `cH_(q,K) subset P_m`, then every member of `P_m`
is at least `m`, while `120c in cH`; hence `m<=120c`. The maximum member of
`P_m` is

```text
floor(41(12m+1)/16)<=3690c+2.                           (22)
```

For `q>=3693`, the label `cq>=3693c>3690c+2`, a contradiction. Therefore
every THM-4174 body with `q>=3693`, at every content, lies outside every
three-band wrapped alphabet `P_m`. The boundary is sharp at content one:
with `m=120`,

```text
q=3692,
K={132,145,168,170,176,190,193}
```

gives `H_(q,K) subset P_120`, since the anchors and `K` lie in the first band
and `3692` is the last point of the third band.

This theorem closes one fixed anchored newcomer lane. It does not classify
arbitrary divisor-complete bodies, treat mixed/even tail parity, make the Haar
threshold necessary, prove any body unsafe, or prove LRC(14).

## 7. Source-to-target contract and replay

```text
source:       d=4,5,6 Haar repair hypergraphs over the THM-4156 pool
target:       every q-outside-P anchored seven-choice body and odd tail pair
map:          no seven-cover -> missed repair -> safe inclusion -> THM-4150
preserved:    q label, anchors, exact 4/63 threshold, content, divisor pins
destroyed:    component addresses, selected repair, body-specific safe phases
sidecar:      exact arity filtration B3 superset B4 superset B5 superset B6
positive:     q=50,d=6 has 46,261 edges and tau=8
hostile:      q=50,d=5 has 2,063 edges but a five-cover
decisive test: 2,093,130 exact comparisons plus exhaustive transversals. (23)
```

Primary replay:

```bash
python3 -B 04-computation/lrc14_multideletion_all_newcomers_exact_thm4174.py
python3 -B -O 04-computation/lrc14_multideletion_all_newcomers_exact_thm4174.py
PYTHONHASHSEED=4174 python3 -B \
  04-computation/lrc14_multideletion_all_newcomers_exact_thm4174.py
```

Replay the independent algorithmic audits with

```bash
python3 -B 04-computation/lrc14_four_deletion_residual_audit_thm4174.py
python3 -B 04-computation/lrc14_five_deletion_residual_audit_thm4174.py
python3 -B 04-computation/lrc14_q50_six_deletion_arity_audit_thm4174.py
```

All four fresh streams byte-match the frozen outputs. The second path uses
different group arrays and optimization witnesses; its exact-recursion check
shares the primary branch order. All paths share the rational wall-lattice
identity and Python runtime.
A warning-clean signed-128-bit C++ reconstruction would be a valuable third
path; it is not represented as an existing audit or a dependency of the
mathematical double-counting proof.
