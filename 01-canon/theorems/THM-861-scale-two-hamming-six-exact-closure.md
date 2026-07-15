---
id: THM-861
title: Scale-two Hamming-six exact closure
status: PROVED + FINITE-EXACT — the c=2 primitive AP-centred H6 chart has exactly one tight packet, the ordinary AP [12]; every sporadic packet is loose; 41,882,982-node logical census with independent closed-danger replay
source: codex-2026-07-15-S14/S15 continuation
depends_on: [THM-815, THM-857, THM-860]
related: [THM-823, THM-858, THM-859, THM-862, HYP-6820]
verification:
  - 04-computation/lrc13_scale_two_hamming_six_component_recursion_codex_S14.cpp
  - 05-knowledge/results/lrc13_scale_two_hamming_six_component_recursion_codex_S14.out
  - 04-computation/lrc13_scale_two_hamming_six_closed_danger_replay_codex_S14.cpp
  - 05-knowledge/results/lrc13_scale_two_hamming_six_closed_danger_replay_codex_S14.out
---

# THM-861 — the scale-two Hamming-six chart closes

Use THM-860's notation.  Thus `R subset [12]` has six elements and

```text
A=2([12] minus R) union {w_r:r in R},
w_r=2r (mod 13),              w_r>0,       w_r!=2r.      (1)
```

Assume `A` is primitive and `M(A)<=1/13`.  THM-860 proves that at common
scale two all six effective orders are two and that `R` is one of exactly 64
signed-doubling common-sheet cycles.  There is only one unit modulo two.  For
each label the complete legal replacement ray is therefore

```text
w_r=u_r+26k,       k>=0,
u_r=2r+13  (1<=r<=6),        u_r=2r-13  (7<=r<=12).     (2)
```

## Theorem

Across all 64 roots and all unbounded heights in (2), there is exactly one
terminal whose strict-safe set is empty:

```text
R={7,8,9,10,11,12},
(r:w_r)=(7:1,8:3,9:5,10:7,11:9,12:11).                (3)
```

Its retained speeds are `2{1,...,6}={2,4,...,12}` and its replacements are
`{1,3,...,11}`.  Their union is exactly `[12]`, so (3) is the ordinary AP
equality and `M(A)=1/13`.  Every other packet (1) has

```text
M(A)>1/13.                                              (4)
```

Consequently the **sporadic** part of the primitive AP-centred Hamming-six
branch at `c=2` is empty.  The word “sporadic” is essential: the literal
chart is not empty, because the AP itself has the proper scale-two
presentation (3).

## 1. Why the 64 roots and the rays are complete

At `c=2`, an effective order is one or two.  THM-860(D) shows that a primitive
common-sheet presentation cannot contain an order-one colour: all six orders
are two.  At owner `o`, its own colour fills one of the two sheets.  A distinct
provider `r` fills the opposite sheet exactly when

```text
r -> o       iff       o/r in {+2,-2} in F_13.          (5)
```

Positive indegree at every owner forces a directed cycle.  No cycle of length
two through five can close because `(+-2)^m!=1 mod 13` there.  Hence every
survivor is one directed Hamiltonian six-cycle.  Since `2^6=-1 mod 13`, it
closes exactly for an odd number of negative edges; a chord would create a
shorter forbidden cycle.  Thus the number of roots is

```text
12*2^5/6=64.                                            (6)
```

Solving the two CRT conditions `u=2r mod 13` and `u=1 mod 2` gives (2), and
adding the common modulus `lcm(13,2)=26` gives every positive replacement.
Therefore the root list and every labelled future ray used below are
exhaustive, including the downward replacements in (3).

## 2. Complete exact component recursion

For a speed set `Q`, write

```text
E(Q)={t in R/Z: ||qt||>1/13 for every q in Q}.           (7)
```

Numerically order the six replacements `x_1<...<x_6`.  At depth `j`, retain
the exact sorted rational components of the prefix set `E_j`, its remaining
labelled step-26 rays, and `x_j`.  With `L_j` the longest component length and
`s=6-j`, THM-815 gives every hypothetical covering continuation the cap

```text
x_(j+1)<=floor(22s/[13(13-2s)L_j]).                     (8)
```

For `j<6` the prefix has at most eleven speeds, so the lower-runner theorem
gives `L_j>0`.  Equation (8) therefore enumerates only finitely many members
of each ray.  Distinct labels have distinct residues modulo thirteen, so each
terminal occurs in one numerical order.  Induction proves completeness with
no height cutoff.

The two exact shortcuts are the ones proved in THM-857.  If a candidate's
full open safe tooth lies inside a parent component, the child's retained
length makes (8) smaller than the candidate, so there is no ordered tight
continuation.  Otherwise the intersection is streamed in endpoint order; as
soon as one child interval makes the next cap strictly smaller than the least
legal future ray member, that child is dead.  At depth six, the first nonempty
child interval certifies looseness.  Every comparison is an integer cross
product; neither shortcut changes the logical node or edge census.

## 3. Frozen census

The complete cap-admissible logical tree is:

| depth | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| nodes | 64 | 6,266 | 641,866 | 40,800,088 | 433,938 | 758 | 2 |
| dead / no possible tight continuation | 0 | 0 | 0 | 40,666,446 | 433,354 | 756 | 0 |
| complete-tooth certificates | 0 | 0 | 0 | 38,293,571 | 0 | 0 | 0 |
| streaming-cap certificates | 0 | 0 | 0 | 2,372,875 | 433,354 | 756 | 1 |

Thus

```text
logical nodes               41,882,982
candidate edges             41,882,918
covering terminals                   1
loose terminals                      1.                (9)
```

The non-shortcut materialized prefixes have these exact extrema:

| depth | prefixes | minimum longest component | multiplicity | maximum next cap |
|---:|---:|---:|---:|---:|
| 0 | 64 | `11/780` | 1 | 720 |
| 1 | 6,266 | `11/9347` | 1 | 2,396 |
| 2 | 641,866 | `11/31109` | 1 | 3,828 |
| 3 | 133,642 | `2779/3780205` | 1 | 986 |
| 4 | 584 | `411/165529` | 1 | 151 |
| 5 | 2 | `1/156` | 1 | 24 |

In particular the actual `c=2` survivor atlas sharpens THM-860's generic
first-ray bound: its root cap is at most 720, hence each ray has at most 28
members before later prefix caps are applied.

## 4. Independent closed-danger replay

The companion implementation does not reuse the primary interval-intersection
routine.  At every expanded prefix it constructs all closed danger teeth from
scratch, sorts and merges touching intervals, takes their open complementary
gaps, and independently subtracts the candidate's closed comb.  It compares
that result endpoint-for-endpoint with a full child reconstruction.  It also
replays the exact parent gap, tooth, child gap, remaining count, least-future
speed, and cap behind every shortcut.

Each root emits a canonical SHA-256 certificate over its labelled nodes and
witnesses; the root certificates are combined in root-index order.  The
all-root replay asserts both tables, (3), and (9).  Frozen source, output, and
manifest digests are recorded with the reproduction commands below.

## 5. Sparse cycle telemetry and the faithful object

For (5), the tournament vertices are challenged before any completion: the
six missing labels form a sparse signed relation, not a tournament.  Every
root has

```text
edges=6, directed cycles=1, SCC size=6, Hamiltonian paths=6,
negative-edge histogram over roots={1:12,3:40,5:12}.   (10)
```

The sign switch exchanges the two sheet routes.  Completing absent pairs by
an increasing-label tie path creates a tournament, but destroys the sheet
predicate.  Even the exact signed cycle decides only common-sheet routing: it
does not distinguish the AP terminal (3) from the 63 loose roots or decide
heights inside a root.  The predicate-preserving object is instead

```text
(literal strict-safe components,
 remaining labelled step-26 rays,
 last numerical speed,
 exact shortcut ancestry).                            (11)
```

This is a useful Kakeya-comb lesson: the six-cycle is the arithmetic stalk,
while metric coverage lives in the recursively eroded interval fibre over it.

## Scope guardrail

This theorem closes the primitive AP-centred Hamming-six chart at common scale
`c=2`, including all unbounded ray heights.  It does **not** evaluate the
`c>=3` ramified presentation bank from THM-860.  In particular THM-862's
scale-three bank has `336`, `672`, and `496` unit contexts.  This theorem also
does not evaluate the finite Hamming-five bank from
THM-858, the seven-comb wall, or global `n=12` sporadic-branch emptiness.

Reproduce the primary census with

```bash
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_two_hamming_six_component_recursion_codex_S14.cpp \
  -o /tmp/thm861-primary
/tmp/thm861-primary > /tmp/thm861-primary.out
cmp /tmp/thm861-primary.out \
  05-knowledge/results/lrc13_scale_two_hamming_six_component_recursion_codex_S14.out
```

Reproduce the independent replay with

```bash
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_two_hamming_six_closed_danger_replay_codex_S14.cpp \
  -o /tmp/thm861-replay
/tmp/thm861-replay > /tmp/thm861-replay.out
cmp /tmp/thm861-replay.out \
  05-knowledge/results/lrc13_scale_two_hamming_six_closed_danger_replay_codex_S14.out
```

The frozen SHA-256 values are:

```text
primary source   2241f57f1c244928dc35e4eecc54c63a8e45f201a6ae66aea1efe304caf616b3
primary output   bd7b25d352d499ff36a52cad2cf36cd77245e6dbbcac383db6d70d20d88bbe28
replay source    6faf0e5dac31c47a5cd60389632b6f7c002f5bbc75c773b97cf8ed02d90fb137
replay output    95c46bd553eff29a49a895370f02da5b77fb8212eac82957f069e6983856f735
replay manifest  d963e886d84c60414ac9eb51c84cde2dc06d8ee35a01758c598071ff5a3d75e0
```
