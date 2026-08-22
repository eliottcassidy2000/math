---
id: THM-3388
title: "Three-sheet phase triangles and the q=3 cover clutter"
status: >
  PROVED analytic phase-triangle criterion + FINITE-EXACT literal q=3
  clutter/atlas + INDEPENDENTLY HOSTILE-AUDITED.  For three sheets, full
  transverse cover is a 3-uniform clutter.  For speeds u,v,w coprime to
  three, define the finite affine gap set
  P(u,v)={p congruent uv mod 3 gcd(u,v):14|p|<3(u+v)}.  The triple covers iff
  some p,q,r in the three cyclic gap sets satisfy w p+u q+v r=0.  This is the
  vanishing of the normalized phase-gap H^1 class; pairwise overlap alone is
  not sufficient.  In the literal pool there are 48 cover edges inside 82
  pairwise-feasible triangles, an independence profile
  (1,10,45,72,38,6,0,...), and exactly 585 globally safe plus three
  core-rescued q=3 body rows, reproducing the THM-3387 total 588.  A ternary
  dilation orbit gives exact word/support/multiplicity/harmonic recurrences.
  The result classifies the q=3 slice of THM-3387 but proves no new
  refined-ledger decrement or LRC(14).
source: codex-2026-08-14-q3-phase-triangle
audit: independent CRT valuation proof, 506598 gluing instances, 47905 triples through speed 100, exact atlas replay, dilation sign hostile, and harmonic audit
depends_on:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
related:
  - THM-3382-fibonacci-ray-dual-index-harmonic-bifurcation-and-ternary-heap-addresses
  - THM-3366-all-sector-complement-clock-completion
script: 04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py
output: 05-knowledge/results/lrc14_q3_phase_triangle_clutter_thm3388.out
script_sha256: 5323346310a9a6b188caa0131b177b2ae8e23c7113808cda8955f89828e62154
output_sha256: 5a32319fb8a91b476d292da292ae3cc9933f5f94aad7eb0e834f49e52252c535
semantic_sha256: 082e97aa25d8019ba7de49c0a76333c7a3a221dd19cb4bc3e8d5b43ef9a42216
hash_basis: LF-normalized bytes
---

# THM-3388 -- three-sheet cover is phase closure, not a triangle of pair tests

**PROVED analytic phase-triangle criterion + FINITE-EXACT literal `q=3`
clutter/atlas + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and the missing coordinate

THM-3387 proves that the exact cyclic obstruction is the locus where the
transverse speeds cover every sheet.  Its `q=2` obstruction is binary and
collapses to a gcd graph.  At `q=3`, every transverse speed is coprime to
three and still blocks at most one sheet, so the first possible obstruction
uses three speeds.  The closest corrected near miss is MISTAKE-382: deleting
a grid can lose isolated pointwise information.  The least-used sidecar is
the signed phase gap between the three sheet owners.

| field | connection |
|---|---|
| source | three transverse speeds on a degree-three fibre |
| target | an affine integral 1-cochain on the cycle of three sheet owners |
| map | take the three signed centre gaps and clear their denominators |
| preserved | strict overlap, cyclic sheet assignment, common-phase existence |
| destroyed by pair tests | gap magnitude, sign, and circulation around the cycle |
| sidecar | one value in each affine gap set plus its weighted closure equation |
| decisive hostile | `(1,4,7)`: every pair can overlap on different sheets, but no phase triangle closes |

## 2. The affine gap sets

Let `u,v,w` be distinct positive integers not divisible by three.  Write a
source point in the three sheets as

```text
t,                  t+1/3,                  t+2/3.       (1)
```

A speed blocks at most one of them, because THM-3385's capacity is
`ceil(3/7)=1`.  Assign `u,v,w` cyclically to the sheets in `(1)`.  The centres
of their strict danger intervals can be lifted to

```text
alpha=a/u,          beta=b/v-1/3,           gamma=c/w-2/3, (2)
```

with radii `1/(14u),1/(14v),1/(14w)`.  Define the cleared oriented gaps

```text
p=3uv(alpha-beta)=3(va-ub)+uv,
q=3vw(beta-gamma)=3(wb-vc)+vw,
r=3wu(gamma-alpha)=3(uc-wa)-2wu.             (3)
```

For `g_uv=gcd(u,v)`, the first line and its cyclic companions imply

```text
p == uv (mod 3g_uv),
q == vw (mod 3g_vw),
r == wu (mod 3g_wu).                         (4)
```

The last congruence uses `-2 == 1 (mod 3)`.  The two intervals for `u,v`
overlap strictly exactly when

```text
14|p|<3(u+v).                                 (5)
```

Accordingly define the finite affine set

```text
P(u,v)={p in Z:p == uv (mod 3 gcd(u,v)), 14|p|<3(u+v)}. (6)
```

It is nonempty exactly when

```text
3(u+v)>14 gcd(u,v).                           (7)
```

This is the sharp **pairwise** different-sheet overlap test.  It is only a
necessary shadow of three-sheet cover.

## 3. Exact phase-triangle criterion

The three differences in `(3)` telescope.  After clearing the common
denominator, their closure law is

```text
w p+u q+v r=0.                                (8)
```

Conversely suppose

```text
p in P(u,v), q in P(v,w), r in P(w,u),
w p+u q+v r=0.                                (9)
```

Put

```text
A=(p-uv)/3, B=(q-vw)/3, C=(r+2wu)/3.          (10)
```

The congruences `(4)` say

```text
gcd(u,v)|A, gcd(v,w)|B, gcd(w,u)|C,            (11)
```

and `(8)` becomes `wA+uB+vC=0`.  The elementary cycle-gluing lemma says
that `(11)` and this relation are equivalent to the existence of integers
`a,b,c` with

```text
va-ub=A, wb-vc=B, uc-wa=C.                    (12)
```

For completeness, solve the first two equations as congruences for `b`
modulo `v/gcd(u,v)` and `v/gcd(v,w)`.  Their generalized-CRT compatibility,
prime by prime, is exactly `gcd(u,w)|C`; after choosing `b`, the first two
equations give `a,c`, and `wA+uB+vC=0` forces the third.  Thus `(9)` really
glues one triple of centres; it is not three independently chosen pair
witnesses.

The bounds in `(6)` make the three selected single-tooth arcs pairwise
intersect.  Three pairwise-intersecting circular arcs with no common point
must cover the circle, but these arcs have total length at most `3/7<1`.
Hence they have a common source time `t`, proving a three-sheet cover.  The
reverse implication was derived in `(2)`--`(8)`.  Therefore

```text
{u,v,w} covers all three sheets somewhere
iff there exist p,q,r satisfying (9).          (13)
```

A cyclic relabelling of the sheets rotates the construction.  Reversing the
sheet order is induced by `t -> -t`, so `(13)` is independent of the chosen
ordering of the three speeds.

Normalize `(3)` back to the rational gaps

```text
delta_uv=p/(3uv), delta_vw=q/(3vw), delta_wu=r/(3wu).
```

Equation `(8)` is `delta_uv+delta_vw+delta_wu=0`.  Thus the allowed affine
edge cochain has zero class in `H^1(C_3;Q)`, or equivalently is the coboundary
of the centre `0`-cochain `(alpha,beta,gamma)`.  This is the precise H1
content: pair feasibility populates three edge fibres; phase closure decides
whether they glue.

## 4. The cover clutter and the literal atlas

For an arbitrary transverse set `U`, a full cover selects one speed blocking
each sheet.  Since one speed cannot block two sheets,

```text
B_3(U)=empty iff U contains no triple satisfying (13). (14)
```

So the intrinsic object is a 3-uniform cover clutter, and globally safe sets
are its independent sets.  On

```text
V={1,2,4,5,7,8,10,11,13,14},                 (15)
```

the exact companion finds `48` cover edges.  The pair graph from `(7)` has
`82` triangles, leaving `34` pairwise false positives.  The smallest clean
hostile and positive control are

```text
(1,4,7): P(1,4),P(4,7),P(7,1) all nonempty, no (8);
(1,4,5): (p,q,r)=(1,-1,-1) satisfies (8).     (16)
```

The independence polynomial is

```text
I_V(z)=1+10z+45z^2+72z^3+38z^4+6z^5.         (17)
```

The six independent five-sets are frozen in the exact output.  For literal
six-speed bodies, choose `c` core clocks from `{1,2,3,4}` and `6-c`
transverse vertices from `V`.  Equation `(17)` gives

```text
4 I_5+6 I_4+4 I_3+I_2
=4*6+6*38+4*72+45=585                         (18)
```

globally transverse-safe rows.  An independent exact integer event sweep of all
`2,793` candidate body rows finds exactly three additional core rescues:

```text
C=(1,2,3),(1,2,4),(2,3,4), always U=(8,11,13). (19)
```

Hence the q=3 slice has exactly

```text
585+3=588                                     (20)
```

pointwise-exact rows, reproducing the q=3 entry of THM-3387 by a structural
classification rather than a black-box interval total.

## 5. Ternary ancestry, subsets, and harmonic weight

Common dilation by any positive integer `s` coprime to three preserves `(13)`.
If `s==1 (mod 3)`, send every gap `p` to `sp`; if `s==-1 (mod 3)`, send it
to `-sp`.  The affine congruences and strict bounds scale correctly, and
`(8)` is multiplied by `+s^2` or `-s^2`.

Start with the edge `E={1,4,5}` and the three multipliers `7,11,13`.  At word
depth `d` there are `3^d` ancestry words, but commutation leaves only
`binom(d+2,2)` scale values `7^i 11^j 13^k`, `i+j+k=d`, with multinomial
collision multiplicities.  The three root orbits are disjoint because their
`2`- and `5`-adic valuations differ.  Every resulting three-element block is
again a cover edge.

Let `H_d` be the harmonic mass of the **distinct integer support** at depth
`d`, and `W_d` the ancestry-word multiplicity-weighted mass.  Since

```text
sum_(e in E)1/e=29/20,
1/7+1/11+1/13=311/1001,
1/(7*11)+1/(7*13)+1/(11*13)=31/1001,          (21)
```

the complete-homogeneous recurrence and multinomial theorem give

```text
1001 H_d=311 H_(d-1)-31 H_(d-2)+H_(d-3),       d>=3,
1001 W_(d+1)=311 W_d,                           d>=0.    (22)
```

The whole integer orbit is a convergent subseries of the harmonic series:

```text
sum_(n in orbit)1/n
=(29/20)(1-1/7)^(-1)(1-1/11)^(-1)(1-1/13)^(-1)
=29029/14400.                                  (23)
```

More generally, for every subset `A subset N^3`, take exactly the blocks
indexed by `A`:

```text
S_A=union_((i,j,k) in A){7^i11^j13^k,4*7^i11^j13^k,5*7^i11^j13^k}.
```

Then `A -> S_A` is a Boolean-algebra embedding onto the block-saturated
subalgebra of the subsets of this integer orbit (complements are relative to
the orbit), and

```text
sum_(n in S_A)1/n=(29/20)sum_((i,j,k) in A)7^(-i)11^(-j)13^(-k). (24)
```

This is an exact structured realization of subsets as harmonic subseries.
It does not identify ancestry words with integer support: word order is the
sidecar destroyed by the commutative monoid quotient.

## 6. Why this is not a tournament

The pair observable `(7)` is symmetric.  Even retaining all three pair
answers loses the circulation, as `(1,4,7)` proves.  Orienting those answers
would add a gauge but not the missing integer gaps.  The tempting size-four
packet--three speed intervals plus a common time--is an incidence witness,
not a tournament: eliminating the existential time produces a ternary
hyperedge with the affine cochain as sidecar.  Missing edges and simultaneous
edges must remain literal.

## 7. Exact verification and scope

The standard-library companion:

- checks `(7)` on `8,778` pairs below `200`;
- compares `(13)` with an independent exact endpoint/mid-cell event sweep on
  all `2,925` triples through speed `40`;
- reconstructs integral centres for all `2,209` positive triples and checks
  all six ordering gauges in the low range;
- enumerates the literal clutter, independence polynomial, all `2,793` q=3
  body rows, and the three core rescues; and
- checks `21` ternary lattice shells, both recurrences, collision counts, and
  the exact total harmonic mass.

The independent hostile audit additionally checked `506,598` admissible
integer gluing instances and every one of `47,905` triples through speed
`100`, including all six orderings.  It found no discrepancy.  Three controls
from that audit are now frozen in the companion:

```text
(1,4,41): every pair gap set is nonempty, but no phase triangle closes;
2*(1,4,5): the lawful scaled gaps are (-2,2,2), not (2,-2,-2);
C=(1,3,4), U=(8,11,13), t=389/2464:
  all sheets fire and only the omitted core clock 2 is dangerous.          (25)
```

It contains no floating literal or optimization-dependent `assert`.  Reproduce
with

```text
python 04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py
python -O 04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py
```

Ordinary and optimized runs LF-normalized-byte-match the stored output.

This theorem classifies the q=3 slice inside THM-3387.  It gives no new
refined-ledger subtraction, physical drift realization, arbitrary-phase
transport, or proof of LRC(14).

**QED.**
