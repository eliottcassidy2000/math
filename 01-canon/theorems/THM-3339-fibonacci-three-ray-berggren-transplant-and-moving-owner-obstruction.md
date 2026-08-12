---
id: THM-3339
title: "Fibonacci three-ray Berggren transplant, Cassini edge-current, and affine moving-owner obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY STRUCTURAL-AUDITED.  The positive
  solutions of |n^2-mn-m^2|=1 are exactly consecutive Fibonacci parameters.
  After the mandatory odd/odd primitive normalization, their Pythagorean
  triples occupy three interleaved exact Berggren rays with parameter words
  (BA)^r, A(BA)^r, and C(BC)^r.  The associated four-entry window carries the
  Pythagorean triple as an asymmetric current on the three perfect matchings
  of K4.  Oriented consecutive Farey flanks give the six total orders of the
  three matching channels, with Cassini sign (-1)^k.  The six edge products
  form a transitive T6 whose only changing arc is the opposite-edge pair
  03/12; that isolated swap is odd and cannot be induced by any K4 relabeling.
  Each adjacent matching reflection has four S4 lifts: among the 16 lift
  pairs, exactly four generate an owner-fixing S3 (one for each owner), while
  twelve generate transitive S4.  Thus the quotient order does not choose an
  affine owner.  No LRC, Jacobian, global Berggren--Farey, or tournament
  equivalence follows.
source: codex-2026-08-12-fibonacci-berggren-branch-transplant
audit: >
  independent structural audit (exact branch identities, period-three
  primitive normalization, residue hexagon, Cassini edge-order hostile,
  S4/V4 owner-lift census, and tournament typing: ACCEPT)
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
related:
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
script: 04-computation/fibonacci_berggren_three_ray_owner_thm3339.py
output: 05-knowledge/results/fibonacci_berggren_three_ray_owner_thm3339.out
script_sha256: 1cef777f9920196f37f33d6208e31b53f2ce1b4aa8980df74ef8c9b446f478cf
output_sha256: 355324c33fcb1b87bc1265776b24049be42cbf5e6b2bc9bb53443a08e40cf279
hash_basis: working-tree bytes (LF)
---

# THM-3339 -- the Fibonacci locus is three Berggren rays

**PROVED + VERIFIED-EXACT + INDEPENDENTLY STRUCTURAL-AUDITED.**

This is a repository synthesis and proof interface; no literature-priority
claim is made.  The theorem keeps four superficially similar objects apart:

1. one internal `K4` whose vertices are four entries of a Fibonacci window;
2. the six edge products of that `K4`;
3. the three perfect matchings, viewed as `V4` directions; and
4. the six residue states, which are six orders of those three directions.

Only item 2 has a natural weighted tournament on six vertices.  Item 4 is a
hexagon of transitive tournaments on three channel vertices, not one
tournament on six vertices.

## 1. The generalized window carries one Pythagorean current

For positive integers `0<m<n`, put

```text
W(m,n)=(e_0,e_1,e_2,e_3)=(n-m,m,n,n+m),                  (1)
w_ij=e_i e_j.                                            (2)
```

The six products satisfy

```text
a=n^2-m^2 = w_03,
b=2mn      = 2w_12,
c=n^2+m^2 = w_02+w_13 = w_23-w_01.                      (3)
```

Thus `(a,b,c)` is a Pythagorean triple.  If `gcd(m,n)=1` and the parities
differ, it is primitive.  Conversely every primitive Pythagorean triple has
such a window, up to the conventional order of its legs.

Equation (3) is not a bare matching statistic.  On the three perfect
matchings of `K4`, it uses four different operations:

```text
03|12 : one edge and twice the opposite edge;
02|13 : a sum;
01|23 : a signed difference.                             (4)
```

The matching quotient retains which edges are opposite but forgets the
coefficient, sum, and sign in (4).  Those operations are the current sidecar.

## 2. The exact golden/Pell slice

Define

```text
D(m,n)=n^2-mn-m^2.                                       (5)
```

For `n>m>0`,

```text
|D(m,n)|=1
iff (m,n)=(F_k,F_(k+1)) for a unique k>=2.                (6)
```

Including the fixed-cusp boundary `(1,1)` extends (6) to `k>=1`.  Moreover,

```text
D(F_k,F_(k+1))=(-1)^k.                                   (7)
```

Here is a complete descent proof.  If `n>=2m`, then writing `n=2m+r` gives

```text
D(m,n)=m^2+3mr+r^2.
```

Its absolute value is one only at `(m,n)=(1,2)`.  Otherwise a solution has
`0<n-m<m`, and

```text
D(n-m,m)=-D(m,n).                                        (8)
```

Repeatedly applying `(m,n)->(n-m,m)` strictly decreases the larger entry and
ends at `(1,1)` or `(1,2)`.  Reversing the descent applies

```text
G=[0 1;1 1],
```

which generates consecutive Fibonacci pairs.  This proves both directions
of (6), not merely a finite Pell search.

On this locus, (1) becomes the four-consecutive-number window

```text
W_k=(F_(k-1),F_k,F_(k+1),F_(k+2)).                       (9)
```

The branch in THM-3334 instead has `n=m+1`.  Since
`F_(k+1)-F_k=F_(k-1)`, the two loci meet only at `k=2,3`, namely parameters
`(1,2)` and `(2,3)`.  They then separate:

```text
U-spine next parameter:       (3,4) -> (7,24,25),
golden next raw parameter:    (3,5) -> (16,30,34).        (10)
```

Consequently THM-3334's fixed-discriminant `-4` quadratic
`m^2+(m+1)^2` and the present norm form of discriminant `5` are distinct
carriers with two initial intersections, not one recurrence in disguise.

## 3. Primitive normalization has period three

Use the ordered Euclid lift

```text
Psi(m,n)=(n^2-m^2,2mn,n^2+m^2).                          (11)
```

THM-3333's content law gives

```text
gcd(Psi(F_k,F_(k+1)))=2  iff k=1 mod 3,
                       =1  otherwise.                    (12)
```

Indeed Fibonacci residues modulo two have period three; both adjacent terms
are odd exactly in the displayed class.  At those indices define

```text
T=1/2[-1 1;1 1],
T(m,n)=((n-m)/2,(n+m)/2).                                (13)
```

The pair in (13) is integral, and

```text
Psi(T(m,n))=(mn,(n^2-m^2)/2,(n^2+m^2)/2).                (14)
```

Thus (14) is exactly the primitive normalization of (11), with the odd leg
restored to the first coordinate.  The first hostile is

```text
(m,n)=(3,5): Psi=(16,30,34), Psi(T(m,n))=(15,8,17).       (15)
```

Dividing without this leg/current sidecar would place the normalized triple
in the wrong Berggren coordinate chamber.

## 4. Exact transplant into three Berggren rays

Use THM-2596's parameter branch matrices

```text
A=[ 0 1;-1 2],       B=[0 1;1 2],       C=[1 0;2 1],      (16)
```

acting on columns `(m,n)^T`.  They induce the standard three Berggren child
matrices on (11).  Direct multiplication gives

```text
G^3=AB,                   T G^3 T^(-1)=CB.                (17)
```

Let `u_k=(F_k,F_(k+1))` and `u_2=(1,2)`.  Equations (12)--(17) imply

```text
u_(3r+2)       =(AB)^r u_2,
u_(3r+3)       =(AB)^r A u_2,
T u_(3r+4)     =(CB)^r C u_2.                             (18)
```

Matrix products act rightmost first.  Therefore the corresponding
root-to-child parameter words are exactly

```text
(BA)^r,                    A(BA)^r,                    C(BC)^r. (19)
```

The third ray is the odd/odd primitive-normalization transplant; it is not a
third power of the same raw positive matrix.  In particular Fibonacci's
period-three content law selects three interleaved ancestry rays.  It does
not define a `C3` action permuting Berggren's three children, consistently
with THM-2596 and THM-2632.

## 5. The intrinsic six-state object

Let

```text
a=(1,0),                    b=(0,1),                    c=(1,1)
```

be the three nonzero elements of `F_2^2`.  For each consecutive Farey flank
`{u_k,u_(k+1)}`, order its columns `(x_k,y_k)` so that
`det(x_k,y_k)=+1`, reduce modulo two, and apply THM-2632's channel-order map

```text
pi_k=(xbar_k,xbar_k+ybar_k,ybar_k).                       (20)
```

Starting at `k=2`, the result is the period-six cycle

```text
(b,c,a) -> (b,a,c) -> (a,b,c) -> (a,c,b)
        -> (c,a,b) -> (c,b,a) -> (b,c,a).                 (21)
```

Every step is one adjacent swap, alternating the two Coxeter generators of
`S3`.  These are all six total orders of the three `V4` channels.  Relative
to `(a,b,c)`,

```text
sgn(pi_k)=(-1)^k
         =F_(k-1)F_(k+2)-F_k F_(k+1).                    (22)
```

Thus Cassini's sign is exactly the sign of the channel order.  A chosen
bijection from `{a,b,c}` to the three perfect matchings turns it into
matching-order sign, but that bijection is an additional `S3` gauge.  Neither
(20) nor (22) selects one of the four affine `K4` vertices as owner.

## 6. What the tournaments do and do not remember

For `k>=3`, the six edge products of (9) have the total order

```text
w_01 < w_02 < {w_03,w_12} < w_13 < w_23,                 (23)
w_03 < w_12  for odd k,
w_12 < w_03  for even k.                                 (24)
```

All inequalities except (24) follow from positivity and monotonicity of the
Fibonacci sequence.  The central gap is exactly

```text
w_03-w_12=(-1)^k.                                        (25)
```

Hence the comparison orientation on the six edge labels is a transitive
`T6`, and one Fibonacci step flips exactly the arc between the opposite edges
`03` and `12`.

That flip is not induced by relabeling the four vertices.  THM-2753 proves
that every permutation induced by `S4` on the six edges lies in `A6`, whereas
the isolated transposition

```text
(03 12)                                                    (26)
```

is odd.  Therefore Cassini alternation is an edge-ranking/XOR sidecar, not
ambient `K4` transport.

There is also a literal four-vertex comparison tournament, but it is too
coarse.  For every `k>=3`, the entries in (9) give the same transitive order
`0<1<2<3`.  In particular

```text
(1,2,3,5),                    (2,3,5,8)                   (27)
```

have the same `T4`, but opposite Cassini sign, different primitive-content
class, and different rays in (19).  At the root window `(1,1,2,3)` there is a
tie, so even this `T4` does not exist without a tie-breaking gauge.

Nor is the matching quotient an ancestry symmetry.  The `V4` relabeling
`(01)(23)` fixes all three perfect matchings but sends

```text
(1,2,3,5) -> (2,1,5,3),                                  (28)
```

destroying the four-consecutive recurrence.  Finally, no tournament on the
four vertices is invariant under the full translation `V4`: translation by
a nonzero direction swaps the endpoints of every edge in that direction.

## 7. Exact affine-owner obstruction

The matching action is the standard exact sequence

```text
1 -> V4 -> S4 -> S3 -> 1.                                 (29)
```

Choose the two adjacent reflections that generate the quotient hexagon
(21).  Each has four lifts to `S4`, so there are sixteen ordered lift pairs.
The exact census is

```text
4 pairs:   two involutions generate S3 and fix one K4 vertex;
12 pairs:  at least one four-cycle is present and the pair generates S4.  (30)
```

The four complements in the first line are precisely the four point
stabilizers, one for each possible affine owner.  The twelve groups in the
second line are transitive and have no fixed owner.  All sixteen pairs project
to the same two quotient reflections and hence close the same six-step
matching loop.  Quotient loop closure cannot distinguish the two lines of
(30), much less choose one of the four owners in the first line.

Equations (26) and (30) give complementary obstructions:

- the intrinsic matching order has a four-way owner-lift ambiguity;
- the extra weighted Cassini flip is not an `S4` edge transport that could
  resolve that ambiguity.

An owner-fixed lift becomes lawful after explicitly choosing one vertex, and
a branch-dependent `V4` cocycle may be studied as extra data.  Neither is
canonically supplied by the Fibonacci values, the bare matchings, or the
six-state quotient.  This is the precise moving-owner boundary.

## 8. Relation to THM-3334 and stopping boundary

THM-3334 has an **external** `K4` when a hypotenuse such as
`1105=5*13*17` has four Gaussian parent representations.  Its vertices are
four distinct Berggren ancestors and its prime-XOR matchings compare those
parents.  The present **internal** `K4` consists of four entries of one window
and its matchings carry one triple by (3).  Equal cardinalities and the same
abstract matching geometry do not identify the two torsors.

The lawful source-target maps are:

| source | target/map | preserved | lost or required |
|---|---|---|---|
| golden Pell pair | primitive Euclid lift | triple and period-three content | raw odd/odd scale and leg order |
| primitive golden triple | word (19) | exact Berggren ancestry ray | no `C3` branch action |
| oriented Farey flank | mod-two order (20) | six channel states and Cassini sign | height, exact fraction, owner |
| ordered window | edge products plus (3) | full Pythagorean current | coefficients/signs after bare matching quotient |
| weighted six edges | comparison `T6` | Cassini middle arc | not an `S4` edge motion |
| matching loop | quotient (29) | `S3` order | `V4` cocycle and one of four owners |

Nothing here identifies the tail discriminant `-4` with a cubic
discriminant, transports a Keller map, supplies an LRC hull owner or phase,
or proves a global Berggren--Farey equivalence.  THM-3335 and THM-3336 are
listed only as **RESERVED** adjacent construction lanes; this theorem does
not depend on their intended statements.

## 9. Exact reproduction

Run

```bash
python 04-computation/fibonacci_berggren_three_ray_owner_thm3339.py
python -O 04-computation/fibonacci_berggren_three_ray_owner_thm3339.py
```

The dependency-free companion checks the symmetric-square branch action;
indices `2..201` on all three rays; the period-three content law and normalized
transplant; every Pell solution below `1000` as a hostile converse census;
the generalized `K4` current; the complete six-state residue cycle; edge
orders through `k=99`; all `24` `S4` edge/matching actions; the `V4` kernel;
the odd isolated-swap no-go; and all `16` owner-lift pairs.  Normal and
optimized runs byte-match the stored transcript and end in
`ALL CHECKS PASSED`.

QED.
