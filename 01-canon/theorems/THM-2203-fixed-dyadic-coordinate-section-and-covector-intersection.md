---
id: THM-2203
title: "Fixed dyadic coordinate section and bounded covector intersection"
status: >
  PROVED relative to LRCUpTo13 through THM-2199. In the actual live
  dyadic-terminal rank-eight lane, not in an arbitrary abstract rank-eight
  character cover, the scalar 5+3 coefficients occupy a fixed nine-coordinate
  section of the original speed row:
  (8H,16q_1,...,16q_5,208s_1,208s_2,208s_3).
  Terminal primitivity forces the scalar evaluation multiplier to be one.
  Thus the undivided and divided scalar lifts have fixed denominators 16 and
  208; no unknown row ratio occurs. More generally, twelve independent
  ambient relations of height H_0 yield eight independent relations supported
  on any chosen nine coordinates, of height at most 80H_0^5. For
  H_0=78*182^13 this supplies eight internal scalar relations, a bounded
  saturated integral basis, unrestricted fibres in every composite radix,
  and an explicit scalar-kernel minor bound. Combining the fixed coordinate
  section with the MSS height ceiling and the elementary largest-speed tooth
  bound forces every actual scalar blocker valuation to be at most 19.
  Hence THM-2192's unique-deepest scalar tail has only 1,140 valuation
  profiles, not an infinite all-depth ledger. THM-2199's direct global
  primitive box is strictly smaller than the internal-minor box, so the
  theorem closes the scalar transport/Smith/radix and infinitary-profile
  debts, not LRC(14), feasible enumeration, or the remaining exact scalar
  noncancellation problem.
source: codex-2026-07-24-fixed-dyadic-scalar-section
depends_on:
  - THM-759-tight-instance-ratio-bound
  - THM-763-strict-finite-height-for-tight-lrc-instances
  - THM-2073-lrc14-dyadic-deletion-tower
  - THM-2076-guard-capacity-terminal-rank-floor
  - THM-2077-terminal-kakeya-needle-and-recursive-quarter-escape
  - THM-2098-mixed-torus-collision-budget-and-vertical-gap
  - THM-2168-three-target-second-depth-majorization
  - THM-2187-height-preserving-saturation-and-universal-radix-carrier
  - THM-2199-effective-positive-subspace-rank-lift
related:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2196-bounded-relation-cone-circuit-atlas-and-carry-lock-rank-ladder
---

# THM-2203 -- fixed dyadic scalar section

The apparent arithmetic lift in the scalar `5+3` tail is not arbitrary. The
dyadic deletion tower already locates all nine scalar characters in fixed
coordinates of the original thirteen-speed row. Pulling **covectors** through
that section, rather than manufacturing a covariant lift from row ratios,
removes the denominator and Plucker-height ambiguity in THM-2192 Section 6.

The provenance qualifier is load-bearing throughout: this theorem concerns
the actual live `n=8` terminal produced by the THM-2073 tower. An abstract
guard/eight-terminal character cover need not come with this dyadic section.

Put

```text
H_0=78*182^13=18750922831149193194381342621696.       (1)
```

## 1. The scalar coefficients are original coordinates

THM-2073 starts from

```text
S=2C disjoint_union {x,y},
C=2^r Q_r disjoint_union {2^i h_i:0<=i<r},            (2)
```

where every `h_i` is odd and every `Q_i`, in particular `Q_r`, is primitive.
THM-2076 gives

```text
|Q_r|=11-r                                             (3)
```

and identifies the last guard as `h_(r-1)`. THM-2098 Section 3 identifies its
guarded `n=8` application with the actual live dyadic-terminal lane. Hence

```text
r=3,                  guard=h_2,                       (4)
```

and unwinding (2) gives

```text
C=8Q_3 disjoint_union {h_0,2h_1,4h_2},
S=16Q_3 disjoint_union {2h_0,4h_1,8h_2,x,y}.           (5)
```

On THM-2168's surviving scalar line, let `alpha` be the positive primitive
generator and write

```text
g=H alpha,                  c_i=q_i alpha,       1<=i<=5,
u_j=s_j alpha,              c_*j=13s_j alpha,   1<=j<=3. (6)
```

Let `d` be the actual integral evaluation cocharacter and put

```text
N=alpha(d)>0.                                           (7)
```

The guard and eight terminal evaluations in (4)--(6) give

```text
h_2=NH,
Q_3=N(q_1,...,q_5,13s_1,13s_2,13s_3).                 (8)
```

All entries in the coefficient tuple are integers. Since `Q_3` is primitive,

```text
1=gcd(Q_3)=N gcd(q_1,...,q_5,13s_1,13s_2,13s_3).      (9)
```

Thus

```text
N=1.                                                   (10)
```

Let `I` be the nine original coordinate labels consisting of the last guard
and the eight terminals. Equations (5), (8), and (10) give the fixed section

```text
S_I=(8H,16q_1,...,16q_5,208s_1,208s_2,208s_3)
   =8(H,2q_1,...,2q_5,26s_1,26s_2,26s_3).             (11)
```

In particular, for

```text
w_*=(H,q_1,...,q_5,13s_1,13s_2,13s_3),                (12)
```

the coordinate injection

```text
J_* e_H=e_(h_2)/8,             J_* e_(q_i)=e_i/16,
J_* e_(13s_j)=e_(*j)/16                                (13)
```

satisfies `S J_*=w_*` and has denominator `16`. For THM-2192's divided row

```text
w=(H,q_1,...,q_5,s_1,s_2,s_3),                        (14)
```

replace the last line of (13) by

```text
J e_(s_j)=e_(*j)/208.                                  (15)
```

Then `S J=w`, the denominator is `208`, and the inverse on the coordinate
image is the fixed diagonal `(8,16,...,16,208,208,208)`. Neither map contains
an unknown ratio of entries of `S` or `w`.

## 2. Coordinate-supported intersection lemma

The following elementary lemma explains why a bounded ambient relation basis
can now be transferred internally.

> **Lemma.** Let `v in Q^13` have no zero coordinate, let
>
> ```text
> B in Z^(12 x 13),       rank(B)=12,       Bv^T=0,
> ||B||_infinity<=L,                                      (16)
> ```
>
> and choose any nine-coordinate set `I`. Then the row space of `B` contains
> eight independent integer relations supported on `I`, each of height at
> most
>
> ```text
> 80L^5.                                                  (17)
> ```

Let `O=I^c`, so `|O|=4`. Since `B` has rank twelve and annihilates the
nonzero row `v`,

```text
row_Q(B)=v^perp.                                        (18)
```

The coordinate projection `v^perp -> Q^O` is onto. Indeed, prescribe the four
coordinates on `O` and correct any one coordinate in `I` by division by its
nonzero entry of `v`. Therefore the restriction `B_O` has row rank four.

Put

```text
C=B_O^T in Z^(4 x 12).                                 (19)
```

Choose four pivot columns with nonzero determinant. For each of the eight
nonpivot columns, the signed `4 x 4` cofactors of that column together with
the four pivots give a five-supported vector

```text
lambda_j in ker(C).                                    (20)
```

Its nonpivot entry is the fixed nonzero pivot determinant. Hence the eight
vectors `lambda_j` are independent. Every cofactor is a determinant of a
`4 x 4` matrix with entries bounded by `L`; Hadamard gives

```text
||lambda_j||_infinity<=(2L)^4=16L^4.                   (21)
```

Now put

```text
x_j=lambda_j^T B.                                      (22)
```

The map `lambda -> lambda^T B` is injective because the rows of `B` are
independent, so the eight `x_j` remain independent. Equations (19)--(20)
make them vanish on `O`, and their five-supported construction gives

```text
||x_j||_infinity<=5*16L^4*L=80L^5.                    (23)
```

This proves the lemma.

## 3. Eight bounded internal scalar relations

Assume that the primitive positive row `S` has zero `1/14`-safe Haar mass.
THM-2199 supplies twelve independent integer relations of height at most
`H_0`. Place them in the rows of `B` and apply the lemma to the fixed set `I`
from Section 1. We obtain eight independent relations supported on the nine
scalar coordinates, all of height at most

```text
K=80H_0^5.                                             (24)
```

Define the primitive integer row

```text
u=(H,2q_1,...,2q_5,26s_1,26s_2,26s_3).                (25)
```

It is primitive: any odd common divisor would divide every entry of the
primitive tuple `Q_3=(q_i,13s_j)`, while the prime two cannot divide the odd
guard `H=h_2`. Since `S_I=8u`, all eight relations from (24) annihilate `u`.
They form an `8 x 9` rank-eight matrix whose rational kernel is `Q u`.
Signed maximal minors and Hadamard therefore give

```text
max(H,2q_i,26s_j)
 <=(sqrt(8)K)^8
 =8^4(80H_0^5)^8.                                     (26)
```

The same relations, expressed on the two scalar rows (12) and (14), have
respective coefficient heights

```text
2K=160H_0^5,
26K=2080H_0^5.                                        (27)
```

Thus THM-2192's requested eight-dimensional scalar kernel is reached by a
fixed coordinate operation. The hostile rational planes

```text
span_Q(e_1,...,e_7,Me_8+e_9)
```

still show that an arbitrary eight-plane can have unbounded lattice height
inside a bounded ambient lattice. They do not obstruct the actual scalar
lane, whose diagonal section (13)--(15) is fixed before the row varies.

### Saturated basis and universal radix

The relation lattice

```text
M_u=u^perp intersection Z^9                             (28)
```

is saturated of rank eight because `u` is primitive. Starting with a
primitive divisor of the first relation from (24), adjoin the other seven
independent relations one at a time and apply THM-2187 equation (25). The
resulting basis-height profile, in units of `K`, is at most

```text
(1,1,3/2,9/4,27/8,81/16,243/32,729/64).              (29)
```

Thus `M_u` has an integral basis of height at most

```text
(729/64)K=(3645/4)H_0^5.                              (30)
```

For `w_*` and the divided row `w`, the same construction gives saturated
basis heights at most

```text
(3645/2)H_0^5,              (47385/2)H_0^5,           (31)
```

respectively. THM-2187 equations (7)--(9) then give an integral right inverse:
modulo every integer `m>=2`, each eight-row digit map from nine coordinates
onto eight relation digits is surjective and every unrestricted fibre has
exactly `m` elements. This closes the Smith and composite-radix interfaces,
but owner-restricted fibres still require the phase/current sidecar.

## 4. The all-depth scalar ledger is actually finite

The fixed factors in (11) expose a consequence which is invisible in the
abstract scalar cover.  Let `S` be a primitive hypothetical LRC(14)
counterexample, so

```text
M(S)<=1/14.
```

THM-763, applied with thirteen nonzero speeds and the settled
lower-dimensional LRC, gives

```text
sum_(v in S) v <= binom(14,2)^12=91^12
 =322475487413604782665681.                           (32)
```

There is also a useful non-strict version of THM-759's danger-tooth
argument.  Write `w` for the largest speed, `w_2` for the second largest,
and `S_0=S\{w}`.  The settled twelve-speed theorem gives

```text
mu=M(S_0)>=1/13.
```

If

```text
w>w_2/(14mu-1),
```

choose a maximizing phase of `S_0` and a circular interval radius `rho`
strictly between

```text
1/(14w)  and  (mu-1/14)/w_2.
```

Every speed in `S_0` stays strictly `1/14`-safe on that interval, while the
interval is wider than one closed danger tooth of `w`.  It therefore contains
a phase at which `w` is strictly safe as well, contradicting
`M(S)<=1/14`.  Consequently

```text
w<=w_2/(14mu-1)<=13w_2,
w<=13/14 sum_(v in S)v <=13/14*91^12.                (33)
```

Return to the actual scalar section.  Its three original blocker coordinates
are `208s_j`.  If

```text
lambda_j=nu_13(13s_j),
```

then positivity gives

```text
208*13^(lambda_j-1)<=208s_j<=w.                      (34)
```

The first excluded valuation is separated by an exact integer gap:

```text
208*13^19
 =304079420398092791020816
 >13/14*91^12
 =598883048053837453521979/2.                        (35)
```

Equations (33)--(35) prove

```text
lambda_j<=19                  for j=1,2,3.            (36)
```

THM-2192 orders the actual blocker depths as

```text
1<=lambda_1<=lambda_2<lambda_3.
```

For a fixed deepest value `c`, the shallower pair is a two-element multiset
from `{1,...,c-1}`, so there are `binom(c,2)` choices.  Therefore the entire
valuation ledger has

```text
sum_(c=2)^19 binom(c,2)=binom(20,3)=1140             (37)
```

profiles.  Exact exclusions such as THM-2198 and THM-2204/2205 remove rows
from this finite ledger.  The number `1140` is not a feasible full
root-capacity enumeration: the largest torsion layer still has order
`13^20`, and the unit labels and owner-current sidecars remain enormous.
The gain is a quantifier correction.  “All deeper profiles” means a finite
depth-`19` list in the actual dyadic LRC lane, not an unbounded family.

## 5. Numerical comparison and exact remaining debt

THM-2077 equation (13), at `r=3`, propagates (26) back through the tower:

```text
max(S)<=(144/5)max(Q_3)
       <=(72/5) 8^4(80H_0^5)^8.                       (38)
```

This is an explicit internal-carrier box, but it is not a new or improved
finite reduction. THM-2199 already gives the much smaller direct bound

```text
max(S)<=12^6 H_0^12                                   (39)
```

for every primitive zero-Haar row, and (11) then bounds `u` immediately.
The exponent `12` in (39) is strictly better than the exponent `40` in
(26)/(38). The value of Sections 1--4 is instead structural:

```text
source:          actual dyadic-terminal n=8 row;
map:             fixed diagonal coordinate section;
preserved:       evaluation, scalar kernel, integrality after clearing 208;
destroyed:       four outer dyadic/tail coordinates;
needed sidecar:  none for scalar height transport;
decisive test:   eight supported covectors from the four-coordinate
                 elimination lemma.                                  (40)
```

Consequently the scalar-lift and denominator/Plucker-height debts in
THM-2192 are closed. The remaining LRC task is not finiteness. It is to empty
the finite primitive locus, or equivalently to prove exact noncancellation
for the `216` surviving residue profiles and deeper unique-maximum valuation
patterns. This theorem does not perform that discharge and is not a proof of
LRC(14). QED.
