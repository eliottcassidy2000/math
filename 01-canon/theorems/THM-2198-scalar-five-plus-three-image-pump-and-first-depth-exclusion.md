---
id: THM-2198
title: "Scalar five-plus-three image pump and first-depth exclusion"
status: >
  PROVED + VERIFIED-EXACT. Multiplication by thirteen has an exact labelled
  root-image law on the scalar five-unit/three-deep cover. On every
  thirteen-root fibre, a unit mask is a moving singleton/chord given by a
  two-integer window, the guard is the complement of a moving
  three/four-integer block, and a 13-divisible mask is one all-or-nothing
  fibre bit after its coefficient is divided by thirteen. The resulting
  insertion/deletion word records labelled ownership, simultaneous event
  packets, full integer winding, and the active bits of the two shallower
  blockers left after the unique deepest blocker becomes safe. Independently,
  the mixed residual expands in measure by at least 13/10 on its first image
  and never contracts thereafter. Division support yields quantitative and
  private-owner transition words through the first two blocker valuations,
  with equality carriers (1,13) and (1,12,13). On the first possible
  unique-depth profile (1,1,2), an exact depth-three audit gives a
  global sheet-time incidence deficit: for every pair of shallower masks,
  the five largest conditional unit-mask capacities miss the residual by at
  least 86 points. Hence (1,1,2) is empty and every surviving scalar 5+3
  branch has unique deepest blocker depth at least three. The signed
  residual is exactly a -1 eigenfunction of the thirteen-root transfer;
  over a divided-deep safe parent, a positive defect has one labelled
  deep-danger anti-defect root and twelve neutral roots. This theorem alone
  leaves deeper profiles; THM-2204/2205/2207 later close all profiles of
  deepest valuation three. This is not a proof of LRC(14).
source: codex-2026-07-24-scalar-five-plus-three-transition
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2137-deep-scalar-tail-boundary-complexity
  - THM-2138-all-depth-unit-annulus-extremal-law
  - THM-2168-three-target-second-depth-majorization
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
related:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2184-two-scale-tail-continuation-profile
  - THM-2196-bounded-relation-cone-circuit-atlas-and-carry-lock-rank-ladder
  - THM-2199-effective-positive-subspace-rank-lift
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-2200-convex-semigroup-and-finite-place-support-hole-trichotomy
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
script: 04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py
output: 05-knowledge/results/lrc14_scalar_five_plus_three_image_pump_thm2198.out
script_sha256: b3f5b37d4e0fb34070851fcedecf9435e08d9f78fafb651c5a6c4b39528e0eee
output_sha256: 5e2130ac4078da4c3e4f600f7c29e1f8290cea2ffcbdc012e144117d9d786d76
hash_basis: working-tree bytes (LF)
---

# THM-2198 -- the scalar `5+3` image pump

Use the scalar-cover notation

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                          (1)
```

Suppose that `H,q_1,...,q_5` are positive thirteen-units, `H` is odd, the
`q_i` are pairwise distinct, the three actual blocker coefficients
`c_1,c_2,c_3` are pairwise distinct positive multiples of thirteen, and

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere.  This is THM-2168's fully scalar `5+3` residual.  Put

```text
lambda_j=nu_13(c_j).
```

THM-2192 proves that a survivor can be labelled so that

```text
1<=lambda_1<=lambda_2<lambda_3.                      (3)
```

We first prove the exact all-depth transition law missing from the static
carrier, then exclude the first possible profile `(1,1,2)`.

## 1. Exact root-image disintegration

Let

```text
rho:R/Z -> R/Z,                   rho(x)=13x.
```

Choose a representative `0<=y<1`.  Its thirteen roots are

```text
x_k=(y+k)/13,                    k in F_13.           (4)
```

Label the root sheets in the **reversed guard trivialization** by

```text
s=-Hk mod 13.                                         (5)
```

This is the sign-reverse of THM-2192's convention in which increasing the
sheet label increases the guard value.  The reversal makes the integer
window itself run in the positive direction; no unoriented chord length or
coverage predicate changes.

For any thirteen-unit `q`, put

```text
I_q(y)=Z intersection (qy-13/14,qy+13/14),
delta_q=Hq^(-1) mod 13,                               (6)
```

where the inverse in `delta_q` is taken modulo thirteen.  Put also

```text
J_H(y)=Z intersection (Hy-13/7,Hy+13/7).              (7)
```

Then the exact fibre incidences are

```text
A_q(y)={s:x_k in D_q}
      ={delta_q d mod 13:d in I_q(y)},                (8)

F_H(y)={s:x_k notin C_H}
      ={d mod 13:d in J_H(y)},                        (9)

B_H(y)={s:x_k in C_H}=F_13\F_H(y).                   (10)
```

Indeed, `x_k in D_q` precisely when some `n in Z` satisfies

```text
|q(y+k)-13n|<13/14.
```

Writing `d=13n-qk` gives

```text
|qy-d|<13/14,             d=-qk mod 13.
```

Equation (5) now gives `s=Hq^(-1)d`, proving (8).  The identical argument
with radius `13/7` and `q=H` proves (9).

A blocker divisible by thirteen behaves differently.  If `c=13v`, then

```text
cx_k=v(y+k)=vy mod 1.                                 (11)
```

Consequently

```text
{s:x_k in D_(13v)}
 =F_13 if y in D_v,
 =emptyset otherwise.                                (12)
```

Thus multiplication by thirteen converts every divisible blocker into one
all-or-nothing fibre bit and divides its coefficient by thirteen.  This is
the **image pump**.

The same formula is exact on every primitive torsion layer.  For `m>=2`,
put `N=13^m`, `Q=13^(m-1)`, and reduce a primitive numerator

```text
z mod N -> r=z mod Q.
```

The thirteen lifts of `r/Q` are exactly

```text
z=r+kQ,                       k in F_13,              (13)
```

and remain primitive.  Equations (8)--(12), with `y=r/Q`, disintegrate the
primitive guard-safe annulus at depth `m` into its depth-`m-1` image phases.
No measure approximation or phase mesh occurs.

## 2. The labelled event word, ownership current, and winding

The integer windows make the continuous transition law explicit.  Away from
boundary events, if `n in Z`, then

```text
qy in (n-1/14,n+1/14)       => I_q(y)={n},
qy in (n+1/14,n+13/14)      => I_q(y)={n,n+1}.        (14)
```

As `qy` increases, the crossing

```text
qy=n+1/14
```

inserts the endpoint labelled `delta_q(n+1)`, and the crossing

```text
qy=n+13/14
```

deletes the endpoint labelled `delta_q n`.  Thus a unit mask evolves by

```text
singleton -> chord -> next singleton,                (15)
```

with toothpick step `delta_q`.  Over one turn of `y`, it has `q`
insertions and `q` deletions.  Its net sheet displacement in the
trivialization (4)--(5) is

```text
q delta_q=H mod 13.                                  (16)
```

The residue step alone therefore does not determine the movie: the full
positive integer `q` is its winding sidecar.

The guard has the parallel law

```text
Hy in (n-1/7,n+1/7)
       =>J_H(y)={n-1,n,n+1},
Hy in (n+1/7,n+6/7)
       =>J_H(y)={n-1,n,n+1,n+2}.                     (17)
```

At `Hy=n+1/7`, guard-unsafe sheet `n+2` is inserted; at
`Hy=n+6/7`, sheet `n-1` is deleted.  Hence `B_H(y)` has nine or ten
sheets, exactly as in THM-2192.

Finally, the bit belonging to `13v` turns on at

```text
vy=n-1/14
```

and off at `vy=n+1/14`.  If several displayed rational event times agree,
they form one simultaneous labelled event packet.  They are not
arbitrarily tie-broken into a tournament.

For the five unit masks define the ownership current and deficiency

```text
omega_s(y)=#{i:s in A_(q_i)(y)},
Z(y)=B_H(y)\{s:omega_s(y)>0}.                         (18)
```

Equations (14)--(17) update `omega`, `B_H`, and `Z` exactly, one labelled
event packet at a time.  Away from the image of the null exceptional set in
(2), equation (12) says that fibre coverage is automatic when at least one
deep bit is on; when every deep bit is off, fibre coverage is exactly

```text
Z(y)=emptyset.                                        (19)
```

This supplies the data absent from THM-2197's static Boolean quotient:
deleting an endpoint updates its labelled ownership count, the cyclic event
word says which update occurs next, and the integer coefficients retain
winding.  Chord lengths or the deficiency set alone do not determine this
transition.

There is an exact functional form for the corresponding `H`-drift.  Let
`P` be any finite generic set of image phases on which every deep bit is
off.  Define

```text
e_i(y)=1_(|I_(q_i)(y)|=1),
o_i(y)=|A_(q_i)(y) intersection F_H(y)|,
chi_C(y)=1_(y in C_H),                               (20)

H_P=sum_(y in P)|B_H(y)|
       -sum_(i=1)^5 sum_(y in P)|A_(q_i)(y)
                                      intersection B_H(y)|.
```

Since a unit mask has `2-e_i(y)` endpoints and

```text
|B_H(y)|=10-chi_C(y),
```

one has the identity

```text
H_P=sum_(i,y)e_i(y)+sum_(i,y)o_i(y)
                         -sum_y chi_C(y).             (21)
```

Thus the functional is exactly

```text
singleton-event loss + guard-outside endpoint loss
                         - nine-sheet C-frame credit.
```

Fibrewise coverage implies `H_P<=0`; the converse need not hold because
incidence overlap is forgotten.  This is the signed event invoice behind
the toothpick self-similarity, not an independence heuristic.

## 3. The unique deepest blocker leaves two shallower processes

Take the primitive layer

```text
N=13^(lambda_3+1).
```

Write `c_3=13^lambda_3 u_3` with `u_3` a thirteen-unit.  At a primitive
point `t=z/N`,

```text
c_3t=u_3z/13.
```

This is a nonzero thirteenth root, whose circle norm is at least
`1/13>1/14`.  The unique deepest blocker is therefore safe throughout this
layer.

Under the first image map, the other two blockers become the all-or-nothing
bits

```text
1_(D_(c_1/13))(y),          1_(D_(c_2/13))(y).       (22)
```

On an image phase where either bit is on, that deep bit covers the whole
thirteen-root fibre.  Outside the null image of the exceptional set, on the
residual phase set where both bits are off, the five labelled unit-mask
movies must cover every sheet of `B_H(y)`.  This is the precise
two-shallower-mask obstruction: it is a global matching movie over the
common residual image phases, not two extra chords inside one static fibre.

## 4. The first unique depth `(1,1,2)` is empty

Suppose for contradiction that

```text
(lambda_1,lambda_2,lambda_3)=(1,1,2).                (23)
```

Set

```text
N=13^3=2197,                     Q=13^2=169.
```

Multiplying primitive numerators by `H` modulo `N` normalizes the guard to
one and replaces every coefficient by its product with `H^(-1)`.  This is a
bijection and preserves all three valuations.  The primitive guard-safe
universe is

```text
U={z mod N:13 does not divide z and 7||z||_N>N},
|U|=1450.                                            (24)
```

Here `||z||_N=min(z mod N,-z mod N)`.  The unit sign classes are

```text
(Z/NZ)^*/{+/-1},                 1014 classes,       (25)
```

and the unit parts of a depth-one blocker run through

```text
(Z/QZ)^*/{+/-1},                 78 classes.         (26)
```

Both possible shallow blockers are retained even when they reduce to the
same sign class: the exact audit uses all

```text
C(78+1,2)=3081
```

unordered pairs including the diagonal.

For a unit label `a mod N` and a shallow unit part `u mod Q`, define

```text
S_a={z in U:14||az||_N<N},
T_u={z in U:14||13uz||_N<N}.                         (27)
```

For every shallow pair put

```text
R_(u,v)=U\(T_u union T_v).                            (28)
```

Let

```text
c_1(u,v)>=...>=c_5(u,v)
```

be the five largest cardinalities obtained by first taking the globally
distinct unit masks `S_a` and then intersecting each with `R_(u,v)`:

```text
S_a intersection R_(u,v).
```

Two globally distinct masks may coincide after this intersection; retaining
both only enlarges the deliberately generous sum bound.  The checker also
verifies that quotienting the sole global duplicate does not change the five
largest entries.

The exact exhaustive result is

```text
|R_(u,v)|-sum_(i=1)^5 c_i(u,v)>=86                  (29)
```

for all `3081` pairs.  The minimum is unique up to the fixed sign-class
representatives:

```text
(u,v)=(14,46),                 |R_(u,v)|=1046,

(c_1,...,c_5)=(204,200,190,186,180),                 (30)

corresponding unit labels=(183,799,599,1000,1007).
```

The image-pump form of the worst residual has `112` active base phases:

```text
74 fibres of size 9 and 38 fibres of size 10,
74*9+38*10=1046.                                    (31)
```

Thus (29) is a global sheet-time incidence, or zeroth Hall, deficit across
the entire matching movie.  Even before overlaps among the five chosen
movies are subtracted, their total conditional capacity is at most `960`,
leaving `86` residual sheet-time vertices.

For the five extremal labels in (30), the rows

```text
(label,capacity,singleton losses,outside-endpoint losses)

((183, 204,  0,20),
 (799, 200,  0,24),
 (599, 190, 10,24),
 (1000,186,  0,38),
 (1007,180, 12,32)).                                 (32)
```

The `74` nine-sheet frames are exactly the `C_H` frames.  Hence (21) reads

```text
H_P=(20+24+34+38+44)-74=86.                          (33)
```

This is the promised exact `H`-drift at the hostile boundary.

To see why (29) proves the theorem, reduce the five actual unit coefficients
modulo `N` and sign.  Repeated residue masks cannot enlarge a union, so at
most five distinct sets `S_a` occur.  Their union inside any `R_(u,v)` has
cardinality at most the sum of the five largest individual cardinalities,
which is strictly below `|R_(u,v)|` by (29).  Hence some `z in U` is safe
for all eight terminal masks.

There are no torsion endpoints: `N` and `Q` are powers of thirteen and are
coprime to seven and fourteen.  The uncovered point therefore satisfies
every guard and terminal inequality strictly and thickens to an uncovered
open interval.  This contradicts the almost-everywhere cover (2), proving
that (23) is impossible.

By (3), depth two was the first possible value of the unique deepest
blocker.  We have proved

```text
lambda_3>=3                                            (34)
```

in every surviving scalar `5+3` branch.

## 5. Quantitative measure pump and private owner words

The labelled root law has a complementary measure consequence.  Define the
mixed residual

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),                 (M1)
A_r=T^r(A_0),                     T(t)=13t.
```

For a Borel set `A subset C_H`, let

```text
n_A(y)=#{x in A:T(x)=y}.                             (M2)
```

Equations (9)--(10) give `0<=n_A(y)<=10`, and, up to the finite endpoint
set,

```text
{y:n_A(y)>0}=T(A).
```

The set `T(A)` is Borel: on the thirteen standard sheets, `T` is a Borel
isomorphism, so the image is a finite union of Borel images.  Haar
disintegration gives

```text
measure(A)=(1/13) integral n_A(y)dy
          <=(10/13)measure(T(A)),

measure(T(A))>=(13/10)measure(A).                    (M3)
```

For every Borel `B`, one also has

```text
measure(T(B))>=measure(B).                           (M4)
```

Indeed `B subset T^(-1)(T(B))`, while Haar invariance gives
`measure(T^(-1)(E))=measure(E)`.  Multiplication by thirteen is Lipschitz
and sends null sets to null sets.

THM-2137 gives the sharp five-unit residual floor

```text
measure(A_0)>=961/6930.                              (M5)
```

Consequently

```text
measure(A_r)>=12493/69300             for every r>=1. (M6)
```

The first factor is the exact `13/10` root-image expansion; all later images
use noncontraction (M4).  In the ownership coordinates, `n_(A_0)(y)` is the
cardinality of the safe-sheet deficiency.  Thus the support image and its
occupancy are both retained.

There is also exact division support.  If `b=13^r c`, then

```text
x in D_b iff T^r(x) in D_c.                          (M7)
```

There is an exact finite-place formulation of the first image. Regard the
integer-valued transfer

```text
T_*1_A=n_A
```

as a function with values in `F_13`. Since `A subset C_H` forces
`0<=n_A<=10<13`,

```text
supp_(F_13)(T_*1_A)=T(A).                            (M7a)
```

Thus reduction modulo thirteen loses no occupied first-image fibre. Together
with

```text
integral n_A=13 measure(A),       ||n_A||_infinity<=10, (M7b)
```

this is exactly the support/`L1` content of (M3). It is the set-transfer
analogue of THM-2022's whole-layer finite-place survival, while THM-2200's
Frobenius law concerns independent polynomial copies.

The analogy stops after one image. A general fibre has occupancy at most
thirteen, and occupancy thirteen becomes zero in `F_13`; moreover the five
unit coefficients are not integrally divisible under `T`. This is why later
steps use only (M4) and the labelled event word rather than iterating (M7a).

Applying (M7) to (2), and discarding only the image of the null exceptional
set, gives

```text
A_r subset union_(j=1)^3 D_(c_j/13^r)
                         for 0<=r<=lambda_1.          (M8)
```

### 5.1. A unique least valuation

Suppose first that the whole valuation ladder is strict:

```text
lambda_1<lambda_2<lambda_3.
```

Put

```text
a_1=c_1/13^lambda_1,
R_1=A_(lambda_1) minus D_(a_1).                      (M9)
```

Since `a_1` is a thirteen-unit and every danger comb has measure `1/7`,
(M6), (M8), and (M9) give

```text
measure(R_1)>=12493/69300-1/7
             =2593/69300.                            (M10)
```

For every `0<=u<=lambda_2-lambda_1`,

```text
T^u(R_1)
 subset D_(c_2/13^(lambda_1+u))
       union D_(c_3/13^(lambda_1+u)),

measure(T^u(R_1))>=2593/69300.                       (M11)
```

At the second valuation, define

```text
a_2=c_2/13^lambda_2,
d_3=c_3/13^lambda_2,
R_2=T^(lambda_2-lambda_1)(R_1).                      (M12)
```

Then `a_2` is a thirteen-unit, `13|d_3`, and the sharp pair floor in
THM-1166 yields

```text
2593/69300<=measure(R_2)
             <=measure(D_(a_2) union D_(d_3))
             <=25/91.                                (M13)
```

Equality in the last inequality occurs only for reduced ratio `(1,13)`.
Thus saturation forces

```text
d_3=13a_2,
lambda_3=lambda_2+1
```

and equality of the remaining unit parts.  This classifies the equality
carrier; it does not claim that a survivor must saturate (M13).

### 5.2. A deepest-private transition

THM-2138 says that five unit masks and any two positive-valuation masks
cannot cover `C_H`.  Hence

```text
P_3=A_0 minus (D_(c_1) union D_(c_2))                (M14)
```

contains a nonempty open set.  The cover (2) makes `c_3` its unique deep
owner:

```text
P_3 subset D_(c_3)                                   (M15)
```

almost everywhere.  Its images retain positive measure by (M4), and (M7)
gives

```text
T^r(P_3) subset D_(c_3/13^r)
                         for 0<=r<=lambda_3,          (M16)

T^r(P_3) intersection D_(c_j/13^r)=empty
                         for j=1,2 and r<=lambda_j,  (M17)
```

up to null boundaries.

If the least valuation is repeated,

```text
lambda_1=lambda_2<lambda_3,
```

then at `r=lambda_1` the blocker-owner word is

```text
(0,0,1).                                             (M18)
```

At that level, (M6), (M8), and THM-1166 give

```text
12493/69300<=measure(A_r)<=36/91.                    (M19)
```

Equality in the upper bound forces the normalized triple to be a scaled
permutation of `(1,12,13)`.  Since exactly two entries are thirteen-units,
this has the oriented form

```text
{c_1/13^r,c_2/13^r}={g,12g},
c_3/13^r=13g,                                       (M20)
```

for a thirteen-unit `g`; in particular `lambda_3-r=1`.

Under the strict ladder, (M17) first gives

```text
(0,0,1) at r=lambda_1.                               (M21)
```

At `r=lambda_2`, the first blocker is no longer integrally divisible and
its old zero bit is genuinely lost.  The surviving word is

```text
(*,0,1).                                             (M22)
```

Indeed `T^lambda_1(P_3) subset R_1`, so
`T^lambda_2(P_3)` lies in `R_2 minus D_(a_2)` and has positive measure.
The quantitative statement (M10) controls the two-carrier remainder; the
private deepest-owned subpiece is topologically positive but has no uniform
mass floor proved here.  The star records a destroyed integral-division
coordinate, not a cosmetic wildcard.

## 6. Signed residual eigenlaw and the one-anti-defect root

The labelled movie has a complementary exact signed law. Put

```text
v_j=c_j/13,
B(y)=sum_(j=1)^3 1_(D_(v_j))(y),
R(x)=1_(C_H)(x)-sum_(i=1)^5 1_(D_(q_i))(x),
(cal L f)(y)=sum_(T(x)=y) f(x).                       (S1)
```

For a thirteen-unit `q`, the root windows (6)--(8) give, away from
endpoints,

```text
cal L 1_(D_q)=2-1_(D_q).                              (S2)
```

The guard-root count (9)--(10) is equivalently

```text
cal L 1_(C_H)=10-1_(C_H).                             (S3)
```

There are exactly five unit masks, so the constant terms cancel:

```text
cal L R=-R.                                           (S4)
```

The cover (2) also gives the almost-everywhere pointwise inequality

```text
R<=B o T.                                             (S5)
```

Indeed, `R` is positive exactly on `A_0`, where some covering blocker
`c_j=13v_j` makes `v_jT(x)=c_jx` dangerous; off `A_0`, `R<=0`.
Applying `cal L` to (S5), and using `cal L(B o T)=13B`, gives

```text
-13B<=R<=B o T.                                       (S6)
```

Thus a positive residual over the divided-deep safe state `{B=0}` forces
its image into `{B>0}`. In particular,

```text
B(y)=B(Ty)=0        ==>        R(y)=0.                (S7)
```

The deepest blocker itself can pay this transition. Write

```text
c_3=13^beta u,        13 does not divide u,
v_3=c_3/13.
```

For `x=z/13^(beta+1)`, `13` not dividing `z`, and `y=T(x)`,

```text
v_3y=uz/13                  (safe),
v_3T(y)=uz                  (maximally dangerous).    (S8)
```

Multiplication by thirteen has erased the nonzero root digit. A recursive
state cannot infer that one of the two shallower blockers pays (S6) unless
it retains this pre-collapse `F_13` digit.

More precisely, let `y` be generic with `B(y)=0`. Equation (S6) and
`R<=1` give `R(y) in {0,1}`. Every root `x` of `y` satisfies
`R(x)<=B(Tx)=0`, while (S4) says that the thirteen integer residuals sum
to `-R(y)`. Hence

```text
R(y)=0  ==>  {R(x):Tx=y}={0^13},
R(y)=1  ==>  {R(x):Tx=y}={-1,0^12}.                  (S9)
```

In the second line the unique `-1` root has `B(x)>0`; otherwise the left
side of (S6) at `x` would force `R(x)>=0`. A positive defect over a
deep-safe parent therefore emits one canonically labelled deep-danger
anti-defect root and twelve neutral roots. Equation (S9) is the exact
one-step desuspended automaton.

### 6.1. The first lifted digit gate

The remaining repeated-minimum profiles have the normalized form

```text
c_1=13a,       c_2=13b,       c_3=13^k c,       k>=2, (S10)
```

with `a,b,c` thirteen-units. Put `v_3=13^(k-1)c` and take a critical
primitive parent

```text
y=z/13^k,             13 does not divide z.
```

Then `v_3y=cz/13` is safe. On its thirteen roots

```text
x_r=(z+r13^k)/13^(k+1),             r in F_13,
```

the deepest divided phase is independent of the root digit:

```text
v_3x_r=cz/169 mod 1.                                (S11)
```

Suppose also that `a,b` are safe at `y` and `R(y)=1`. Then `B(y)=0`, so
(S9) supplies a unique anti-defect root `x_kappa` with `B(x_kappa)>0`.
There are exactly `156` primitive classes modulo `169`. On the `132`
classes satisfying

```text
14||cz||_169>169,
```

the deepest divided comb is safe on all thirteen roots, and therefore
`x_kappa` must lie in `D_a union D_b`: one of the two shallower blockers
is forced to shelter the exceptional digit. On the remaining `24` classes,
the deepest comb is dangerous on every root and can shelter the digit by
itself. The count is exact: these are the twelve positive and twelve
negative primitive residues of norm below `169/14`.

More generally, on every `r`-fold descendant fibre with `1<=r<=k-1`,

```text
x=(z+n13^k)/13^(k+r)
       ==> v_3x=cz/13^(r+1) mod 1.                   (S12)
```

The primitive danger count and density at that digit depth are

```text
B_(r+1)=12(13^(r+1)+13(-1)^(r+1))/91,

B_(r+1)/phi(13^(r+1))
 =1/7+(-1)^(r+1)/(7*13^r).                           (S13)
```

The smallest density is `24/169`, at `r=2`. Since

```text
24/169-961/6930=3911/1171170>0,                      (S14)
```

THM-2137's base `A_0` floor alone cannot empty the deepest digit gate.
This does not obstruct the stronger image pump (M6), whose first-image
mass exceeds `1/7`; it says that advancing the private defect selector
requires its label/carry or a new uniform floor after owner deletion.

Full winding is genuinely necessary even for one terminal. Take

```text
H=1,       y=1/50,       q=5,       q'=18.
```

The two coefficients agree modulo thirteen and have the same current
integer window `{0,1}`, hence the same current chord. At the root

```text
x=(y+2)/13=101/650,
```

their integer windows differ because

```text
floor(5x)=0,                 floor(18x)=2.            (S15)
```

Thus the next chords are distinct. Modulo-thirteen chord data cannot update
the movie; `q mod 169`, or the equivalent full carry/winding coordinate, is
already forced at the first lift.

These statements hold away from the finite band boundaries and images of
the cover-null set. On the `13`-power torsion annuli they are exact: all
inequalities are strict, and a failure of (2) would thicken to an open
failure.

## 7. Exact audit and hostile controls

The companion uses two independent descriptions:

1. direct strict inequalities on residues modulo `2197`;
2. the root-image integer windows (6)--(12), phase by phase modulo `169`.

It verifies equality of all `1014` unit masks and all `78` depth-one masks
between the two descriptions.  It also checks:

```text
|U|=1450;
root-fibre histogram=(110 fibres of size 9, 46 of size 10);
all 2028 guard-unit normalizations;
129792 coefficient-normalization identities;
twelve nonzero depth-two unit parts (six sign classes), all empty on U;
all 3081 shallow pairs, including 78 diagonal pairs;
1013 distinct unit masks;
the sole sign-representative duplicate {1,1098}, both empty;
the exact unique minimum (30);
the event-loss decomposition (32)--(33);
the full conditional table digest
d31e5c874d8b5893ff33fc35095c18dcdf865c7f8296285c1b3441a8d8d679d9.
```

Three pairs show why the tempting repeated-maximum bound

```text
|R|-5 max_a |S_a intersection R|
```

is insufficient: its margins are `-22,0,-2` at shallow pairs
`(7,14),(14,61),(23,46)`.  Keeping the five fixed labelled movies and using
their five separate capacities repairs the estimate, with margins
`94,102,98` on those hostile rows.  The inherited one-active-mask direction
also reappears with minimum margin `30`, and a constructed five-mask union is
included as a positive-direction control.

Normal and optimized Python transcripts are byte-identical.  Reproduce with

```text
python3 04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py
python3 -O 04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py
```

## 8. Connection and loss ledger

The exact carrier is

```text
source:
  scalar 5+3 cover on a thirteen-root tower;
map:
  x -> 13x, with occupancy and sheets labelled by -Hk mod 13;
target:
  a cyclic movie of rooted, labelled mask--sheet incidences;
preserved:
  fibre coverage, owner multiplicity, event ties, coefficient labels,
  full winding, division support, quantitative image mass, and the two
  shallower deep-mask bits;
destroyed by the static chord quotient:
  endpoint identity through deletion, event order, winding, and deep bits;
needed sidecar:
  the labelled ownership current plus the simultaneous event word;
decisive test:
  the depth-three global incidence deficit (29).      (35)
```

THM-2199 places every primitive zero-Haar row in an explicit finite box, so
an exhaustive extension of (29) is possible in principle. Its bound is vast
and no such audit is performed. The rank theorem makes the remaining search
finite but does not replace coefficient winding by a smaller state.

THM-2201 rewrites the labelled sheet component as five Hasse jets in

```text
F_13[C_13]=F_13[epsilon]/(epsilon^13).
```

In that basis deck translation is triangular and endpoint insertion/deletion
is affine. The occupancy-thirteen failure of (M7a) is exactly the top term
`epsilon^12`. This is a faithful representation of the movie, not a rule for
deriving its event word; the guard, deep bits, and winding in (35) remain
load-bearing.

Stacking the successive chord matchings resembles a labelled
Brauer/matching cobordism.  That analogy is useful only at the carrier level.
Without chirality and over/under crossing data it is not a braid, and without
Reidemeister-compatible crossing data it is not a knot invariant.  The
winding and ownership sidecars are load-bearing; a knot or tournament
quotient that forgets them cannot recover (29).

The theorem identifies the recursive state, supplies measure and private-owner
transitions, and excludes the first unique-depth profile. By itself it does
not exclude `lambda_3>=3`; THM-2204/2205/2207 later close the depth-three
rows. It gives no uniform mass floor after owner deletion, does not bound unit
coefficient winding, classify the remaining `216` residue-length profiles,
or prove LRC(14). QED.
