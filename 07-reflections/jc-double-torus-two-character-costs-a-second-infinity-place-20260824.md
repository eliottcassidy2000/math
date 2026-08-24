# Two cubic characters are cheap; keeping one infinity place is the expensive part

**Status (2026-08-24): CURRENT SCOUT SYNTHESIS WITH PROVED CANONICAL
SUCCESSORS.**  Sections 2--4 retain the original noncanonical exact near miss,
replayed by
`04-computation/jc2_double_torus_two_character_two_infinity_scout.py`.
The Mordell--Weil and localization deductions use standard elliptic-surface
arguments spelled out here.  The subsequent canon now sharpens its boundary:
THM-3942 proves the complete affine-linear whole-factor obstruction;
THM-3944 proves the repeated-factor conductor collision; and THM-3947 proves
the full scalar-weighted repeated-square trichotomy.  THM-3946 now proves the
complete affine one-factor split and has passed independent hostile audit;
THM-3949 proves the all-degree coprime one-variable Newton obstruction; and
THM-3950 classifies every nonconstant `A1` pullback and exposes its fixed
equianharmonic residual.
This reflection reserves no theorem ID; its affine deformation is now an
independent precursor to the `m=1` slice of THM-3946.

## 1. The rank-two gate before searching

THM-3937 explains one failed design language completely.  Its generic affine
quadratic-resolvent curve is

```text
C=E minus {I_+,I_-},
```

so, after choosing `I_+` as elliptic origin,

```text
Pic(C)=E(K)/<I_--I_+>.                                  (1)
```

This gives an immediate Smith-normal-form gate.  If `E(K)` is torsion-free,
then quotienting a free abelian group of any rank by one vector produces at
most one cyclic torsion factor.  Therefore

```text
dim_F3 Pic(C)[3] <= 1                                    (2)
```

even if the Mordell--Weil rank is two, ten, or larger.  Merely raising rank is
not the missing move.  A second cubic direction requires at least one of:

1. intrinsic rational `E(K)[3]` in addition to a three-divisible boundary
   difference;
2. more deleted boundary sections and an honest rank-two relation matrix;
3. nonprincipal vertical classes that survive localization; or
4. a higher-genus generic curve whose Jacobian has two rational three-lines.

The cheapest option is the first.  It can be forced by giving the same branch
two different torus decompositions.

## 2. An exact double-torus identity

Over an algebraically closed field `k` of characteristic zero, put

```text
p_0=X,                         p_1=X+t,
q_0=3X^2+3tX+t^2-t,
q_1=3X^2+3tX+t^2+t.                                  (3)
```

Then the difference-of-cubes factorization gives the exact identity

```text
H=q_0^2-4p_0^3=q_1^2-4p_1^3.                           (4)
```

Indeed `q_1-q_0=2t`, `q_1+q_0=2(3X^2+3tX+t^2)`, and

```text
(X+t)^3-X^3=t(3X^2+3tX+t^2).
```

Thus the normal double cover

```text
B=k[t,X,W]/(W^2-H)                                      (5)
```

has two Cardano radicands `q_0+W` and `q_1+W`, supported over the distinct
linear divisors `X=0` and `X+t=0`.

The common branch polynomial is the irreducible quartic

```text
H=9X^4+(18t-4)X^3+(15t^2-6t)X^2
  +(6t^3-6t^2)X+t^4-2t^3+t^2.                          (6)
```

Its discriminant as a quartic in `X` is

```text
disc_X(H)=6912 t^4(t-1)^3(t+1)^3.                       (7)
```

The affine singular support of `H=0` is exactly

```text
(t,X)=(0,0),(1,0),(-1,1).                               (8)
```

The eliminated singular ideal has Groebner basis

```text
12X^2+2t^4+t^3-8t^2+5t,
12Xt-t^3+6t^2-5t,
t(t-1)^2(t+1)^2.                                        (9)
```

Hence `(5)` has only three isolated singularities and is normal.  Every
closed fibre over `t` is integral.  Outside `0,+/-1`, equation `(7)` rules
out a square.  At the three exceptional values,

```text
H_0=X^3(9X-4),
H_1=X^2(9X^2+14X+9),
H_-1=(X-1)^2(9X^2-4X+4),                               (10)
```

and both residual quadratics have discriminant `-128`, so none is a square.
Thus every vertical prime is the single principal fibre and Weil
localization gives

```text
Cl(B)=Pic(C),                                            (11)
```

where `C` is the smooth generic affine quartic over `K=k(t)`.

## 3. The branch is rational, but it has exactly two ends

The first torus presentation gives a rational normalization.  With parameter
`s`, set

```text
h=2(1-s)/(s^2+3),
t=(s-1)^3(s+3)/(s^2+3)^2,
X=4(s-1)^2/(s^2+3)^2=h^2.                               (12)
```

Then exactly

```text
q_0(t,X)=2h^3,                  H(t,X)=0.                (13)
```

This map is birational.  On the branch recover

```text
h=q_0/(2X),
v=(2t-1+3h^2)/(h-1),
s=(v-1)/h.                                               (14)
```

Homogenizing `(12)` gives the basepoint-free degree-four map

```text
[S:R] |-> [
  (S-R)^3(S+3R) : 4(S-R)^2R^2 : (S^2+3R^2)^2
].                                                       (15)
```

The inverse image of the target infinity line is therefore

```text
S^2+3R^2=0,                                              (16)
```

two distinct points.  Neither cancels from the other coordinates.  Hence the
affine normalization of `(6)` is

```text
P1 minus {two points} ~= G_m,                            (17)
```

not `A1`.  This is not a presentation artefact: `p_0 o nu=h^2`, and `h`
has a pole at both points in `(16)`.  Any affine target realization retaining
this polynomial torus coefficient must send both points to infinity.

The conic hidden in `(12)` makes the cost transparent.  Solving
`q_0(h^2,t)=2h^3` gives

```text
t^2+(3h^2-1)t+3h^4-2h^3=0,
disc_t=-(h-1)^3(3h+1).                                  (18)
```

After removing the square `(h-1)^2`, the normalization is

```text
v^2=-(h-1)(3h+1),                                       (19)
```

and the degree-two function `h` has two poles.  The same conic, up to scaling
and sign, occurs throughout the normalized translated-coefficient ansatz
`p_1=p_0+c`: changing which complementary cube factor is assigned to
`q_1-q_0` only exchanges the two ends.  This is a sharp obstruction for that
ansatz, not a claim about every possible double-torus curve.

## 4. Why this near miss really has two cubic characters

Let `I_+` and `I_-` be the two `K`-rational quartic infinities, labeled by

```text
W/X^2 -> +3 and -3,
```

and choose `I_+` as origin.  Let `P_0` be the point on `X=0` where
`q_0+W=0`, and `P_1` the point on `X+t=0` where `q_1+W=0`.  Both radicands
have the same pole divisor: order two at `I_+` and order one at `I_-`.
Consequently

```text
div(q_0+W)=3P_0-2I_+-I_-,
div(q_1+W)=3P_1-2I_+-I_-.                               (20)
```

In the elliptic group this says

```text
3P_0=I_--I_+=3P_1,
T=P_0-P_1 in E(K)[3].                                   (21)
```

The point `T` is nonzero because `P_0` and `P_1` have distinct generic
`X`-coordinates.  Thus the second torus decomposition creates intrinsic
rational three-torsion; it is not merely a second formula for the same affine
class.

The binary-quartic invariants of `(6)` are

```text
I=9t^2(t^2+8),
J=54t^2(t^4-20t^2-8),
4I^3-J^2=186624 t^4(t-1)^3(t+1)^3.                       (22)
```

The relatively minimal rational elliptic surface has fibres

```text
IV at 0,             I3 at +1 and -1,             I2 at infinity. (23)
```

Their Euler numbers sum to twelve and their root ranks sum to
`2+2+2+1=7`, so Shioda--Tate gives rank one.  Torsion injects into the
component group `Z/3` of the additive `IV` fibre; `(21)` supplies a nonzero
element, hence

```text
E(K)=Z G direct_sum Z/3<T>                              (24)
```

for some free generator `G`.

Put `R=I_--I_+=3P_0`.  The two infinity points are distinct, and the torsion
in `(24)` has exponent three, so `P_0` cannot be torsion.  Write
`P_0=mG+aT` with `m!=0`.  Equations `(1),(11),(21),(24)` now give

```text
Cl(B)=Pic(C)=E(K)/<R>
             ~= Z/(3m) direct_sum Z/3,

Cl(B)[3] ~= (Z/3)^2.                                    (25)
```

This conclusion does not require knowing whether `m=1`.  One three-line is
the boundary-divisibility class of `P_0`; the other is the intrinsic torsion
line `T` measuring the difference between the two torus decompositions.

Finally the resolvent units are still scalar.  A generic-curve unit would
have divisor supported on `I_+,I_-`, hence make a nonzero multiple of `R`
principal.  But `R=3P_0` has infinite order.  Thus `C^*=K^*`, and the complete
principal vertical ledger in `(10)-(11)` reduces this to

```text
B^*=k^*.                                                 (26)
```

So the near miss satisfies scalar units and two independent cubic characters
at once. Among those and the one-place invoice, it fails exactly the last,
and minimally: `G_m` instead of `A1`. It has not discharged the separate
descent, nonmonogenic-order, source, or Keller/Jelonek gates.

## 5. What the experiment changes

The live concept board after THM-3937 is now:

```text
single torus decomposition  -> boundary divisibility -> one C3 line;
second torus decomposition  -> intrinsic E[3]        -> second C3 line;
translated cube factors     -> conic normalization    -> two infinity ends;
one affine branch place     -> h has one pole          -> blocks this ansatz.
```

THM-3942 makes the translated example part of an exact affine-linear
classification.  For independent affine-linear `p_0,p_1`, whole-factor UFD
allocation has only a `1|2` conic with two ends and a `0|3` Fermat cubic with
three.  Both retain two independent characters.  Thus the failure is the
place divisor, not character scarcity, and any nonlinear `A1` escape must
split a factor internally or exploit gcd/multiplicity overlap.

THM-3944 and THM-3947 then audit the cheapest internal split.  For
`p_1-p_0=alpha G^2`, the generic scalar weighting produces three distinct
smooth one-place parabolas: componentwise geometry improves, but the full
branch is reducible.  The only two scalar collisions double `p_0` or `p_1`
and recover the nonnormal square-conductor geometry.  On that conductor seam
one Cardano row is an exact cube, while the other survives only as the class
`(2,1)` on the original order's `G_m^2` regular locus and ramifies across the
normalization.  A complementary exact affine slice in
`jc2_double_torus_nonlinear_balanced_partial_split_scout.py` separates the
factor as `t(t-c)`: `c!=0` gives an irreducible two-ended conic, whereas
`c=0` is precisely the doubled-conductor collision.

THM-3946 identifies that scout as its equality-slope `m=1` slice and closes
every affine slope: the other exceptional slope `m=-omega^2` also has two
ends, while all remaining coprime affine rows have three.  Coincident zeros
are reducible by THM-3947, exactly one constant factor returns to THM-3942,
and the MISTAKE-472-repaired two-constant boundary is a reducible cubic in
`P`.  THM-3949 removes the affine-degree restriction for every coprime
nonassociate one-variable pair.  Its last Newton row is blocked by a residual
cubic with zero quadratic coefficient and nonzero constant, so it cannot have
one triple direction.  Thus the one-variable one-factor lane is closed in the
standard chart, though arbitrary target lines are not.

THM-3950 changes the board again. On any nondegenerate `A1` pullback, the two
cusp parameters `r,s` have a reduced ratio whose denominator debt factors
exactly through `S+omega^2R` and `S-omega R`. The assigned-factor ratio is a
degree-three rational map with branch values `{0,1,-omega,infinity}`--the two
scalar collision seams of THM-3947 together with the zero/infinity endpoints.
Its `S3` closure is a fixed `j=0` elliptic cover. The explicit degree-one-ratio
row simultaneously has a reduced one-place graph component, a normal
quadratic surface, and two independent extendable Cardano characters. What it
does not have is a one-place full discriminant: an irreducible genus-one
residual remains, and no Keller map or same-field cubic atlas is supplied.

The next search should not maximize Mordell--Weil rank blindly.  The cheapest
decisive experiments are:

1. construct or obstruct the same-field cubic atlas and full source-boundary
   incidence of the THM-3950 survivor, including its genus-one residual;
2. search genuinely bivariate factors and simultaneous internal splits in two
   or three cube factors for an irreducible reduced full branch;
3. compute the generic three-parabola quadratic normalization and its complete
   class/Cardano lattice before deciding whether recombination is possible;
4. replace the line boundary in the multiple-torus sextic families by an
   explicitly normalized non-line or birational boundary;
5. keep one torus structure but move to a higher-genus generic fibre and
   search for a second Jacobian three-line not arising from a translated cube
   factor; or
6. create at least three deleted boundary sections, or a controlled split
   vertical fibre, and compute the full localization relation matrix.

For each lane the cheapest hostile gate is the normalization pole divisor.
The present sequence shows why this must be checked before any expensive
class-group computation: the desired second character appeared immediately;
whole factors forced extra ends; separated internal factors retained two
ends; the first one-place collision paid nonreduced conductor debt; and the
first reduced `A1`/normal/two-character survivor paid a genus-one companion.
A survivor must control the full discriminant and same-field cubic/source
attachment, not merely pass those componentwise gates.

## Reproduction

```bash
python3 04-computation/jc2_double_torus_two_character_two_infinity_scout.py
python3 -O 04-computation/jc2_double_torus_two_character_two_infinity_scout.py
python3 04-computation/jc2_double_torus_nonlinear_balanced_partial_split_scout.py
python3 -O 04-computation/jc2_double_torus_nonlinear_balanced_partial_split_scout.py
```

Each normal/optimized pair must byte-match its corresponding frozen output in
`05-knowledge/results/` after LF normalization.
