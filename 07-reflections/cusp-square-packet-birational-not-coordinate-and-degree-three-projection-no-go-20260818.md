# The cusp-square packet is birational, not a coordinate embedding

**Status: VERIFIED-EXACT + FINITE-EXACT.  For the explicit `U_*` in
THM-3556, the packet subring is proper but has the full source function
field.  Neither source coordinate is polynomial in the packet, while both
are rational in it.  No packet-polynomial projection of packet degree at
most three has nonzero constant Jacobian.  Degree four and higher remain
open.  THM-3556 is unchanged, and no `JC(2)` conclusion is claimed.**

The object is the explicit specialization in THM-3556 --
`cusp-square-packet-marked-root-kummer-owner`:

```text
U=1+y-y^2/2-(3/2)vy(y-3),
T=y^2-6vU,
S=y^3-9vUy,
L=v^2(8vU-y^2),
Z=(L,T,U,S).
```

The exact companion is
`04-computation/cusp_square_packet_subring_projection_audit_20260818.py`;
its matching transcript is
`05-knowledge/results/cusp_square_packet_subring_projection_audit_20260818.out`.

## Inheritance pass and concept board

The closest proved mechanism is THM-3556 itself: the six Jacobian minors of
`Z` generate `(1)` in `Q[v,y]`, but no constant combination is one.  Its new
dual cubic

```text
Y^3-3TY+2S
```

is the correct starting representation for source recovery.

The canonical hostile is the gap already isolated there: arbitrary minor
Bezout coefficients are not automatically packet-descending coefficients,
and packet-descending coefficients are not automatically the Pluecker
coordinates of one exact decomposable form `dA wedge dB`.  The corrected near
miss is MISTAKE-421: a full etale/immersive map does not make a selected
coordinate eliminant primitive or injective.  MISTAKE-424 gives the parallel
image-typing warning: geometric fibres, rational-point fibres and a dense
scheme-theoretic image are different assertions.

The least-used relevant sidecar is not another discriminant.  It is the
literal defining equation for the specialized `U`, transported through
`6vU=y^2-T`.  The live board was:

1. packet subring versus packet function field;
2. the cubic/quadratic common-root eliminant and its first subresultant;
3. exceptional packet fibres on the subresultant divisor;
4. arbitrary six-minor Bezout coefficients;
5. exact decomposable projection pairs in a bounded packet-degree filtration.

The niche lane (subring versus field) overtook the initial projection search:
it supplied both the rational inverse and the hostile that rules out
polynomial recovery.  The wildcard was the `y=0` slice of the exceptional
divisor; it produced an exact two-point collision over a cubic number field.

## 1. The extra quadratic makes the visible cubic birational

Every packet point satisfies the visible cubic

```text
f(Y)=Y^3-3TY+2S,                   f(y)=0.             (1)
```

Using `vU=(y^2-T)/6` in the displayed formula for `U` gives a second relation

```text
r(Y)=(U+T)Y^2-(2U+S+3T)Y+2U^2-2U+3S,   r(y)=0.       (2)
```

Put

```text
a=U+T,
b=-(2U+S+3T),
c=2U^2-2U+3S.
```

The first Euclidean remainder is the exact identity

```text
a^2 f(Y)+(b-aY)r(Y)=D Y+C,                            (3)
```

where

```text
D = S^2+3ST+SU-3T^3-6T^2U+9T^2-5TU^2+14TU-2U^3+6U^2,
C = -3S^2+2ST^2+4STU-9ST-4SU-6TU^2+6TU-4U^3+4U^2.
```

Direct pullback gives `D(Z)y+C(Z)=0`, while `D(Z)(0,0)=4`.  Therefore `D`
does not vanish in the packet function field and

```text
y = -C/D.                                              (4)
```

The identity `y^2-T=6vU` then gives

```text
v = (C^2-TD^2)/(6UD^2).                               (5)
```

Equations (4)--(5) prove

```text
Q(L,T,U,S)=Q(v,y),
[Q(v,y):Q(L,T,U,S)]=1.                                (6)
```

Thus the visible cubic's degree three is not the packet field degree.  The
specialized `U` cuts it with a quadratic, and their generically linear gcd
selects the marked root rationally.

As an elimination control, exact lexicographic Groebner computation over
`Q` for `(f,r)` in `Q[Y,T,U,S]` has 13 basis elements and exactly one
`Y`-free generator.  It agrees, up to a rational unit, with the 34-term
resultant `Res_Y(f,r)` printed in the transcript.  Direct substitution makes
that resultant zero.  This is a polynomial-identity computation, not a
rational-point inference.

## 2. The subring is nevertheless proper

Let `alpha` satisfy the irreducible cubic

```text
p(alpha)=324alpha^3+54alpha^2-27alpha+2=0.             (7)
```

In `K=Q(alpha)`, define

```text
P=(alpha,0),
Q=(-2alpha, 2(1-9alpha)/(1-6alpha)).                   (8)
```

Exact remainder reduction modulo `p` gives

```text
Z(P)=Z(Q)=(8alpha^3,-6alpha,1,0).                      (9)
```

The denominator in (8) is nonzero in `K`, and both the `v`- and
`y`-coordinates differ.  Hence `Z` is not geometrically injective.  Any
polynomial expression for either `v` or `y` in `L,T,U,S` would take the same
value at `P` and `Q`, contradicting (8).  Consequently

```text
Q[L,T,U,S] proper-subset Q[v,y],                       (10)
v,y notin Q[L,T,U,S].                                  (11)
```

Together, (6) and (10)--(11) give the exact answer: both source coordinates
are rational packet functions, neither is a polynomial packet function, and
the exhibited failure lies on the exceptional inverse denominator `D=0`.
The collision is geometric over a cubic extension; no rational collision is
claimed or needed for the subring obstruction, and no classification of all
exceptional fibres is asserted.

The connection contract is now explicit.  The packet map preserves the
generic function field and hence generic source coordinates.  It destroys
point separation on `D=0`.  A marked-sheet or normalization sidecar restores
that information.  The cheapest decisive test is exactly the collision (8).

## 3. An arbitrary minor certificate is not a projection

Write the six packet minors in coordinate order as

```text
(M_LT,M_LU,M_LS,M_TU,M_TS,M_US).
```

An exact source-degree search for arbitrary `q_ij(v,y)` has rank ledgers

```text
coefficient degree     matrix       rank / augmented rank
0                      26 x 6       6 / 7
1                      37 x 18      18 / 19
2                      49 x 36      36 / 37
3                      62 x 60      50 / 50.           (12)
```

The degree-three solve produces a displayed five-term coefficient tuple in
the transcript with

```text
sum q_ij M_ij=1.                                      (13)
```

This is a positive control for the known unit minor ideal.  It is not a
legal projection certificate:

- the first five nonzero `q_ij` take different values at the equal-packet
  points (8), so this tuple does not descend through `Q[L,T,U,S]`;
- its Pluecker expression
  `q_LT q_US-q_LU q_TS+q_LS q_TU` is nonzero already at `(v,y)=(0,0)`.

Thus (13) is neither a descended coefficient tuple nor a decomposable
bivector.  Integrability is also not inferred from (13).

## 4. Exact projection no-go through packet degree three

For packet degree `d`, let `W_d` be the pullback to `Q[v,y]` of packet
polynomials of total degree at most `d`, modulo constants.  If
`A,B in W_d`, then

```text
Jac(A,B)=sum_(i<j) (a_i b_j-a_j b_i) Jac(w_i,w_j).     (14)
```

The coefficient bivector in (14) is decomposable.  The computation makes a
strict relaxation: it allows an arbitrary coefficient on every pairwise
Jacobian `Jac(w_i,w_j)`.  Failure in this larger linear span therefore proves
failure for every legal exact decomposable pair without having to solve the
Pluecker equations.

The exact `Q` ranks are

```text
d   formal / pullback dim   pairs   coefficient matrix   rank / augmented
1          4 / 4              6          26 x 6               6 / 7
2         14 / 14            91         139 x 91             67 / 68
3         34 / 33           528         336 x 528           187 / 188.       (15)
```

The single degree-three pullback dependence is the first packet-filtration
appearance of the cusp identity.  At every row of (15), adjoining the target
constant `1` raises rank by one.  Hence `1` is outside even the full arbitrary
bivector span.  Two prime-field replays, modulo `1000003` and `1000033`, give
the same ranks.  At degree two, the script additionally constructs and
verifies a 68-support dual functional that annihilates all 91 Jacobian
columns and evaluates to one on the constant target; its canonical digest is
printed in the transcript.

Since any nonzero constant can be rescaled to one, (15) proves:

```text
There are no A(Z),B(Z) of packet degree <=3
with Jac_(v,y)(A(Z),B(Z)) in Q*.                       (16)
```

No injectivity test for a projection is therefore applicable in this box.
The packet itself is generically recoverable by (4)--(5) but fails global
injectivity by (8)--(9).

## 5. Boundary and reusable lesson

Direct progress is (6), (10)--(11), and the bounded no-go (16).  The niche
progress is the exact exceptional collision.  The reusable method is to
place legal `dA wedge dB` forms inside the full pairwise-Jacobian span and
test that linear superset first; a rank obstruction there is stronger and
cheaper than a Pluecker-ideal solve.

The honest remaining frontier is packet degree at least four.  Failure
through degree three does not exclude a higher-degree nonlinear projection,
another packet, or a planar Keller counterexample.  Equality of function
fields also does not turn this packet into a coordinate embedding: the
exceptional divisor has already recorded the lost sheet.

Reproduce from the repository root with

```text
python3 04-computation/cusp_square_packet_subring_projection_audit_20260818.py
python3 -O 04-computation/cusp_square_packet_subring_projection_audit_20260818.py
```

Normal and optimized transcripts are byte-identical after newline
normalization.
