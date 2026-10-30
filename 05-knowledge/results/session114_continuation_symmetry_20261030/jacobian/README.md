# One-coordinate symmetry and the first aliasing shell of Keller maps

**Session 114, 30 October 2026. PROVED special-class statements.** The planar
Jacobian conjecture remains open. These are elementary proofs in characteristic
zero, not claims of literature priority. We use no higher-dimensional Jacobian
claims from other repository files.

The previous session lifted **full-map** finite symmetry to continuous symmetry
under a degree bound, then forced linearity. Here one genuinely symmetric
coordinate already forces an arbitrary polynomial Keller companion to be
triangular. Full finite symmetry also remains harmless at the first aliasing
boundary: the first allowed nonlinear terms are shears, not a possible
counterexample.

## 1. One coordinate is enough; the other has no degree bound

Work over an algebraically closed field K of characteristic zero. Let
F=(P,Q) in K[x,y]^2 satisfy {P,Q}=P_x Q_y-P_y Q_x=c in K*. Fix coprime positive
integers p,q and the action

    (x,y) -> (lambda^p x, lambda^(-q) y).

**Theorem 1.** If Q has exact integer weight -q, that is,

    Q(lambda^p x,lambda^(-q)y)=lambda^(-q) Q(x,y)

for every lambda in K*, then

    Q=b y,                 P=a x+f(y),                 ab=c,

where a,b are nonzero scalars and f is any polynomial. No symmetry or degree
hypothesis is imposed on P. Conversely every such pair has the stated symmetry
of Q and determinant c.

If instead P has exact weight p, the symmetric conclusion is

    P=a x,                 Q=b y+g(x),                 ab=c.

**Proof.** Set t=x^q y^p. The nonnegative solutions of pi-qj=-q give
Q=yB(t). Write P as its finite direct sum of integer-weight components. A bracket
with Q shifts weight r to r-p, so only the weight-p part of P can contribute
to the weight-zero constant c. That part is xA(t). Consequently

    {xA(t),yB(t)} = AB+q t A'B+p t AB' = c.             (1)

In particular A and B are nonzero. If their degrees are alpha and beta, the
coefficient of t^(alpha+beta) is

    (1+q alpha+p beta) lead(A) lead(B).

It is nonzero in characteristic zero. Thus alpha=beta=0. Now Q=b y and the
full original Keller equation is b P_x=c. Integrating in x in characteristic
zero gives P=(c/b)x+f(y). The argument with P fixed uses the weight -q part of
Q and the same identity. This proves both assertions. QED.

This is not Reynolds averaging of a map. We select one coefficient channel
of an identity already satisfied by the original map. No assertion that the
projected pair inherits all the original equations is required.

Characteristic zero is essential: in characteristic ell>0, the pair
(x+x^ell,y) has determinant 1 and exact weight -q in its second coordinate,
but is not of the asserted form. The derivative/integration step would fail.

The inverse is explicit: for (X,Y)=(a x+f(y),b y),

    y=Y/b,                 x=(X-f(Y/b))/a.

## 2. A finite certificate using only the symmetric coordinate

Let zeta have order m. Suppose only

    Q(zeta^p x,zeta^(-q)y)=zeta^(-q) Q(x,y),

and deg Q <= d, where d>=1. The integer weight of a supported monomial lies
between -qd and pd. Its difference from -q lies in

    [-q(d-1), pd+q].

**Corollary 2.** If

    m > max(q(d-1), pd+q),                             (2)

then Q has exact weight -q, and Theorem 1 gives F=(a x+f(y),b y).
The degree of P is unrestricted. Neither full finite equivariance of F nor
a degree bound on F is needed.

For a prescribed P of degree at most e, the symmetric sufficient bound is

    m > max(qe+p, p(e-1)),                             (3)

giving F=(a x,b y+g(x)).

The bounds are sufficient no-alias bounds, not asserted optimal Keller
thresholds. A sharper finite support test simply asks whether every supported
monomial has the specified integer weight, rather than merely its residue.
Simultaneous tests with orders m_1,...,m_k in these same coordinates and weights
combine using lcm(m_1,...,m_k).

**Strict-bound hostile.** At p=q=1, d>=2 and m=d+1, Q=y+x^d is cyclic
semi-invariant of weight -1 and (x,Q) is Keller, but Q is not b y. Thus one
cannot replace the strict inequality in (2) uniformly by a weak inequality
and retain its specified triangular direction. This is an automorphism, not
a Jacobian counterexample.

## 3. Complete classification at the old sharp boundary

Return to weights (1,-w), w>=1. Here **full equivariance** means
F(zeta x,zeta^(-w)y)=(zeta P,zeta^(-w)Q) in these coordinates.

**Theorem 3.** Let D>=2 and m=wD+1. Every fully equivariant Keller map of degree
at most D is a triangular automorphism. More precisely:

- If w>1, F=(a x+alpha y^D,b y), with ab=c.
- If w=1, F is either (a x+alpha y^D,b y) or
  (a x,b y+beta x^D), with ab=c.

The diagonal case belongs to both alternatives. Thus the old linearity cutoff
is sharp, but its very first nonlinear alias cannot produce a Keller
counterexample.

**Proof for w>1.** We have

    m-w(D-1)=w+1>0,
    m-(D+w)=(w-1)(D-1)>0.

Corollary 2 applies to Q and forces Q=b y, P=a x+f(y). Equivariance requires
-wj=1 modulo m for every nonzero term y^j of f. For 0<=j<=D the quantity wj+1
ranges from 1 to m, so the only solution is j=D. This proves the claimed form.

**Balanced proof.** When w=1,m=D+1, the support condition gives exactly

    P=xA(t)+alpha y^D,     Q=yB(t)+beta x^D,     t=xy.

Indeed the only weights congruent to 1 within [-D,D] are 1 and -D; weight -D
is attained only by y^D. The statement for Q is symmetric. Direct calculation
splits the determinant into three distinct integer-weight channels:

    weight 0:    AB+t(A'B+AB')-alpha beta D^2 t^(D-1),
    weight -m:   -alpha D y^(D+1) B'(t),
    weight +m:   -beta D x^(D+1) A'(t).                 (4)

Hence alpha B'=beta A'=0. If both aliases are nonzero, A and B are constants,
and the nonzero term -alpha beta D^2 t^(D-1) contradicts a constant determinant,
since D>=2. If just alpha is nonzero, B is constant and nonzero, and the
weight-zero equation B(A+tA')=c forces A constant. The beta case is symmetric.
If both aliases vanish, (1) forces A and B constant. QED.

**Degree-one boundary.** For D=1, m=w+1, the two action weights are congruent.
The action is scalar, and all invertible linear maps are allowed. They need
not be triangular in the displayed coordinates. There are no translation
terms because m>=2. This is why Theorem 3 states D>=2.

## 4. The balanced shell extends one degree further

**Corollary 4.** For weights (1,-1), m>=3, every fully equivariant Keller map
of degree at most m is of one of the forms

    (a x+alpha y^(m-1),b y),  (a x,b y+beta x^(m-1)).    (5)

**Proof.** In the degree-m support triangle, weights congruent to 1 can only
be 1 or 1-m: weight 1+m is beyond the triangle. A monomial of weight 1-m is
y^(m-1)(xy)^j; its next degree after m-1 is m+1, so only j=0 is available.
The same holds for Q. Formula (4), with D=m-1>=2, proves the result. QED.

At m=3 the degree bound is sharp: the degree-four automorphism

    (x+y^2, y+(x+y^2)^2)

is fully equivariant modulo 3 but is not triangular in the displayed
coordinates. It has determinant 1. This does not assert sharpness at every m.
For every m>=3, composing the two allowed shears of exponent m-1 gives a
fully equivariant nontriangular automorphism of degree (m-1)^2.

## 5. Precisely what remains open

These criteria can certify that a map with actual, sufficiently strong
coordinate symmetry is invertible. They do not manufacture such symmetry in
an arbitrary Keller map. In particular, projecting a nonsymmetric coordinate
may change the constant determinant: the previous session's omitted-bracket
identity and the concurrent weighted two-shear resonance remain the obstruction.

A general counterexample would have to evade these certificates in every
applicable displayed coordinate system. No ability to find a suitable
coordinate change, classify all multi-alias supports, or solve JC2 is proved.

Sources within this research sequence:
- Session 113 `session113_symmetry_consumers_20261030/jacobian/README.md`, the
  full-map no-alias theorem and exact omitted bracket equation.
- Concurrent `phase_budget_shared_orbits_20261030/GEOMETRY_JC.md`, section 3,
  the weighted two-shear resonance criterion.
Both are prior results, not claims established by this packet's finite tests.

## 6. Reproduction

Run `python check_one_coordinate.py` and `python -O check_one_coordinate.py`.
All gates explicitly raise exceptions and remain enabled under optimization.
The script checks exact integer-weight support, the bracket channel projection,
the polynomial identity (1), boundary support and determinant formulas,
finite-certificate examples with arbitrarily larger companion degree, and
hostiles at the strict boundary and one degree past the balanced shell.
`check_results.json` records the deterministic output. The algebraic proofs
above establish all parameters; finite controls are regression tests.
