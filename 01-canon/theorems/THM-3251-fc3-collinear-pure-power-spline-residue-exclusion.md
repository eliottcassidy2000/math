---
id: THM-3251
title: "FC(3) collinear pure-power phases: B-spline residue blocks cannot cancel"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.  Let ell be
  a nonconstant algebraic affine polynomial on the coordinate two-simplex
  whose three vertex values are collinear in C.  For every integer d>=3 and
  algebraic A,h,C with A nonzero, q=A(ell-h)^d+C satisfies
  int_Delta exp(q)dA !=0 and !=1/2.  For three distinct ordered knots, the
  pushforward is a two-piece linear B-spline and its primitive weights obey
  `sum alpha_j z_j^2=-1`; the same two residue blocks and endpoint-collision
  obstruction as THM-3250 follow.  When exactly two knots agree, either
  residue block can vanish (both sharp cases occur), but they cannot both
  vanish because their sum has value 1/2 at s=0.  Every surviving block is
  independent from the endpoint exponentials, and the two surviving blocks,
  when both occur, are independent from each other because their residues
  are 1/d and 2/d.  Beukers Corollary 1.4 at the ordinary point s=1 excludes
  both forbidden values.  Together with THM-3250 this closes pure powers in
  every nonconstant algebraic affine coordinate, but not arbitrary
  one-variable polynomials, general cubics, nonalgebraic phases, or FC(3).
source: codex-2026-08-02-fc3-cubic-frontier
depends_on:
  - THM-3250-fc3-noncollinear-pure-power-turn-current-exclusion
related:
  - THM-3039-the-FC-n-exponential-period-bridge-forced-level-is-the-simplex-volume
  - THM-3202-fc3-real-affine-coordinate-quadratic-two-piece-spline-exclusion
  - THM-3203-fc3-complex-affine-coordinate-quadratic-cauchy-green-cycle-exclusion
external:
  - "Riesz--Markov--Kakutani representation theorem."
  - "F. Beukers, A refined version of the Siegel--Shidlovskii theorem, Annals of Mathematics 163 (2006), 369--379, Corollary 1.4."
script: 04-computation/fc3_collinear_pure_power_spline_residues_thm3251.py
output: 05-knowledge/results/fc3_collinear_pure_power_spline_residues_thm3251.out
script_sha256: 89082ea37498aaaf1daebd8cea552103e110834135d087443d882d39eba46c95
output_sha256: 1a823697d5a5fbe80acb8bd5e2e0cce502209c1b5de9832f0cdf50c12dc03daa
---

# THM-3251 -- collinear pure-power phases on the triangle

## 1. Statement and knot multiplicities

Let

```text
Delta={(u,v):u>=0, v>=0, u+v<=1},        area(Delta)=1/2.     (1)
```

Let `ell in Qbar[u,v]` be a nonconstant affine polynomial whose three
vertex values are collinear in `C`.  Fix

```text
d in Z, d>=3,        A,h,C in Qbar,        A!=0,              (2)
z=ell-h,             q=A z^d+C,
K(s)=int_Delta exp(sq)du dv.                                  (3)
```

Then

```text
K(1)!=0,                       K(1)!=1/2.                     (4)
```

There are exactly two nonconstant knot types:

1. three distinct vertex values, ordered along their common real affine
   line; or
2. one value repeated twice and one distinct value.

If all three values agree, `ell` is constant and is excluded by hypothesis.
Sections 2--4 handle the first type and sections 5--7 the second.  The proof
does not assume that the common line is the real axis.  Since all values are
algebraic, an algebraic affine reparameterization makes its line coordinate
real-valued; the formulas below are written directly in the original
complex knot coordinate so the phase remains the pure power `Az^d`.

## 2. Three distinct knots: the complex-line B-spline

Order the three distinct values along their line as `a,b,c`, and put

```text
p=b-a,             q_0=c-b,             L=c-a=p+q_0.         (5)
```

The three differences are nonzero and have one common complex direction.
The ordinary real-coordinate pushforward of area on `Delta` is the
two-piece linear B-spline of THM-3202.  Changing from its real parameter to
the complex coordinate `z` gives, for every entire `f`,

```text
int_Delta f(z)dA
 =1/L [ 1/p   integral_a^b (z-a)f(z)dz
       +1/q_0 integral_b^c (c-z)f(z)dz ].                    (6)
```

For `f=1`, the right side is `(p/2+q_0/2)/L=1/2`, so both the density and
the complex-line orientation are normalized.  Reversing the ordering swaps
`a,c`, reverses both line integrals, and leaves (6) unchanged.

Define the entire pure-power primitives

```text
J_(k,x)(s)=integral_0^x t^k exp(sAt^d)dt,       k=0,1,       (7)
```

and the second-divided-difference weights

```text
alpha_a=-1/(Lp),       alpha_b=1/(p q_0),
alpha_c=-1/(L q_0).                                       (8)
```

Direct algebra gives

```text
sum_x alpha_x=0,
sum_x alpha_x x=0,
sum_x alpha_x x^2=-1.                                      (9)
```

Expanding the two line integrals in (6) therefore yields

```text
K(s)=e^(Cs)(G_0(s)+G_1(s)),                                 (10)
G_1=sum_(x in {a,b,c}) alpha_x J_(1,x),
G_0=-sum_(x in {a,b,c}) alpha_x x J_(0,x).                   (11)
```

The sign in the last identity is load-bearing.  At `s=0`, (9) gives
`G_0(0)+G_1(0)=1/2`, as required.

## 3. Distinct-knot source and collision audit

As in THM-3250, differentiation of `t^(k+1)exp(sAt^d)` gives

```text
dsJ_(k,x)'+(k+1)J_(k,x)=x^(k+1)exp(Ax^d s).                 (12)
```

Thus

```text
dsG_1'+2G_1=S,
dsG_0'+ G_0=-S,                                               (13)
S=sum_x alpha_x x^2 exp(Ax^d s),
S(0)=-1.                                                       (14)
```

Group equal endpoint powers by `xi=Ax^d`, with source coefficient

```text
tau_xi=sum_(x:Ax^d=xi)alpha_xx^2.                            (15)
```

Equation (14) says `sum_xi tau_xi=-1`, so some grouped source is nonzero.
The group `xi=0` has zero source because every knot in it is zero.  Hence
some group satisfies

```text
tau_xi!=0,                         xi!=0.                     (16)
```

This includes opposite-endpoint collisions for even `d`, root-of-unity
collisions for odd `d`, and a zero knot.  No genericity of the endpoint
powers is used.

## 4. Distinct-knot functional independence

Let `mathcal E` be the `Qbar(s)`-span of the distinct functions among the
endpoint exponentials, `exp(-Cs)`, and `1`.  If `G_k`, `k=0,1`, belonged to
`mathcal E`, comparison at the source class (16) after applying (13) would
give a rational solution of

```text
dsR'+(k+1+d xi s)R=epsilon tau_xi,                            (17)
```

where `epsilon` is `-1` or `1`.  There is none.  A nonzero finite pole gains
one order under the derivative; a pole of order `n` at zero has multiplier
`k+1-dn`, which is never zero because `1<=k+1<=2<d`; and a polynomial gains
one degree through `d xi sR` because `xi!=0`.

Both blocks therefore survive modulo `mathcal E`.  They are jointly
independent: modulo `mathcal E`, the operator `T=ds d/ds` has eigenvalues
`-1,-2` on `G_0,G_1`.  A relation `[G_0]+R[G_1]=0` would give `dsR'=R`,
requiring the nonintegral order `1/d` at zero.  This proves functional
independence in the three-distinct-knot case.

## 5. One doubled knot: the triangular marginal

Now suppose the repeated value is `a` and the value at the remaining vertex
is `c!=a`.  Put `L=c-a`.  The barycentric coordinate `t` of the unique
vertex has marginal density `1-t` on `[0,1]`.  Therefore

```text
int_Delta f(z)dA=1/L^2 integral_a^c (c-z)f(z)dz.             (18)
```

Its mass is `1/2`, and it is invariant under relabeling the two equal
vertices.  With the primitives (7), equation (18) becomes

```text
K(s)=e^(Cs)(G_0(s)+G_1(s)),                                  (19)
G_1=[J_(1,a)-J_(1,c)]/L^2,
G_0=c[J_(0,c)-J_(0,a)]/L^2.                                  (20)
```

The residue equations are now

```text
dsG_1'+2G_1=S_1,
S_1=[a^2 exp(Aa^d s)-c^2 exp(Ac^d s)]/L^2,                  (21)

dsG_0'+G_0=S_0,
S_0=c[c exp(Ac^d s)-a exp(Aa^d s)]/L^2.                     (22)
```

The two sources are no longer opposites.  Requiring both to be nonzero
would be false: for symmetric endpoints and even `d`, `S_1=0`; if the
unique value is zero, `S_0=0`.

## 6. Why one surviving block is enough

If `S_k=0` identically, (21) or (22) is a homogeneous half-residue equation.
The corresponding `G_k` is entire; evaluation at zero and then its Taylor
recurrence force `G_k=0` identically.

The two sources cannot both vanish.  Otherwise both blocks would vanish,
and (19) would give `K(0)=0`, contradicting

```text
K(0)=area(Delta)=1/2.                                        (23)
```

Thus at least one block is nonzero.  In every nonzero grouped source, a
nonzero coefficient occurs at an exponent `xi!=0`: every coefficient in
(21)--(22) attached to a zero endpoint contains that endpoint as a factor
and vanishes.  The rational obstruction (17) therefore proves that each
nonzero block survives modulo `mathcal E`.

If both blocks survive, the eigenvalue argument following (17) proves they
are jointly independent.  If exactly one survives, independence is already
proved.  Since every zero-source block is the zero function, in all cases

```text
F:=G_0+G_1                                                   (24)
```

is a nonzero functional class independent from `mathcal E`.  This is why
one surviving residue block suffices for the period itself; no vanished
block is silently retained in the value relation.

## 7. E-functions, Beukers, and the two values

The primitive expansion is

```text
J_(k,x)(s)=sum_(n>=0)
 A^n x^(dn+k+1)s^n/((dn+k+1)n!).                             (25)
```

It has algebraic coefficients, exponential conjugate and denominator
growth, and the holonomic equation (12), so it is an E-function.  Every
surviving block and every endpoint exponential is therefore an E-function.
Equations (13), (21), and (22) give a first-order system with denominator
only `s`; `s=1` is ordinary.

Beukers Corollary 1.4 transfers the functional independence proved in
sections 4 and 6 to algebraic linear independence of the corresponding
values at `s=1`.  In either knot type, equations (10) or (19) give

```text
K(1)=e^C F(1).                                                (26)
```

If `K(1)=0`, then `F(1)=0`, an impossible value relation among the surviving
blocks.  If `K(1)=1/2`, then

```text
F(1)-(1/2)exp(-C)=0,                                         (27)
```

an impossible relation after adjoining `exp(-Cs)`.  Coincidence of that
exponential with an endpoint exponential or with `1` only merges basis
elements; it does not remove the nonzero coefficient of `F`.  This proves
(4).

## 8. Combined consequence and sharp scope

THM-3250 handles noncollinear vertex images; the present theorem handles
every nonconstant collinear image.  Hence, together:

```text
For every nonconstant ell in Qbar[u,v], every d>=3, and
A,h,C in Qbar with A!=0,

int_Delta exp(A(ell-h)^d+C)dA is neither 0 nor 1/2.           (28)
```

This is the complete algebraic affine-coordinate **pure-power** branch.  It
does not cover:

* a general depressed cubic `A(ell-h)^3+B(ell-h)+C` with `B!=0`;
* sums of powers in different affine coordinates;
* nonalgebraic coefficients;
* the projective-leading-form reduction required for full `FC(3)`; or
* arbitrary homogeneous Factorial-Conjecture inputs.

At `d=2`, the second residue is integral and can become elementary; the
all-quadratic result is instead THM-3203.  The constant-coordinate boundary
is excluded here because it has no knot geometry and reduces directly to a
constant exponential.

The connection contract is

```text
source:    area measure on a triangle with collinear affine image
target:    a two-piece or doubled-knot linear B-spline
map:       affine pushforward, then two pure-power primitives
preserves: total period and the forbidden values 0 and 1/2
loses:     transverse level-set position
sidecar:   ordered knot multiplicity and residue-block survival
hostiles:  endpoint-power collision, zero knot, and either one-block death. (29)
```

## 9. Exact verification

Run

```text
python3 04-computation/fc3_collinear_pure_power_spline_residues_thm3251.py
python3 -O 04-computation/fc3_collinear_pure_power_spline_residues_thm3251.py
```

The exact companion checks:

* 99 direct simplex / B-spline pushforwards, including orientation reversal;
* the symbolic weight moments `(0,0,-1)`;
* 216 direct pure-power moment / primitive identities;
* opposite, cubic-root, and zero endpoint collisions for distinct knots;
* the sharp doubled-knot controls where `G_1` alone or `G_0` alone dies; and
* the nonintegral residue boundary for `3<=d<=40`, with `d=2` planted as the
  failure control.

Normal and optimized transcripts are byte-identical.  The verifier checks
the algebra and normalizations; the E-function value step is the cited
Beukers theorem.

**End of proof.**
