# Rational six-cycles, lost phase, and changing-generator Pell growth

2026-10-15. **PROVED elementary:** the rational neutral-cycle lemma,
fixed-parameter preperiodic classifications, the exact-period interpretation
of rational sixth dynatomic points, and the changing-generator identities.
**CITED / CONDITIONAL:** Stoll's general exclusion of rational six-cycles.
No general rational six-cycle exclusion or Collatz convergence theorem is
proved here. [Exact controls](quadratic_cycle_phase_20261015b.json) and
[script](../../04-computation/experiments/quadratic_cycle_phase_20261015b.py).

## 1. Inheritance: three different meanings of six

[THM-4139, rational-three-cycle-order-six-lift-and-horizontal-carrier](../../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md)
already proves the complete rational preperiodic graph at c=-29/16 and
constructs an order-six linear lift of a projective three-cycle.
[THM-4146, the horizontal-divisor/fibre firewall](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md)
keeps that lift distinct from a scalar polynomial orbit.
[Arithmetic braids, geometry section 3](arithmetic_braids_20260917_geometry.md)
already identifies the irreducible parabolic three-cycle at c=-7/4.
We reuse these stronger statements rather than reclassify a few sample points.

The live objects are a scalar polynomial orbit, a marked projective cycle,
a linear lift, and a quadratic observable of a linear state. A period in
one is not automatically a period in another. The missing coordinates are
respectively the parameter, marked point, sign sheet, and source vector.

## 2. Rational cycles cannot be neutral except at two fixed points

**Lemma R1 (PROVED).** If a rational point of exact period d for
f_c(x)=x^2+c, c rational, has a multiplier which is a root of unity, then
d=1. The only possibilities are (x,c)=(1/2,1/4), multiplier1, and
(-1/2,-3/4), multiplier-1.

**Proof.** The multiplier is rational, so a root of unity means +1 or -1.
First, all rational points in a periodic orbit have the same reduced
positive denominator b. Indeed, at a prime p with v_p(c)<0, a bounded
orbit must have v_p(x)=v_p(c)/2: a smaller valuation doubles and escapes;
a larger valuation maps to v_p(c), then escapes. If v_p(c)>=0, a negative
valuation of x similarly doubles and escapes. This proves the assertion
and also that the denominator of c is b^2.

Write the cycle points a_i/b in lowest terms. Its multiplier is
2^d product(a_i)/b^d. Any odd prime dividing b has strictly negative
valuation in that multiplier, so b is a power of2. If b=1 the multiplier
has positive 2-valuation and cannot be a root of unity. If b=2^h, its
2-valuation is d(1-h), hence h=1. The multiplier is now product(a_i)=+/-1,
so each a_i is +1 or -1. But f(1/2)=f(-1/2), and a map is injective on
a periodic cycle. Thus the cycle has one point. Substitution gives the
two stated pairs. This argument covers every rational parameter. ∎

The denominator condition applies to preperiodic points too, since their
forward orbits are bounded. For c=-a<0, the real escape bound is
|x| <= (1+sqrt(1+4a))/2. Combining these two bounds gives complete finite
candidate sets at the two requested parameters:

| Parameter c | Every rational preperiodic point | Unique rational cycle | Multiplier |
|---|---|---|---|
| -7/4 | +/-1/2, +/-3/2 | 1/2 -> -3/2 -> 1/2 | -3 |
| -29/16 | +/-1/4, +/-3/4, +/-5/4, +/-7/4 | -7/4 -> 5/4 -> -1/4 -> -7/4 | 35/8 |

The bounds restrict the numerators to the listed odd integers; evaluating
the map verifies that every candidate really is preperiodic. Consequently
neither parameter has a rational exact six-cycle. The second row is
inherited from THM-4139, not a new generalization of that theorem.

At parameter -7/4 the *algebraic* three-cycle instead has polynomial
P=x^3+x^2/2-9x/4-1/8, irreducible over Q, and

    f^3(x)-x=(f(x)-x) P(x)^2.

Its multiplier is1. There is no conflict with R1: those three points are
not rational. A six-step return also occurs on every three-cycle, but its
**exact** period remains three.

## 3. What a rational sixth dynatomic point would mean

Define the polynomial

    Phi_6(x,c) = (f^6(x)-x)(f(x)-x)
                    / ((f^3(x)-x)(f^2(x)-x)).

This quotient is interpreted by polynomial cancellation, not by evaluating
an undefined 0/0. Its degree in x is54.

**Lemma R2 (PROVED).** For rational x,c, Phi_6(x,c)=0 if and only if x
has exact period6 under f_c.

Any such zero satisfies f^6(x)=x, so its exact period divides6. At an
exact three-cycle the residual factor is proportional to lambda_3+1;
at an exact two-cycle it is 1+lambda_2+lambda_2^2; at a fixed point it is
lambda_1^2-lambda_1+1. The latter two quadratics have discriminant-3
and no rational root. R1 excludes lambda_3=-1. These exhaust the proper
divisors of6. For exact period6 neither denominator factor vanishes, so
the converse is immediate. The multiplier factors follow by expanding
f^d(x+h)=f^d(x)+lambda_d h+O(h^2), or polynomial division. ∎

Thus the six-cycle problem asks for rational points on a particular
arithmetic curve, not for a six-order matrix. In
[Stoll, *Rational 6-cycles under iteration of quadratic polynomials*](https://arxiv.org/html/0803.2836v2),
nonexistence is conditional on the asserted analytic properties of the
L-series and weak BSD for the Jacobian of the genus-four quotient curve
X_0^dyn(6). The paper obtains a rank-three subgroup and uses the rank bound
to complete its rational-point argument. That conditional theorem is not
made unconditional by the elementary fixed-parameter checks above.

This citation concerns the verified paper; it is not a claim that a
comprehensive literature search has excluded every later advance.

## 4. An exact bridge: two six-order lifts produce the Q--Pell matrix

Let

    B = [[0,-1],[1,1]],   U_u = [[1,u],[0,1]],
    C_u = U_u^(-1) B^(-1) U_u.

Both B and C_u have cube -I and sixth power I. With
K=[[2,1],[1,2]], one has B^T K B=K; C_u preserves the positive form
U_u^T K U_u. Nevertheless,

    B C_u = [[1,u],[u,u^2+1]],
    det(B C_u)=1,    tr(B C_u)=u^2+2.                     (1)

For real u!=0 the product is hyperbolic, with an eigenvalue greater than1.
At u=2, with V=[[1,-1],[0,1]],

    V^(-1) (B C_2) V = [[3,4],[2,3]].                    (2)

This is exactly the Pell matrix in
[golden/nested trees, section 6](golden_nested_trees_20261013c.md), the
conditional two-generation Q covariance continuation of
[THM-4630, Hofstadter Q and the silver-ratio ARMA fixed point](../../01-canon/theorems/THM-4630-hofstadter-q-jitter-is-the-silver-ratio-arma-fixed-point-rho1-equals-1-minus-sqrt2.md).
Equation (2) is an algebraic identification of matrices. It does not
identify the polynomial orbit with Q or establish Q's sampling hypotheses.

**Consequence (PROVED).** Separate bounded cyclic generators need not have
bounded changing-generator compositions. There is no common positive
definite quadratic form making both B and C_u nonexpansive when u!=0:
such a form would make their expanding product nonexpansive. Each form's
change of frame must be paid. The golden-tree comparison and Q projective
normalization controls are therefore complementary, not interchangeable.

For Collatz, this does not refute the terminal norm theorem: its actual
signed identification maps already share an entrywise-l1 bound. It rules
out replacing those maps by an attractive two-coordinate recurrence while
silently dropping the switching cost or physical scale.

## 5. Quadratic invariance can conceal an impossible pointwise state

For the signed three-cycle A=-P_3,

    A^3=-I, A^6=I, A^T I A=I, det(I-A)=2.

The nonzero covariance Q=I is stationary under congruence, but the
pointwise relation z=Az forces z=0. Stationarity of a distribution and
a pointwise native relation are different predicates.

The correct second-moment condition for a pointwise linear relation
Lz=0 is

    L Q L^T=0, Q positive semidefinite,                  (3)

not merely covariance invariance. Factoring Q=RR^T shows that (3) implies
LR=0 and hence LQ=0. A signed-cycle compiler retains the relation itself;
it does not manufacture independent terminal signs from stationary
moments. The marked ROOT coordinate and its first-moment cross terms
must also be retained.

This is the precise map from the cycle analogy back to the
[Gram construction](collatz_gram_slack_20261015b.md): a quadratic projection
can forget phase; authenticated native edges restore constraints. The
six-cycle and Pell examples are hostile controls on a proposed closure,
not extra Collatz orbit edges.

## 6. Reproduction

    python3 -B 04-computation/experiments/quadratic_cycle_phase_20261015b.py

The script uses exact SymPy and rational arithmetic. It verifies both
six-order generators, their positive forms and Pell conjugacy, the
covariance hostile, all12 rational preperiodic candidates and their full
edges, the two degree54 dynatomic polynomials, and the parabolic cubic
factorization. The all-parameter R1/R2 proofs are symbolic arguments above;
the finite checks are controls, not their quantifiers.
