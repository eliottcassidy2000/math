# Independent audit of the all-order moving collision period

**Status: ACCEPTED ANALYTICALLY, with independent FINITE-EXACT controls.**
The accepted conclusion is the fixed-source-volume obstruction for every
member of the specified three-parameter compensated Russell family and
every nonzero clock in `t k[[t]]`, over every characteristic-zero field.
The formal source-coordinate escape is also accepted. No polynomial
termination, polynomial primitive pair, or planar JC consequence follows.

Audited report:
[planar_jc_long_20260906_collision.md](planar_jc_long_20260906_collision.md).
Producer source:
[planar_jc_long_20260906_collision.py](../../04-computation/planar_jc_long_20260906_collision.py).
Its frozen source/output SHA-256 values are respectively
`7320aecbc79c7c53c0118889a4c307f53102f6aa488a35641faac16873bf1e82`
and
`069811ba7e42bc8b83d0a0bb76c589b0f673406f12ecfde8fa4cd83b96d62d6c`.
No mathematical repair to the report was required by this audit.

## 1. All-order proof and the actual target module

Differentiating the common-image equation gives
`f_t+z_i'f_x=gamma'`. At the same target section, a two-form therefore
evaluates as `J_i=Omega(f_x,gamma')`; the tangential correction vanishes by
alternation. Every formal relation `sum lambda_i f_x=0` yields
`sum lambda_i J_i=0` identically. This proof uses all form coefficients
at their common target value. It requires neither closedness nor a
two-dimensional target, and asserts necessity only.

For the Russell surface times the `w` line, `(C,E,w)` are local coordinates.
The implicit surface derivative in `B` is a unit at the collision. Hence
a relation between the `(C,E)` tangents also annihilates their `B`
components, because those components are obtained with the same two
coefficients at the common target point. The `w=t` coordinate has
`w_x=0`. Thus the relation is valid in the full target; no extra two-form
slot is lost by using the surface plane.

The compiler determinant is the exact polynomial identity
`det d(C,E)/d(x,q)=6(B+2)`. Pullback by
`q=Q+s+chi(s)x` gives density `6(B+2)(1+chi'(s)x)`.
The common `B` value has constant term zero, so its factor is a unit and
may be cancelled from the period identity. This proves

```text
Pi(s)+chi'(s) Xi(s)=0,
Xi(0)=-5/18+13/18=4/9.
```

The normalized relation denominator has constant term `54`, so the
relation itself is well-defined. Since `Xi` is a unit and the coefficient
field has characteristic zero, nonzero `chi` has
`ord Pi=ord chi-1` and leading coefficient `-4n chi_n/9`. The implication
`chi'=0 => chi=0` also uses the inherited condition `chi(0)=0`.

## 2. Exhaustiveness over all parameter triples and all clocks

I checked the literal four coefficient formulas against
[THM-4424 / russell-constant-normal-debt-discriminant-contact-correspondence](../../01-canon/theorems/THM-4424-russell-constant-normal-debt-discriminant-contact-correspondence.md).
They are successive-stratum identities, not four numerical samples. The
first zero equation is exactly the affine plane parametrized by `(p,r)`;
the second solves uniquely for `p` as a polynomial in `r`. On that
parabola the next coefficient is a nonzero rational multiple of `F(r)`.
If it vanishes, the following coefficient is a nonzero rational multiple
of `G(r)`. The exact Bezout identity for `F,G` makes simultaneous vanishing
impossible in any characteristic-zero extension field. Irreducibility of
`F` is useful for the quotient-field control, but is not needed for this
nonvanishing implication. Thus the four first orders `2,3,4,5` exhaust
every triple `(a,b,c)` in the stated family.

For any nonzero `phi(t)` in `t k[[t]]`, composition of a nonzero series
of order `n-1` has exact order `ord(phi)*(n-1)`. Its leading coefficient
is `(-4n chi_n/9)*lead(phi)^(n-1)` and is nonzero. The tangent relation
is composed with the same clock and still has last coordinate zero.
Therefore the hypothetical constant density has a nonzero weighted
period, contradicting the necessary identity. No factor `phi'` belongs
in this constant-density period. In contrast, the *specific surface
form* used to construct a compatible density does acquire `phi'` by
the chain rule. The report keeps these two claims separate and fixes
`w=t` throughout. The zero clock and clocks with nonzero constant term
are outside its quantifiers.

The argument excludes every formal target two-form in the completed
common target germ. It consequently excludes regular target forms with
that constant pullback as well. It does not identify later retained
target-lift debts with periods, and does not prove sufficiency of a zero
period in a three-dimensional target.

## 3. The source-volume boundary is sharp in the stated formal category

For `phi=t`, the normalized regular form
`Omega_0=(dC wedge dE)/(6(B+2))` pulls back with density
`rho=1+chi'(t)x`, a formal unit. It has zero weighted period by the exact
identity above. The change

```text
X=x+chi'(t)x²/2, T=t
```

has determinant `rho` and is the identity modulo `t`. Its inverse exists
in `k[X][[T]]`. One explicit construction is the Catalan series inverse
of `X=x+h x²/2`, with `h=chi'(T)` in `T k[[T]]`; each coefficient is a
polynomial in `X`. This works simultaneously in the labelled formal
neighborhoods. In the new coordinates the same target form has constant
density. This does not contradict the fixed-coordinate theorem.

More generally tangents transform by `T_i -> T_i/a_i`, where
`a_i=psi_x(z_i,t)`, so the transported relation is `lambda_i -> lambda_i a_i`.
The sum of relation coefficients becomes `Pi+chi' Xi=0` for the displayed
change. The actual weighted density identity is unchanged. Thus omitting
the source Jacobian would give a false coordinate-invariance statement;
the report explicitly retains it. Only `phi=t` is used for this unit
escape. For a ramified clock, `phi'` need not be a unit, and this argument
does not silently treat it as one. Polynomiality of the change, its
inverse, the compensator, or a primitive target pair is not claimed.

## 4. Independent exact engine

Audit source:
[planar_jc_long_20260906_collision_audit.py](../../04-computation/planar_jc_long_20260906_collision_audit.py).
Output:
[planar_jc_long_20260906_collision_audit.out](planar_jc_long_20260906_collision_audit.out).
This engine imports neither the producer nor any repository computation.

The producer uses hand-written Fraction arithmetic for its exceptional
quartic quotient and solves one coefficient at a time. The audit instead
uses sparse SymPy polynomial rings and its native algebraic field domain.
It derives the entire six-equation Jacobian from the actual compiler with
symbolic `(a,b,c)`, then solves by full truncated-series contraction
`state <- state-J(0)^(-1) residual`. Each iteration gains one formal
order; it does not read the producer's coefficients or Jacobian literal.
The three rational strata and the actual quartic field are solved through
degree five. The normalized relation and all three target two-form basis
slots are independently evaluated through degree four, respecting the
one-order loss under differentiation.

The audit additionally checks the surface equation, compiler determinant,
all-root Bezout certificate, each actual first-live period, all five
retained clocks by composition of the complete computed period, the
compatible density at all three sections, the exact source Jacobian and
relation covariance, and the Catalan inverse through degree six in its
auxiliary coefficient. These **133 live gates** are controls of the
analytic all-order theorem, not finite evidence substituted for its proof.

```bash
python3 04-computation/planar_jc_long_20260906_collision_audit.py
python3 -O 04-computation/planar_jc_long_20260906_collision_audit.py
```

Normal and optimized runs agree byte for byte. The audit source SHA-256 is
`c548ccc11c05903e1cfb315c68fb0c93fc0aa2183dae36288c0890ece60c9c36`;
the output SHA-256 is
`1edd32e9922f196c411a6f2d6ceb308f9fff747feeb587863305f7d87feea4f5`.
