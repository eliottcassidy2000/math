# Independent audit of nonrational scalar times in the invariant Hamiltonian family

**Status: ACCEPTED ANALYTICALLY + independently FINITE-EXACT.**
The audited [report](planar_jc_long_20260906_nonrational.md) proves that,
over every characteristic-zero field `K`, for every integer `a>=2`, every
nonconstant polynomial `f in K[z]`, and every nonzero scalar `lambda`,
the completed Hamiltonian flow of `S=f(p^a(p³-y²))` cannot have **both**
coordinate images `E_lambda(p),E_lambda(y)` in `K(p,y)`.
No repair to its final mathematical statement was required. The Laurent
field, finite-automorphism, and source-completion comparison justifications
discussed during review are explicit in the audited report.

## 1. The carrier, ambient field, and actual iteration

Writing `I=p^a(p³-y²)` and `f(z)-f(0)=z g(z)` gives
`S=f(0)+p²D[p^(a-2)g(I)]`, where `D=p³-y²`. The exponent restriction
`a>=2` is therefore sufficient for membership in the inherited universal
source carrier. Its completed source and depth properties are audited
separately in
[the all-R analytic audit](planar_jc_long_20260906_hamiltonian_analytic_audit.md).

For rationality one needs the complete **field** `K(s)((tau))`, containing
`K(p,y)` by `p=s²+tau`, `y=sp`. The coefficient ring `K[s]((tau))` alone
is not a field. The report correctly separates this valued-field
construction from the entire `D`-adic completion of `B_2`; it does not
use an unproved map between those two completions.

If `f_k z^k` is the first nonconstant term, then `I=tau p^(a+2)` gives
`ord_tau(delta s)>=k`, `ord_tau(delta tau)>=k+1`. Differentiation in `s`
preserves the valuation of a Laurent series with coefficients in `K(s)`,
and differentiation in `tau` lowers it by at most one. Thus

```text
ord_tau(delta F)>=ord_tau(F)+k
```

for every nonzero Laurent series, unless its derivative is zero. This
includes negative starting valuations and rational functions of `s`.
Consequently the exponential converges on the whole field. At every
prescribed coefficient cutoff, only finitely many operator iterates can
contribute; Leibniz and binomial identities are therefore legitimate
coefficientwise. They give multiplicativity, the group law
`E_lambda E_mu=E_(lambda+mu)`, and inverse `E_(-lambda)`. In particular
the later use of positive iterates is justified and not merely formal
notation for an unproved scalar substitution.

The exact velocity has the especially simple form

```text
delta p=2s I f'(I).
```

Its leading coefficient is
`2k f_k s^((2a+4)k+1)` in degree `tau^k`; subsequent iterates of `p`
have order at least `2k`. Therefore `E_lambda(p)-p` has that same leading
coefficient multiplied by `lambda`. Substituting `n lambda` is valid by
the group law, and it stays nonzero for every positive integer `n` in
characteristic zero. This proves infinite order on the actual field
coordinate `p`, rather than only non-local-nilpotence of the derivation.

## 2. Generic curve and the contradiction from rationality

The invariant `c=I` is transcendental over `K`. In fact `p,c` are
algebraically independent and `K(p,y)` is quadratic over `K(c)(p)`.
For even `a=2r`, the substitution `v=p^r y` gives
`v²=p^(a+3)-c`; for odd `a=2r+1`, the substitution `v=p^(r+1)y` gives
`v²=p^(a+4)-cp`. Both substitutions are birational on the function
fields, even though powers of `p` appear in their denominators inversely.

The even radicand has no common root with its derivative because `c!=0`.
For the odd radicand, zero is a simple root; at a nonzero root its
derivative is `(a+3)c!=0`. These facts remain true over an algebraic
closure of `K(c)`. Squarefreeness and odd degree also make each radicand
a nonsquare in that rational function field, so the curve is geometrically
integral. Its smooth proper model has the claimed constant field and
geometric genus. Counting the simple finite branch points and the one
at infinity gives respectively `r+1` and `r+2`, hence genus at least two.
This is an application of the degree-two Riemann--Hurwitz formula, not
an inference from the singular special fibre `c=0`.
[Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B).

Suppose both images were rational. The completed automorphism then
restricts to an injective `K(c)`-homomorphism of `K(p,y)` into itself,
because it fixes `I` and sends nonzero denominators to nonzero elements.
The corresponding nonconstant rational map extends to a morphism of the
smooth proper generic curve. The correspondence does not require that
the inverse field map already be rational.
[Stacks, Theorem 53.2.6](https://stacks.math.columbia.edu/tag/0BY1).

Characteristic zero makes this selfmap separable. If its degree were
`d`, Riemann--Hurwitz would give
`2g-2=d(2g-2)+deg(Ram)`. Since `g>=2` and the ramification divisor is
effective, `d=1`, so the map is an automorphism.
[Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B).

For clarity, the finite-automorphism step includes its genus hypothesis:
the tangent line bundle has degree `2-2g<0`, hence has no nonzero global
sections. The automorphism group scheme is finite type and has zero
tangent space, so it has finitely many geometric points; in characteristic
zero it is also reduced. The `K(c)`-automorphism group injects into this
finite group. Some positive iterate would consequently fix `p`, contrary
to the explicit leading coefficient above.
[Stacks, Section 109.7](https://stacks.math.columbia.edu/tag/0DST).

This argument uses the invariant `I` itself. Replacing the constant field
by `K(f(I))` without checking geometric components could obscure the
argument; the report correctly avoids that step. Critical special fibres
where `f'(c)=0` do not invalidate a statement about the generic invariant
curve or the original two-variable rational function field.

## 3. The necessary bridge to actual source specializations

Equality of the two rational function fields alone does not identify the
two completions. The additional comparison supplied during the audit is
valid and supplies the polynomial-source consequence. Both logarithmic
images of `p,y` belong to `K[s][[tau]]`, hence to `K[[s,tau]]`. Define

```text
sigma:K[[s,tau]] -> K[x][[t]], s -> xt, tau -> t.
```

This map is well-defined: the coefficient of `t^n` is a finite sum over
`i+j=n`. It is injective because `s^i tau^j` maps to the distinct monomial
`x^i t^(i+j)`. It is also continuous from the `(s,tau)`-adic topology to
the `t`-adic topology. The logarithmic iterates of `p,y` are polynomials
in `s,tau`, and their `tau` orders tend to infinity. Their series therefore
converge in the bivariate completion as well. The exact bracket change of
coordinates shows term by term that their images under `sigma` are the
convergent source-coordinate exponentials of `p,y`. Passing to the limits
is justified by this continuity, rather than by an asserted embedding
of the entire `D`-adic completion in a Laurent field.

Suppose a source image were a rational function in `x,t`. Substitute
`x=s/tau`, `t=tau` in that rational function and clear powers and other
denominators to write it `g/f` with `f,g in K[s,tau]`, `f!=0`.
If `P` is the corresponding logarithmic image, then
`sigma(fP-g)=0` in `K[x][[t]]`. The element `fP-g` belongs to
`K[[s,tau]]`, so injectivity gives `P=g/f`. A rationality claim for both
source images would consequently contradict the logarithmic theorem.
This uses only injectivity of `sigma` on the displayed bivariate ring.

Finally `D=t³(1+x²t)²` makes the map from `B_2` to its source ring
continuous from the `D`-adic topology to the `t`-adic topology:
`D^N B_2` maps into `t^(3N)K[x,t]`. It extends to a continuous comparison
from the full `D`-adic completion to `K[x][[t]]`. A polynomial
specialization of the former completed operation would give polynomial
source images under this comparison, which the preceding argument
excludes. No injectivity of this larger comparison map is asserted or
needed. This resolves the topology issue for the stated consequence.

## 4. Scope and the decisive hostile

For constant `f` or `lambda=0`, the flow is the identity, so these
exclusions are necessary. The theorem does not assert that each separate
image is nonrational, that either is transcendental over `K(p,y)`, or
that no other completed carrier can give a rational specialization.
By the explicit comparison above, it does exclude a polynomial source map
inducing this particular completed operation. It says nothing about all Keller
pairs or compositions of different operations.

The hostile `S=u=x²t` is valid. Its derivation has
`delta^n x=n!x^(n+1)`, but its flow is

```text
x -> x/(1-lambda x), t -> t(1-lambda x)².
```

This is a rational symplectic map, with inverse time `-lambda`, preserving
the invariant `u`. Thus failure of local nilpotence cannot establish the
new conclusion. This hostile is outside the stated invariant family and
does not assert convergence in its `D`-adic topology. Its generic
invariant curve is rational, so the genus argument has exactly the
intended missing hypothesis.

## 5. Independent finite controls and frozen manifest

Standalone audit
[source](../../04-computation/planar_jc_long_20260906_nonrational_audit.py)
and [output](planar_jc_long_20260906_nonrational_audit.out):

```bash
python3 -B 04-computation/planar_jc_long_20260906_nonrational_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_nonrational_audit.py
```

The source imports neither the producer nor any repository mathematical
module. One engine works in symbolic rational function fields to verify
the two generic models and their squarefreeness for every `2<=a<=32`.
It also checks the singular `c=0` and genus-one `a=0` boundaries.
A separate sparse polynomial-ring engine constructs the completed
derivation and its truncated exponential for 45 cases: `2<=a<=10`, five
polynomials `f` with first nonconstant orders one, two, and three, and
three nonzero scalar times per case.

It verifies the leading displacement, all generator velocity valuation
bounds, Laurent monomial bounds including negative valuations, actual
component substitution into the invariant, inverse and additive-group
composition, and six actual positive iterates in four cases. The rational
hostile is checked by its differential equations, invariant, group law,
source Jacobian, and ten nonzero iterated derivatives. These are exact
controls of the universal analytic proof; they do not enumerate rational
selfmaps or replace its curve-theoretic step.

The added completion-comparison controls use a second polynomial ring
and the literal `(x,t)` bracket. They test the monomial injectivity and
bracket Jacobian of `sigma`, then compare three nontrivial flows and
their iterates through source row thirteen, at two nonzero times. In
each comparison the source displacement is explicitly nonzero; the
control is not confined to a prefix before the first response.

Normal and optimized output agree byte for byte in **2,268 live gates**.

```text
audit source SHA256:
86e64c180b3c42c40b892e144a38713c0a6e5268183ba361e72214c0c2226578
audit output SHA256:
1bc00cad425e823961d2690ab6f35b486921203832396bb03bc4d3fab6706f0b
```
