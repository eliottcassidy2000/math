# The boundary-linear family realizes an exact order-two unit response

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; controls FINITE-EXACT.**
This supplies an actual smooth polynomial and its full component data to
the session's intrinsic torsion connection. It is a further consequence
of the proved boundary-linear exclusion, not a new case of JC(2).
The h=0 case is already the r=3 one-arm control in
[the earlier torsion note §5](planar_jc48_sep06_torsion.md), up to output
scaling and sign. The present supplier is the full h in b^3 C[b]
thickening with the same actual labelled arm.

## 1. The precise bridge

Use b=-x^2(1+x^2t), v=x(1+x^2t), rb=-v on the specified DG surface.
For every a in C*, c in C and h in b^3 C[b], put

    F=a rb+c+h(b)=-a v+c+h(b),    g=F-c,
    D_F=F_x partial_t-F_t partial_x,
    C_F=C[x,t]/D_F(C[x,t]),       theta=[1].

The claim is

    g^2 theta=0,             g theta!=0.                (1)

Thus the unit response has exact c-primary order two, for the whole
unbounded family including h=0. Let nabla be the canonical transverse
connection from [the proved component-jet theorem](planar_jc48_sep06_torsion.md).
Then

    ord_c(nabla^j theta)=j+2,            j>=0.           (2)

In particular every canonical divergence class nabla(theta) has exact
order three and is nonzero. The sequence is linearly independent over C.

The closest proved mechanisms are [the entire boundary-linear
carrier](planar_jc48_sep08_boundary_linear.md) and the complete torsion
isomorphism in the cited component-jet theorem. The former pays the
smooth polynomial and rational constants field that were previously
missing for this wildcard lane. The canonical hostile is its quartic
F=b^4-v: a rational mate exists, yet its principal parts are unequal
on two components of F=0. The corrected near miss is to identify a
generically exact unit response with a zero polynomial response class.

The map sends theta, via its actual rational primitive, to the full tuple
of scalar principal parts at the components of the special fibre. It
preserves the original polynomial ring, F, its derivation, and component
labels. Passing modulo the common diagonal loses only the allowable
rational repair H(F), exactly as required by the constant-field theorem.
The finite computation below checks a polynomial witness for g^2 theta=0;
the nonzero order and its unbounded derivative growth are analytic.

## 2. Both inherited hypotheses now hold

The boundary-linear proof shows that F has no critical point on A2_(x,t).
On Uinf its derivatives are a b and a r+h'(b), so their only common
zero is r=b=0, on the added boundary D. On E={x=0}, F_x=-a. These
charts cover the affine source after adjoining E, so (F_x,F_t)=(1)
in C[x,t]. No smoothness on all of W is asserted.

Moreover C(x,t)=C(F,b), with r=(F-c-h)/ab. On this rational function
field,

    D_F=(ab/r^2) partial_b |_F,

whose kernel is exactly C(F). Thus the two stated hypotheses of the
torsion theorem are both satisfied, without an unproved chart-entry
or geometric-integrality assertion.

## 3. The full principal part of the unit

Set L(b)=integral h/b^3 db and M(b)=integral h^2/b^3 db, with zero
integration constants. The proved rational mate is

    G0=a^(-3)[-g^2/(2b^2)-2g L(b)+M(b)],
    D_F G0=1,
    G0=-r^2/(2a)+Q(b,rb),              Q polynomial.   (3)

It has no affine pole except along E, where its order is exactly two.
The special fibre F=c has E and Gamma={1+x^2t=0} as distinct reduced
components. It may have further components; all are retained in the
following argument. Because F is smooth, they are disjoint in the
affine source. The function G0 is regular on every component other than
E, since its only affine pole is x=0.

At E,

    g=-a x-a x^3t+O(x^6),
    g+a x is divisible by x^3.

Consequently a^2/g^2-x^(-2) is regular at its generic point. Formula
(3) therefore has the exact scalar principal part

    pp_E(G0)=-a/(2g^2),
    pp_i(G0)=0                  for every i!=E.        (4)

There is no g^(-1) term. This conclusion uses the whole h in b^3 C[b]
condition; the absence of a quadratic x term in g is explicit.

Under the canonical module isomorphism

    tors(C_F)=direct_sum_d
       (C^(components of F=d)/C diagonal)
          tensor (F-d)^(-1) C[(F-d)^(-1)],

theta maps to

    -(a/2) [e_E] tensor g^(-2).                        (5)

Here [e_E] is nonzero because Gamma is a different component with zero
principal part. Formula (5) proves both directions of (1). An actual
polynomial witness is P2=g^2 G0: its only former affine pole is cleared,
so it is polynomial by normality, and D_F P2=g^2. No polynomial
primitive for g can exist, since its component tuple is the nonzero
multiple -(a/2)[e_E] tensor g^(-1).

## 4. The canonical derivative retains the obstruction at every order

For a polynomial Bezout field V with V(F)=1, the inherited operator is

    nabla[q]=[V(q)+(div V)q].

It is independent of V and differentiates the component principal parts.
Applying it to (5) gives the explicit formula

    nabla^j theta maps to
       -(a/2)(-1)^j(j+1)! [e_E] tensor g^(-j-2).       (6)

All coefficients are nonzero in characteristic zero. This proves (2)
and linear independence by the distinct highest pole orders. The
smallest C[F]-submodule stable under nabla and containing theta is the
single full principal-part arm [e_E] tensor g^(-1)C[g^(-1)]: multiplying
theta by g reaches its first level, and differentiation reaches every
higher level. It is not a finite response packet.

This is an actual connection with the previous wildcard, not merely a
similarity between two pole arguments. It still supplies no map from
the moving-source Hamiltonian carrier or the collision quadrics to this
fixed smooth F. Those other entry problems remain distinct and OPEN.

## 5. Reproduction and scope

The matching source checks the exact primitive and polynomial g^2
witness in both charts, the scalar coefficient in (4), the absence of
a simple principal part, and the intrinsic arm differentiation rules.
Its finite h-bank is explicitly bounded; equations (3)–(6) prove the
full unbounded statement.

    python3 -B 04-computation/planar_jc48_sep08_unit_torsion.py
    python3 -B -O 04-computation/planar_jc48_sep08_unit_torsion.py

The principal-part theorem credits the classical relative-exactness
antecedents; no priority claim is made here. The [independent analytic
and source audit](planar_jc48_sep08_unit_torsion_audit.md) passes. Both
normal and optimized runs agree with the frozen output on all 86 gates.

Source SHA256: `f6f9f99d62216f78e366a30a1452a392ee046966ed914aa1cee88c55268d1ba2`.

Output SHA256: `6aec5ecea6bebbf93d59dbc70b5bcc117db3b41cf7e91320343368eafcec4d2b`.
