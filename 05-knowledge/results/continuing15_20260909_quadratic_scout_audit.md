# Independent referee: quadratic pole corridor and the exact eigenvalue lattice

**Status: INDEPENDENT AUDIT ACCEPTED, analytically and by exact controls.**
This accepts the frozen [quadratic logarithmic scout](continuing15_20260909_quadratic_logarithmic_scout.md),
report SHA256 `e969113c9048dcc9ea30195aadaa1993e445235a4f00ed57014f988c510e04cc`.
No producer engine is imported or executed by this referee. Its claimed
general statements are established by the arguments below, not by the
finite sample. The nonsquarefree global W2 problem remains open.

## 1. Normalization and raw squarefree existence

The pole orders were recomputed from the actual local maps. At a fixed
root of odd discriminant order m, x-alpha=u^2 gives
dx/y~2u^(1-m)du. In particular m=3 is a double pole. An earlier informal
root suggestion of order 2e-1 for m=2e+1 was wrong; the frozen report
correctly uses 2e. Even order 2e gives order -e. At a common root,
generic v prevents cancellation and ord D=min(ord A,2 ord B). All
orders >=3 fail, and order two has constant integer residues only when
A has order >=3 and B is simple, with reciprocal derivative integral.
The Jacobian of the Hamiltonian vector field at its fixed point then
forces B'=+/-1 if F is a source submersion.

At infinity the complete normalized orders are d/2-2 for even degree d
and d-3 for odd d. Degrees zero and one are impossible. In the
squarefree case degrees >=3 make dx/y holomorphic on the proper curve,
also impossible. Degree two leaves integer residues and exactly
deg A<=1, E quadratic with leading coefficient 1/n^2. These statements
include geometric normalization and do not silently discard a square factor.

The raw converse has an independent two-variable proof. On y^2=D,
use delta=y partial_x+(D_x/2) partial_y. For
D=kappa^2 x^2+d1 x+d0, the element
U=y+kappa x+d1/(2kappa) satisfies delta U=kappa U.
Replacing d1 by e1+4a1 H makes U a polynomial in the original x,t.
For kappa=1/n, U^n is an eigenfunction of eigenvalue one. The norm U V
is a nonzero element of C(H), so U parametrizes the generic conic;
this pays geometric integrality and the rational constant field.

## 2. The unbounded spectrum is exact

For H_m=x^(m+2)t^2+xt the normalized generic curve has
Y^2=1+4v x^m and R=(Y-1)/(Y+1). The identities
dlog R=m dx/(xY) and x^m=R/[v(1-R)^2] were independently recomputed.
The degree-m cyclic cover is irreducible by the simple valuation at R=0.
For a rational eigenfunction of eigenvalue lambda, residues first force
lambda to be an integer. The deck quotient sigma(F)/F is a constant
m-th root of unity, so F=x^j f(R). Its two valuation congruences are
lambda=j modulo m and 2j=0 modulo m. Equivalently m divides 2lambda.
Thus the spectrum lies in [m/gcd(m,2)]Z even over algebraic constants.

The positive and negative generators in the producer are literal
polynomials over C. Their brackets and products are checked separately:
their product is H_m^(m+2) for odd m and H_m^((m+2)/2) for even m.
Their powers establish both inclusions in the rational AND polynomial
spectrum, including zero. Eigenvalue one remains only at m=1,2.
The omitted K0=0 component of H_m=0 forbids any pole in the scalar
coefficient K(H_m), so these eigenfunctions retain a repeated x factor.
The stated whole-family source-critical consequence is valid.

The m=3 genus-one explanation also stands independently: a logarithm
with the stated two residues would have divisor P_plus-P_minus and
give a degree-one map to P1. Its triple is principal, and the divisor
itself is not. This is a genuine failure of local-to-global integration.

## 3. Minimum witness degree and the W2 sidecar

For F=U P(H), H=-UV/4, U=x+2t, V=x-2t, the displayed inverse formulas
give C(x,t)=C(F)(H). Thus every rational primitive differs by C(F).
The zero components have distinct H labels, with nonzero labels on all
hyperbolas and label zero on U=0. A scalar of F cannot simultaneously
cancel their differing simple principal parts. The original unit is
nonzero and killed by F.

Every polynomial repair differs by C(F) intersect C[x,t]=C[F]. This
intersection uses Bezout for coprime numerator and denominator polynomials
in one variable; it does not assume a general Jacobian centralizer theorem.
Since deg_t F>=3, no nonconstant polynomial in F can cancel the degree-two
part of H. Minimum witness degree is indeed two. Direct chart substitution
gives the claimed pole of degree 2deg(P)+1, so none is a global W2 example.
The admissible-double-root gauge example is source-submersive but also has
a boundary pole. Both hostiles are necessary to delimit the conclusion.

## 4. Reproduction

The separate engine uses normalized local substitutions, the conic
derivation directly in (x,y), polynomial generators and products for
m=1..10, a valuation-congruence lattice check through m=25, shifted
nonlinear witness families, and the actual W2 Laurent chart. It uses
raising gates under normal and optimized Python and records exact output.
Run `python 04-computation/continuing15_20260909_quadratic_scout_audit.py`
or the same command with `python -O` after relocation.

No mathematical repair of the frozen producer is required. Promotion
from its pending-audit status to independently audited status, with a
link to this referee, is authorized; no other producer bytes need change.
