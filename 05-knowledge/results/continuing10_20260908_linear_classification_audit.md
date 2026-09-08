# Independent audit of the complete global source-linear classification

**Verdict: PASS, with no mathematical repair.** The frozen primary
`continuing10_20260908_dg_linear_all_m.md` and its source/output have been
compared against the independent argument below. Work over C, m>=1, and the
actual W_m charts inherited from `planar_jc48_sep08_dg_genus.md`.
Source-linear means degree at most one in t on the original A2_(x,t).
The function F is assumed nonconstant.

## The full coefficient space and rational integrability

Write F=A(x)t+B(x), A=sum a_i x^i. Substitution in the second chart gives

    F=-r^m A(1/r)-r^(2m)A(1/r)b+B(1/r).

Since B cannot cancel a b-dependent pole, global regularity forces deg A<=2m.
The remaining negative powers force the complete condition

    B=B0+sum_(i=m+1)^(2m) a_i x^(i-m).

Conversely it gives the polynomial second-chart formula

    F_inf=B0-sum_(i=0)^m a_i r^(m-i)
         -b sum_(i=0)^(2m) a_i r^(2m-i).

No basis element or constant term is omitted. In particular A=0 would make
F constant, and nonconstant F therefore has A!=0. At D={r=0},

    F|D=B0-a_m-a_(2m)b,
    F_r|D=-a_(m-1)-a_(2m-1)b, F_b|D=-a_(2m).

The positive-degree-part convention makes B0 the actual constant coefficient
of B. It differs by a scalar shift from a full-polynomial-part convention.

Because C(x,t)=C(F)(x), a rational mate for F exists exactly when the rational
differential dx/A(x) has a rational primitive over C(F). Here is a complete
degree argument. Suppose A is nonconstant, has degree n and r distinct roots.
A simple root is impossible: a derivative of a rational function has zero
residue, whereas1/A would have a nonzero simple-pole coefficient there.
At a root of multiplicity n_i>=2, a rational primitive has pole order n_i-1.
It has no other finite pole, so its total finite pole degree is n-r.
At infinity the primitive has a finite limit and differs from that limit by
a nonzero term of order x^(-(n-1)). Thus its zero order at infinity is n-1.
Equality of total zero and pole degrees forces n-1<=n-r, hence r<=1.

This reasoning works over C(F) directly; no specialization or assumed
descent of a primitive is needed. Consequently the complete possibilities are

    A=a in C*, or A=a(x-h)^n with a!=0 and n>=2.

Both converses are explicit: -x/a or1/[a(n-1)(x-h)^(n-1)] has D_F derivative1.
The excluded n=1 case would require a logarithm. This is an all-degree proof,
not an inference from a finite rational-antiderivative search.

## Exactly which rationally integrable functions are globally submersive

For constant A=a, globality makes B=B0. The original chart is smooth, and
F_inf=B0-a r^m-a r^(2m)b has a nonzero boundary derivative precisely when m=1.

For A=a(x-h)^n, 2<=n<=2m, the only possible affine critical line is x=h.
If n<=m, B is constant and the entire line is critical. If n>m, direct
alternating-binomial evaluation gives

    B'(h)=a(-1)^(n-m-1)binom(n-2,m-1)h^(n-m-1).

This is a nonzero monomial when h!=0. The general weighted identity follows
from the same partial-binomial sum used in the all-m maximal-power theorem;
there is no further exceptional complex value of h.

The boundary formulas above now finish the classification. If n=2m,
F_b|D=-a!=0. If n=2m-1, F_b|D=0 and F_r|D is a nonconstant affine function
of b, so there is exactly one boundary critical point. This statement does
not deny additional affine critical points in a degenerate h=0 case.
If n<=2m-2, the boundary condition is the single constant a_(m-1)!=0.
Combining it with the affine condition yields exactly:

| Surface | Complete globally submersive, rationally integrable list |
|---|---|
| m=1 | A constant nonzero, or A=a(x-h)^2 with arbitrary h |
| m>=2 | a*h!=0 and either m+1<=n<=2m-2 or n=2m |

All listed functions use the unique global B above and arbitrary B0. Empty
integer intervals are interpreted literally. In particular W2 admits only
the quartic pure-power branch; its cubic middle branch has a boundary critical
point. At h=0,n=m+1 the original source can be smooth, but the boundary still
prevents a global submersion when m>=2.

## Unit response orders and fibre scope

For every admitted pure-power function let c0=B(h), g=F-c0 and e=n-1. Its
special affine fibre has the exact factorization g=(x-h)K, with
K(h)=B'(h)!=0. The residual K is irreducible, linear in t and primitive;
its component is disjoint from x=h. The two-component valuation proof in
the independently audited all-m theorem therefore applies without change:

    Ann_(C[F])([1])=((F-c0)^(n-1))

in the ORIGINAL quotient C[x,t]/D_F(C[x,t]). Explicitly,
g^e/[a e (x-h)^e]=K^e/(a e) is a polynomial primitive of g^e. A rational
fibre-constant correction cannot cancel a pole on x=h without introducing one
on the other component, so no smaller order can kill the unit class.
The constant-A branch has a polynomial mate and unit class zero.

The primary's canonical-derivative claim is also accepted. For every admitted
pure power, the source gradient is a unit ideal and the rational constant field
is exactly C(F), so the proved `planar_jc48_sep06_torsion.md` applies in its
stated ring. The only reducible source fibre is c0; its two components give
one principal-part arm. The unit's top coefficient is
(B'(h))^(n-1)/[a(n-1)], nonzero modulo the common component diagonal. The
canonical connection differentiates its principal part, giving exact
c0-primary order n-1+j after every j>=0. This does not transfer a connection
to an unexamined global coordinate ring.

Every admitted pure-power source fibre is generically G_m on U0: its x-line
omits h. If n=2m, the nonconstant boundary restriction adds exactly one point,
giving A1 on W_m, as in the accepted alternative-chart theorem. If n<2m, D
maps to one constant value, so a generic fibre misses it and remains G_m.
For constant A on W1, the generic fibre is already A1. No alternative
A2_(F,1/(x-h)) is asserted for the lower-power G_m branches.

The resulting positive orders of the distinguished unit class are

    m=1: {1};
    m>=2: {m,m+1,...,2m-3} union {2m-1}.

Across all surfaces these are exactly{1} union{3,4,5,...}: for any e>=3 take
m=e,n=e+1. **Only the unit class** lacks order2. The response module can contain
order2, for example the canonical derivative of an order1 unit class. The
classification must not be paraphrased as absence of all order-two torsion.

This connects directly to the earlier boundary-linear order-two supplier:
v=x+x^3t is the m=2,n=3,h=0 middle branch, whose missing hypothesis is global
submersivity at its boundary critical point. The higher-t thickenings in that
supplier are outside this source-linear classification.

## No global mate for any nonconstant global source-linear function

This additional corollary does not require F to be globally submersive.
A hypothetical G in O(W_m) with D_F G=1 restricts to a polynomial mate on U0,
and in particular to a rational mate. Thus A must be constant or a pure power.
For a pure power, B'(h)=0 would make the gradient vanish on x=h. Otherwise
the same two reduced special components preclude a polynomial mate. Hence
only A=a constant could remain.

In that case F=a t+B0 and every rational mate is -x/a+H(F). Polynomiality
on U0 forces H polynomial, since (x,F) are polynomial coordinates there.
On Uinf, -x/a=-1/(ar) has a pole and H(F_inf) is regular. It cannot cancel
the pole. Thus no nonconstant global source-linear F has a globally regular
unit-Jacobian mate. This statement retains the original source two-form and
does not claim a new general plane Jacobian theorem.

## The m=1 constant branch really changes when the ring changes

The optional ring-change corollary is accepted by a direct argument. On W1,
omega=dr wedge db is nowhere zero. For F=a t+B0 the Hamiltonian derivation
is -a partial_x on U0 and

    D_F=a r^2 partial_r-a(1+2rb) partial_b

on Uinf. Both expressions are polynomial and compatible, so D_F actually
preserves O(W1). The global response quotient is therefore defined here.

In the original C[x,t] quotient, the unit is zero via q=-x/a. Globally, let
g=F-B0. The function gq=-xt is regular on both charts, with second expression
1+rb, and D_F(-xt)=g. But q has a genuine simple pole on D={r=0}. The other
global special-fibre component T is the closure of t=0 and is disjoint from D:
its second-chart factor is1+rb, comaximal with r. A rational correction H(F)
which cancels the pole of q on D introduces one on T, where q is regular.
The already proved constant-field equality excludes any other correction.
Thus the global unit is nonzero and has exact annihilator(g).

This shows source theta=0 and global theta of order1 for the same function F.
The restriction map on response quotients kills that nonzero global class.
The proof does not apply the affine-plane component theorem to O(W1). For
m>=2 the two-form can vanish on D and the derivation need not preserve the
global ring; no corresponding global module is inferred there.

## Hostiles and independent exact controls

The rational-integrability assumption is essential even for globally smooth
functions. On W2,

    F=(x^2-1)^2 t+x^2

is globally submersive: the affine derivatives at its two possible critical
lines are2 and-2, and its boundary b derivative is-1. But1/(x^2-1)^2 has
residues-1/4 at1 and1/4 at-1, so no rational mate exists. This is a control
of the recovered residue mechanism, not a priority claim.

The independent source reconstructs the global coefficient space by polynomial
division for m=1,...,8, checks its entire boundary jet, then checks every pure
power2<=n<=2m in both h=0 and h!=0 branches. It independently matches the full
classification and each surface's positive unit-order spectrum, checks the
rational primitives and special-fibre witness identities, and retains the
degree2m+1 globality hostile, the middle-exponent boundary critical point,
the constant-A polynomial positive control and the smooth two-root residue
obstruction. These controls validate the formulas; the analytic arguments
above prove the complete unbounded quantifiers.

The source also independently verifies the actual m=1 global Hamiltonian,
both expressions of the global annihilator witness, its source bracket,
the two comaximal global components and the nonzero boundary principal part.

After filing, run:

    python 04-computation/continuing10_20260908_linear_classification_audit.py
    python -O 04-computation/continuing10_20260908_linear_classification_audit.py

Both modes pass583 always-active exact gates, produce identical raw LF stdout,
and regenerate identical certificate bytes. The source imports SymPy but no
other mathematical producer. It writes its certificate to results when filed.
Final frozen pins are:

    Audit source d394b6c047e09f3805c71898d82dda31ef655cfb038942f08d8da34fd025f275
    Audit output 608d8bbecb1619a822dcfce11065775bccadafbf4f30ff2771573101982e8b8d
    Audit certificate 76ef3545c535f60e0c7a50aee90368c86b10307360e7c3c8b9add820b16e0ed7
    Primary source 48ac8b527ba5e88dbb38ab1ca287824e8748574e4927c37bfeb8b2051784de3d
    Primary output 9a4eda039ccb20db4d68322818fb21ea331a38c3c1111816ca94ca34862afb1b

No primary source or frozen output was changed. Root owns promotion,
maintained routing, integration and Git.
