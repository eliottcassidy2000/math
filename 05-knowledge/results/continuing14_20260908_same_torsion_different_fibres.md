# The generic fibre recovers a polynomial that the full pointed torsion misses

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.**
For the entire polynomial mutation family from continuing13, the generic
original affine fibre over `C(T)` determines its primitive polynomial `Q`
up to an affine change of the fibre coordinate. Two explicit degree-fifteen
members nevertheless have isomorphic full ambient torsion Weyl modules
with their distinguished units, identical target supports and component
multiplicities, the same intrinsic rational-pair degree, and the same
ordinary puncture count, while their generic affine fibres are not
isomorphic over `C(T)`. No claim is made after arbitrary algebraic base
extension, about an entire response module beyond its torsion, or about
a polynomial constant-Jacobian mate.

## 1. Inheritance and the recovered coordinate

The direct supplier is the [global polynomial mutation theorem](continuing13_20260908_quintic_mutation.md)
and its [independent audit](continuing13_20260908_quintic_mutation_audit.md).
The [fixed-support iteration theorem](continuing13_20260908_fixed_support_iteration.md)
already shows that the pointed unit can hide unbounded component
multiplicities even when all its supports and source parameters stay fixed.
The present test restores those multiplicities and asks what remains lost.

A forgotten adjacent mechanism is THM-2800,
[two-pole bitangent eliminant and complete Nielsen corridor](../../01-canon/theorems/THM-2800-two-pole-two-double-zero-stieltjes-recurrence-and-first-nielsen-pair.md):
there a response passport loses the choice of a bitangent, already giving
two distinct normalized maps in the first nontrivial case. That is a
different response layer and is not imported as a global-surface theorem.
Here the mutation supplier realizes a similar loss on the actual two charts
of W2, and a separate puncture argument distinguishes the generic fibres.

The external terminology is classical: polynomials with two finite critical
values are related to bicolored plane trees; vertex-degree data form a
coarser combinatorial record. See Betrema and Zvonkin,
[Plane trees and Shabat polynomials](https://doi.org/10.1016/0012-365X(95)00127-I),
1996, and Dupont et al.,
[Pairs of tree dessins, their Shabat polynomials, and monodromy groups](https://arxiv.org/abs/2510.10192),
2025. These references motivate the missing-configuration test. No
classification, existence result, monodromy computation or priority claim
from them is needed for the explicit proof below.

The concept board is full pointed torsion, target multiplicities, polynomial
composition, critical-point configuration, degrees of closed punctures,
and actual source-boundary derivatives. The map sends `Q` with squarefree
square root of its derivative to a global mutation. It preserves the
critical-value multiplicities and loses the affine class of `Q` under the
torsion quotient. The sidecar restored below is the generic punctured curve
as a curve OVER `C(T)`, not only its completed genus or number of punctures.

## 2. A generic-fibre reconstruction theorem for the complete supplier family

Let `Q_i'=f_i^2`, with each `f_i` squarefree of degree `q_i>=2` and
`f_i(0)=Q_i(0)=0`. Choose any admissible source parameters `h_i,lambda_i`
as in the mutation theorem, and put `K=C(T)`, with the SAME abstract target
coordinate in both members. Their generic ORIGINAL affine fibres are

    U_i = Spec K[z, f_i(z)^(-1), (T-Q_i(z))^(-1)].                 (1)

Then

    U_1 is K-isomorphic to U_2
    iff Q_2(a z+b)=Q_1(z) for constants a!=0,b in C.             (2)

To prove this, complete each curve to `P1_K`. Its missing closed points are
`q_i+1` distinct constant degree-one points (the roots of `f_i` and infinity)
and ONE further closed point of degree `2q_i+1`, defined by `Q_i(z)-T`.
The latter polynomial is irreducible in `K[z]`: `Q_i(z)-T` is irreducible
in `C[z,T]` since its quotient ring is `C[z]`, and Gauss's lemma applies.
It is separable in characteristic zero and avoids every constant puncture.

A `K`-isomorphism of the affine curves extends uniquely to their smooth
projective completions. It preserves the residue degrees of the missing
closed points. Thus it sends the constant punctures to the constant
punctures and the single higher-degree puncture to its counterpart.
There are at least three constant punctures; prescribing the images of
three of them determines a Mobius transformation with CONSTANT coefficients.
Consequently the extended map is some `phi in PGL_2(C)`.

At the generic point of the higher-degree puncture, `T=Q_1(z)`, and its
image satisfies `T=Q_2(phi(z))`. Hence `Q_2 o phi=Q_1` as rational
functions in `C(z)`. A polynomial has its sole pole at infinity, so
`phi` fixes infinity and is affine. This proves necessity in (2).

Conversely, differentiating the affine identity gives
`a*f_2(a z+b)^2=f_1(z)^2`. It therefore identifies the constant punctures
as well as the moving puncture, and restricts to the required isomorphism.

This is an iff for the displayed GENERIC curves. The parameters `h,lambda`
do not occur in (1), although they matter to the total fibration and the
extra source component. The argument uses degrees of closed points over
`K`; it must not be silently replaced by a statement over an algebraic
closure of `K`, where those degrees all become one.

## 3. Two explicit primitive polynomials with the same complete passport

Let

    S(w)=6w^5-15w^4+10w^3,
    6r^2-15r+10=0,
    d=-1-r, w_r(z)=d z^3+r,
    Q_r(z)=S(w_r(z)).                                           (3)

Use the two distinct complex roots `r_+,r_-` of the quadratic. Its
discriminant is `-15`; neither root is `0,1,-1`. Choose either square root
`s_r` satisfying `s_r^2=90d`, and put

    f_r(z)=s_r z w_r(z)(w_r(z)-1).                              (4)

Direct differentiation gives

    S'(w)=30w^2(w-1)^2,
    Q_r'=f_r^2, deg Q_r=15, deg f_r=7.                         (5)

The three factors in (4) have disjoint simple zeros: `z=0` is disjoint
from the cubic preimages of `0,1`, since `r!=0,1`; each of those cubics
has three simple roots since `d!=0` and its value at zero is nonzero.
Thus `f_r` is squarefree, its zero at zero is simple, and `Q_r(0)=S(r)=0`.

The roots of `f_r` split into FOUR critical points of value zero
(`z=0` and the three roots of `w_r=0`) and THREE of value one
(the three roots of `w_r=1`). All have local degree three. In fact

    S(w)=w^3(6w^2-15w+10),
    S(w)-1=(w-1)^3(6w^2+3w+1).                                (6)

After removing the four cubic zeros of `Q_r`, three simple roots remain;
after removing the three cubic zeros of `Q_r-1`, six simple roots remain.
The primitive polynomials thus have the same full ramification partitions
over `0,1,infinity`, namely `(3^4,1^3)`, `(3^3,1^6)`, and `(15)`.
These partitions concern `Q_r`, not critical points of the submersion below.

Fix the SAME source parameters

    h=-1/2, lambda=1.                                          (7)

At the corresponding source point `z=-2h=1`,

    w_r(1)=-1, Q_r(1)=-31, f_r(1)=2s_r.                       (8)

The value `2s_r` is nonzero and never `-1`: the latter would imply
`360d=1`, whereas `d=-1-r` is nonreal. Both square-root choices therefore
satisfy every global mutation hypothesis, with no varying source parameter.

## 4. The actual global first functions and all equal invariants

In the original source ring set

    X=x+1/2, z=X+1+X^3 t,
    A_r=X f_r(z), T_r=A_r+Q_r(z), G_r=1/(2A_r^2).              (9)

The independently audited mutation theorem proves that each `T_r` is a
global submersion on BOTH charts of W2 and `J_(x,t)(T_r,G_r)=1`.
Its original `t`-degree is fifteen with leading coefficient `6d^5 X^45`.
The inherited actual-source derivative on `X=0` is
`2s_r(1+2s_r)!=0`; the derivative on the added divisor is `-f_r'(0)!=0`.

The full original-source torsion has FOUR arms at target `0`, THREE at
target `1`, and ONE at target `-31`. At each of these special values
there is also the supplier's one regular component. The total special
fibre component counts are therefore respectively `5,4,2`.

Use the same abstract target coordinate `tau` for both modules. The
complete primitive principal parts are

    root components: 1/[2(tau-c)^2], c=0 or1,
    X=0 component: (1+2s_r)^2/[2(tau+31)^2],
    regular component at each support: 0.                    (10)

There are no simple-pole terms. The full torsion Weyl module is the
direct sum of four standard towers at `0`, three at `1`, and one at
`-31`. Its distinguished unit generates one direction at each support,
with exact scalar annihilator

    [tau(tau-1)(tau+31)]^2.                                  (11)

The first two coefficient vectors in (10) are already identical for the
two examples; rescaling the final one-dimensional block by a nonzero
constant matches the third. These block maps commute with multiplication
by `tau` and with the canonical derivative. Thus the FULL ambient
TORSION modules, together with their distinguished unit and its embedded
generated submodule, are isomorphic. This is stronger than only matching
the unit-generated modules or the scalar annihilator.

For both examples the intrinsic rational-pair degree is thirty and every
rational constant-Jacobian mate gives the same embedded pair field for
its own first function. Every ordinary original source fibre has twenty-three
punctures; the projective completion is rational. For each example these
punctured fibres vary non-isotrivially as the target value varies, as paid
by the supplier. No equality of the two embedded pair fields is asserted.

## 5. Unequal affine classes, hence unequal generic source fibres

Suppose `Q_-(a z+b)=Q_+(z)` with `a!=0`. Each polynomial is centered:
its degree-fourteen coefficient is zero. Comparing degree fourteen forces
`b=0`, because the next nonzero possible degree below fifteen is twelve.
After writing `w=d_+ z^3+r_+`, the remaining identity becomes

    S(L(w))=S(w),
    L(w)=(d_- a^3/d_+)(w-r_+)+r_-.                           (12)

An affine self-map preserving `S` must permute its two critical points
`0,1`. Their critical values are distinct, so each point is fixed; hence
`L` is the identity. Equation (12) then forces `r_-=r_+`, a contradiction.
By (2), the two generic affine source fibres are not `C(T)`-isomorphic.

A cheap exact distinguishing coordinate is also available. For a centered
degree-fifteen polynomial with coefficients `a_15,a_12`, the ratio

    I=a_12^5/a_15^4
     =15^5(2r-1)^5/6^4                                      (13)

is invariant under nonzero scaling of `z`. Its two values are distinct
in the exact quadratic field; the certificate records their nonzero
difference. The proof using critical points in (12) explains the
distinction without relying on this computation.

The common support set `{0,1,-31}` has no nontrivial affine self-symmetry.
In particular a fixed-target comparison cannot quietly exchange `0,1`:
the tempting outer reflection satisfies `S(1-w)=1-S(w)`, not `S(w)`.
Nor may one pull back through an outer critical point instead of the
chosen regular root `r`: taking `r=0` makes the derivative square root
have a repeated zero and destroys the submersion hypothesis.

## 6. Verification and precise survivor

The [standalone producer](../../04-computation/continuing14_20260908_same_torsion_different_fibres.py),
[raw transcript](continuing14_20260908_same_torsion_different_fibres.out),
and [certificate](continuing14_20260908_same_torsion_different_fibres_certificate.json)
check both parameter conjugates in `Q[r]/(6r^2-15r+10)`, with no numerical
roots or supplier-engine import. There are 95 always-active exact gates.
The universe comprises the two formulas (3), their complete critical
groups and polynomial fibre multiplicities, three declared ordinary
values `-3,2,5`, both admissible square-root choices through their squares,
the scaling control, and the collision/target-swap hostiles. The generic
reconstruction iff and nonisomorphism have analytic proofs above; finite
specializations are not their proof.

Reproduce in the repository:

    python 04-computation/continuing14_20260908_same_torsion_different_fibres.py
    python -O 04-computation/continuing14_20260908_same_torsion_different_fibres.py

An outside source writes its certificate beside itself; a filed source
writes under `05-knowledge/results`. Frozen normal and optimized executions
must reproduce the same raw LF stdout and certificate.

The refuted inference is that restoring full component multiplicities to
the pointed unit determines the generic original source curve. The first
failed step identifies a ramification passport with an affine class of
polynomials. The strongest survivor is the reconstruction iff (2): the
generic curve does retain the missing class when the base field and
degrees of its closed punctures are retained. Whether the same pair
remains nonisomorphic after arbitrary algebraic base extension is left
OPEN here. General planar JC and LRC(14) remain OPEN.

Independent acceptance: [referee and reproducible controls](continuing14_20260908_passport_independent_audit.md). Only this audited status and routing paragraph were added after the referee-pinned producer report; the proof and computational bytes are unchanged.
