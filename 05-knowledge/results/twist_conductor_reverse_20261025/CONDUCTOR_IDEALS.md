# Conductor ideals without unique factorization

2026-10-25. **PROVED** by the arguments below. Internal proof review;
finite monomial controls are separate. No new classification of normal
orders, divisor classes, or torsors is asserted.

## I1. The exact coefficient-ideal formula

Let R be a normal Noetherian domain with fraction field K. Let V be a
finite abelian group, and let L be a finite field extension with a marked
K-basis u_z, z in V, satisfying u_0=1 and

    u_z u_y = c(z,y) u_(z+y),       c(z,y) in K^*.

Suppose the finite R-algebras

    B = direct_sum_z J_z u_z  subset  A = direct_sum_z I_z u_z

have nonzero finitely generated fractional coefficient ideals J_z subset I_z,
with I_0=J_0=R. Assume A is normal. Both span L over K, so have the same
fraction field L. A is finite over B and is its normalization. Their conductor is

    (B:A) = direct_sum_z C_z u_z,
    C_z = intersection_(y in V) (J_(z+y) : c(z,y) I_y).       (I1)

Here (J:M)={f in K : f M subset J}; the y=0 condition already places
C_z inside J_z. To prove (I1), multiply an arbitrary sum by every
homogeneous I_y u_y. For each fixed y, different z have distinct target
characters, so cancellation cannot hide any component. These conditions
are necessary and sufficient. This proof uses ideals, not prime elements.

If every J_z is reflexive (divisorial), (I1) has the equivalent valuation form

    C_z = {f in K : v_P(f) >= M_P(z) for every height-one P},
    M_P(z) = max_y [v_P(J_(z+y))-v_P(I_y)-v_P(c(z,y))].       (I2)

The valuation of a fractional ideal is the minimum valuation of its
generators. A divisorial ideal is the intersection of its height-one
localizations. Apply this to each target J_(z+y) and to every generator of
I_y; finite intersections give (I2). The maxima are taken separately at
each prime. They cannot be replaced by a sum of maxima when several
radicands share that prime.

This is the same-field order problem. The warning in
[THM-3397, sections 4 and 7](../../../01-canon/theorems/THM-3397-torsor-killing-versus-effective-boundary-valuations.md)
still applies: the conductor of a genuine nontrivial field-extension
quotient is zero, and killing a divisor class does not erase poles.

## I2. Non-UFD alone does not break the old one-step theorem

In the previous binary construction A_0=R[u_i]/(u_i^2-b_i), assume the
nonzero b_i have pairwise disjoint height-one divisor support. R may be
nonfactorial. Use **two complementary character packets V=C direct_sum D**
and B=A_C A_D, exactly as in the preceding one-step theorem. This paragraph
does not assert that test for three or more packets. The coefficient ideals
of B are still the principal ideals
product_i b_i^e_i(z) R. At each height-one prime, at most one b_i has
positive valuation, so the original B1 and D2 valuation arguments apply,
with that valuation as a multiplicity. Intersecting over height one gives

    (B:A_0) = (product_(shared i) u_i) A_0,
    (B:A_0) = {x in B : u_i x in B for every i}.            (I3)

No normalization claim is made unless A_0 is normal. The hypotheses are
disjoint divisor support and the specified coefficient lattices, not
unique factorization. To obtain genuinely nonprincipal conductors one
must allow nonprincipal coefficient ideals or more general orders.

## I3. A nonprincipal conductor and a codimension-two hostile

Let k=Q, R=k[x,y,z]/(xy-z^2), and A=R[w]/(w^2-2).
The ring R is the invariant ring k[s,t]^{+/-1}, with
x=s^2,y=t^2,z=st. It is normal: an invariant rational function integral
over R is integral over k[s,t], hence polynomial and invariant. The same
argument over k(sqrt(2)) proves A normal.

For any nonzero ideal I of R, put B_I=R direct_sum Iw. It is an algebra,
since (Iw)^2=2I^2 subset R. It is birational to A and finite, with

    (B_I:A)=I A.                                          (I4)

Indeed r+sw and its product with w lie in B_I exactly when both r and s
belong to I. Thus the one-step test is exact here, with the actual generator w.

Take I=(x,z). This is a height-one prime and hence divisorial. Its extension
to A is not principal: at the vertex (x,y,z), x and z are independent
minimal generators modulo (x,y,z)I, over k(sqrt(2)). A principal ideal has
at most one such generator. More explicitly,

    I^2=x(x,y,z),       (I^2)^(**)=(x),

where ** is reflexive hull. Its class is the nontrivial two-torsion class,
consistent with [THM-3397's A1 class-group computation](../../../01-canon/theorems/THM-3397-torsor-killing-versus-effective-boundary-valuations.md).
The conductor is a genuine ideal, not a chosen scalar valuation.

For the hostile take I=m=(x,y,z). Every height-one localization has
m_P=R_P, because the maximal ideal m has height two. Thus all height-one
comparisons say B_m and A agree. Nevertheless (B_m:A)=mA is proper:
1 belongs to B_m but w does not. A valuation-only reconstruction would
give the reflexive hull A and miss this obstruction. This is exactly why
the reflexivity hypothesis in (I2) is required. Normality of A alone is
not enough. Formula (I1) remains valid without that hypothesis.

## I4. Normalization can enlarge the necessary operator bank

Let R=Q[pi]_(pi), w^2=2, t^2=pi, and A=R[w,t], a normal DVR with
uniformizer t and residue field Q(sqrt(2)). Set

    u_1=t,  u_2=wt,
    B=R[u_1,u_2]=R direct_sum Rt direct_sum Rwt direct_sum Rpi w.

Then A is the normalization of B and

    (B:A)=tA.                                           (I5)

Multiplication by w forces the constant coefficient into pi R; the
displayed basis then verifies both containments in (I5). The element 1
passes multiplication by each old port root u_1 and u_2, yet fails
multiplication by the integral ratio w=u_2/u_1. So the old-root one-step
test does not survive this normalization unchanged.

This example is still a product of complementary normal quotient orders:
let H flip w and K flip both w,t. Then A^H=R[t], A^K=R[wt], and their
tensor product has image B. The two original radicands pi and 2pi share
a divisor, which the previous disjoint-support theorem forbids.
The port-preserving repair is to retain u_1,u_2 and add the integral ratio
w, together with its relation u_2=wu_1. Erasing that ratio erases residue
field information. Checking {1,w,t,wt}, an R-basis of A, detects the full
conductor in one finite bank. This last basis test is general linear
algebra; the small single-root theorem is the special stronger result.

## I5. Relation to ancestry

The [valuation/LCA packet](../session56_rmq_lca_connections_20261025/valuation-lca/REPORT.md)
recovers equality depth from a low-digit trie. Such a tree records
valuations of differences, not unit squareclasses or coefficient-ideal
membership. I4 supplies a precise additional coordinate: both radicands
have the same positive valuation while their ratio is the nonsquare 2.
I3 supplies another: identical height-one data can hide a codimension-two
ideal. LCA remains an exact index for its declared predicate; it is not
an exact rational-lift or conductor observer without these extra data.
