# Six roots, seven canonical forms, and actual quotient maps

**Status:** PROVED over C, with the inherited surface theorem explicitly stated
below. The comparison with the classical diagonal-quadrics K3 model is CITED;
the character, local-resolution, fixed-point, and deletion calculations are
proved here. No historical novelty claim is made.

## 0. Inheritance and scope

Inherited input: Session59 `toric-shared-roots/TORIFICATION_OBSTRUCTION.md`,
Theorem T1. For distinct spanning linear forms L_0,...,L_(m-1) on C^3, with
no four branch lines concurrent, the inverse image X_m of the parameter plane
under coordinate squaring is a normal complete intersection of m-3 quadrics.
It maps to P^2 with degree 2^(m-1), has only A1 singularities, and its minimal
resolution Y_m is crepant. In particular Y_6 is K3; Y_7 has p_g=7, q=0,
K^2=16. The hourglass forms are

    E, E+x, E-x, E+z, E-z, E+x+z, E-x-z.

Closest mechanism: the square-root monic basis, not an analogy of integers.
Hostile example: a single simultaneous pair flip is NOT the subgroup generated
by two separate flips. Corrected near miss: an equivariant isomorphism of
underlying canonical line bundles is not the naive coordinate linearization.
Least-used sidecar: the actual triple-incidence sets, needed when lifting
projections to resolutions and counting fixed points.

All assertions retain characteristic zero and the no-fourfold condition.
Statements about Hodge structures are over C and do not assert integral
lattice isometries, rational-point transfers, or arithmetic descent to Q.

## 1. The deck group and its eigensheaves, proved from local bases

Write roots r_0,...,r_(m-1). The effective projective deck group is

    G_m = F_2^m / <(1,...,1)>.

The common minus sign is the identity in projective space. Its character group
is therefore the even-weight code

    E_m = {I subset {0,...,m-1}: |I| even},
    chi_I([epsilon]) = product_(i in I) epsilon_i.

Let p:X_m -> P^2 be the root cover. For |I|=2k, the root product
z_I=product_(i in I) r_i is a section of
O_X(2k)=p^*O_P2(k). Multiplication by it defines a character-eigenmap
O_P2(-k) -> p_*O_X. These maps give the isomorphism

    p_*O_X = direct_sum_(I even) O_P2(-|I|/2)_[chi_I].       (1)

**Proof of surjectivity as well as injectivity.** On the base chart L_j != 0,
normalize r_j=1. The algebra has free basis all squarefree monomials in
u_i=r_i/r_j, i!=j: each u_i satisfies the monic relation
u_i^2=L_i/L_j, and the full degree is 2^(m-1). For an even I the normalized
z_I is exactly the basis monomial with support I\{j}. Every subset A of the
other m-1 indices occurs once: take I=A if |A| is even and I=A union {j}
otherwise. Thus the displayed maps are an isomorphism on each chart. These
charts cover P^2 because the linear forms span. This proves (1), including
at branch crossings, without assuming an abelian-cover formula.

Finite flat duality gives p_*omega_X=Hom(p_*O_X,omega_P2). Every character is
self-inverse, so (1), omega_P2=O(-3), and crepancy give

    H^0(Y_m,K_Ym)_[chi_I] = H^0(P^2,O(|I|/2-3)).             (2)

Here is also a direct differential check, so the regularity at a triple
point is not hidden in duality. In homogeneous base coordinates [a:b:c], put
Omega=a db wedge dc - b da wedge dc + c da wedge db. For a polynomial F of
degree |I|/2-3 the form

    F(a,b,c) Omega / product_(i in I) r_i                  (3)

is homogeneous of degree zero on the root surface and has character chi_I.
At a single or double branch crossing, ds=2u du cancels each possible root
denominator. At a triple point the local equation is
w^2=alpha u^2+beta v^2, alpha beta != 0. If I contains local roots with
indicators a,b,c in {0,1}, (3) is a unit times

    u^(1-a) v^(1-b) w^(1-c) (du wedge dv / w).

The last parenthesis generates the dualizing sheaf of the A1 hypersurface,
and crepancy extends it to the resolution. Conversely, a rational
chi_I-eigenform is a base rational function times (3). A pole of order d
in its coefficient on a branch line would have order 1-a-2d < 0 upstairs
for every d>=1; an ordinary base divisor also forbids a pole. Thus its
coefficient has no divisorial poles and is precisely a global section of
O(|I|/2-3). This independently verifies completeness of (2).

Consequently

    p_g(Y_m) = sum_(2k<=m) binom(m,2k) h^0(P^2,O(k-3)).     (4)

For m=5,6,7,8,9 this gives 0,1,7,31,111, agreeing independently with the
complete-intersection Hilbert-series calculation.

## 2. The canonical sign twist is essential

Let Q_1,...,Q_(m-3) be the defining diagonal quadrics. The residue expression
for a canonical section with numerator P of degree m-6 transforms under a
deck sign change epsilon as

    P(r) Res(Omega_(P^(m-1))/(Q_1...Q_(m-3)))
       |-> (product_i epsilon_i) P(epsilon r) Res(...).      (5)

The quadrics themselves are fixed. Thus the adjunction isomorphism
omega_X ~= O_X(m-6) has a determinant twist in its natural linearization.
Under epsilon=(-1,...,-1), the two factors contribute
(-1)^m(-1)^(m-6)=1, exactly as projective descent requires.

For m=6, the one-dimensional H^{2,0} has the total-parity character.
For m=7, the coordinate numerator r_j gives the character

    c_j = chi_({0,...,6}\{j}),
    c_j(epsilon) = product_(i!=j) epsilon_i.               (6)

All seven characters occur once. They span E_7, have rank six over F_2,
and have precisely one linear dependence: c_0+...+c_6=0. In particular,
these seven canonical characters are NOT the seven nonzero vectors of
F_2^3 and are not by themselves a Fano plane.

The representation is faithful. Equivalently choose the unique even-weight
representative of every element of G_7 (m is odd); the seven canonical
characters are then its seven coordinate signs. Hence the canonical image
of G_7 is exactly the determinant-one diagonal sign subgroup of GL_7.

If a pure root-coordinate permutation is induced by a linear transformation
A of the base forms, the canonical representation is the coordinate
permutation representation twisted by det(A). Indeed the exact sequence
relations -> C^m -> C^3 shows that det(permutation)/det(relations)=det(A).
This is the extra sign required for a reflection of the hourglass arrangement.

## 3. Six roots give fifteen Nikulin involutions, with eight fixed points each

For m=6 set

    H = ker(chi_all:G_6 -> {+1,-1}) ~= (C2)^4.              (7)

It is exactly the symplectic subgroup inside the deck group, not a claim
about the full automorphism group of the K3 surface. Its fifteen nonidentity
elements are represented uniquely by two-coordinate flips {i,j}, since
weight-four supports are their complements modulo the global sign.

Every such involution fixes exactly eight points on the smooth K3 Y_6.
Here is a direct count valid for our nodal specializations.

* At an ordinary intersection L_i=L_j=0 there are eight points upstairs;
  all are smooth and fixed by the {i,j} flip.
* If L_i,L_j meet on a third branch line L_k, there are four A1 points
  upstairs, not eight smooth points. At each one the involution acts by
  (u,v,w)->(-u,-v,w). On the exceptional conic
  w^2=alpha u^2+beta v^2 in P^2 it fixes precisely the two points w=0.
  Therefore its resolution has 4*2=8 fixed points.
* A fixed projective point must be supported entirely in one sign
  eigenspace. The other possibility would require all four coordinates
  outside {i,j} to vanish, forbidden by the no-fourfold assumption.

The blow-up of each A1 point is its minimal resolution. Every deck
transformation preserves the singular set and its ideal, so the action
lifts to this blow-up and the above count is genuinely on Y_6.

No point is fixed by two different elements of H: the stabilizer acts
faithfully on the tangent space by finite-order linearization, and a
nontrivial order-two matrix in SL_2(C) is -Id. An elementary abelian
2-subgroup of SL_2(C) has order at most two. Hence Y_6/H has exactly
15*8/8=15 A1 points. This does not say that its double-plane contraction
has fifteen A1 points when the branch lines have triple crossings.

Topological Lefschetz gives trace(g|H^2)=8-2=6 for each g!=1. Averaging,

    dim H^2(Y_6,Q)^H = (22+15*6)/16 = 7.                  (8)

The rank-fifteen orthogonal complement is rational of type (1,1), hence
algebraic by Lefschetz (1,1). If the six-line arrangement has t triple
points, the sums of the four exceptional curves over each triple are
independent H-invariant algebraic classes, orthogonal to the hyperplane
class. Thus

    rho(Y_6) >= 16+t.                                    (9)

This is a lower bound, not an equality assertion for special arrangements.

## 4. The genuine six-line double-plane quotient

The invariant monomials for the even-sign lift of H have parity either
all zero or all one. Therefore X_6/H is explicitly

    D_6: t^2 = product_(i=0)^5 L_i
    in weighted P(1,1,1,3),    t=product_i r_i.             (10)

The quotient has degree sixteen. At an ordinary branch intersection D_6
has an A1 point; at a triple crossing its equation is t^2=uv(u+v), a D4
rational double point. Its minimal resolution is K3. For t triple points
it has (15-3t) A1 and t D4 singularities. In contrast, Y_6/H has fifteen
A1 points: over each triple it retains a rational curve carrying three
of them; contracting that curve gives the D4 model. Confusing these two
quotients loses the resolution sidecar.

**Primary-source recovery.** Garbagnati--Sarti, *Kummer surfaces and K3
surfaces with (Z/2Z)^4 symplectic action*, Section 10.2, printed pp.54--55,
explicitly gives the diagonal three-quadric model, even-sign action,
six-line double-plane quotient, and a genus-two Kummer subfamily.
https://sarti.pages.math.cnrs.fr/home/GarbagnatiSarti_RMJ.pdf
https://arxiv.org/abs/1305.3514
Our nodal fixed-point verification above is supplied explicitly rather than
silently replacing the smooth model by a singular one.

The symplectic involutions here are Nikulin involutions in the standard
sense, and (C2)^4 is the corresponding classical abelian symplectic group.
We do not claim it is one of Mukai's maximal groups or identify the whole
automorphism group. For historical primary records see Nikulin (1979),
https://www.mathnet.ru/eng/mmo366, and Mukai (1988),
https://doi.org/10.1007/BF01394352. Their classifications are not needed
for the explicit fixed-point proof or the formula (8).

## 5. Seven roots: seven actual K3 deletion maps

Let sigma_j be the class of the j-th coordinate flip. Forgetting r_j gives

    X_7 -> X_6^(j) = X_7/<sigma_j>, degree 2.              (11)

There is no projection center on X_7: such a center would require all
six retained forms to vanish. The square of the omitted coordinate is a
linear combination of retained squares, and adjoining its root gives the
quadratic extension and the stated invariant algebra. The resulting
six-root arrangement still spans and has no fourfold crossing.

The map lifts to a proper generically finite morphism

    f_j:Y_7 -> S_j,    S_j = minimal resolution of X_6^(j).

To check this rather than assume functoriality of minimal resolutions:
a target node comes from a triple not containing j. L_j is nonzero there,
so the double cover is etale locally a disjoint pair of that node and its
blow-up lifts. A source node involving j instead maps to a smooth target
point; its exceptional curve is contracted by f_j. These are all cases.
Consequently f_j need not be finite, and Y_7/<sigma_j> need not itself be
the smooth K3 S_j. Its minimal smooth model is S_j.

The pullback of the unique K3 two-form spans the c_j line of H^{2,0}(Y_7).
Indeed c_j(sigma_j)=+1 and c_k(sigma_j)=-1 for k!=j. Thus

    H^{2,0}(Y_7) = direct_sum_(j=0)^6 f_j^*H^{2,0}(S_j).  (12)

For the hourglass, omitting E removes all three triple points; this S_0
has a smooth three-quadric model. Omitting any outer line leaves two
triples and eight A1 points in the six-root model.

## 6. The minimal hostile test: two different meanings of 'pair flip'

Let i!=j. The single involution sigma_i sigma_j fixes five canonical
characters and negates two. Its resolved quotient therefore has p_g=5,
NOT zero. The order-four subgroup <sigma_i,sigma_j>, on the other hand,
fixes none of the seven characters; its quotient is the five-root model,
whose smooth resolution has p_g=0. These are different maps and degrees.

More generally, in the unique even-weight representative for G_7, an
involution of weight w=2,4,6 has p_g of its resolved quotient 7-w, namely
5,3,1. There are respectively 21,35,7 such involutions. Thus NO nontrivial
involution quotient of Y_7 has p_g=0. The invariant-form assertion follows
either by the explicit field quotient or by descent of invariant forms
through finite quotient singularities and their resolutions.

A rank-three quotient G_7 -> F_2^3 with seven distinct nonzero inertia
labels is a different operation. Its pulled-back character space is the
simplex [7,3,4] code: all seven nonzero words have weight four, so (2)
shows that the kernel-invariant H^{2,0} is zero. This is the exact
character mechanism behind the Campedelli bridge in the separate
seven-line/Fano packet; it is not evidence that (12) was a Fano plane.

## 7. Connection ledger and failure boundaries

| Source -> target | Actual map | Preserved | Lost / sidecar | Cheapest test |
|---|---|---|---|---|
| labelled roots -> E_m characters | sign pairing modulo global sign | deck action, eigenspaces | line incidence not encoded | odd-weight word fails global-sign descent |
| six-root K3 -> six-line double plane | invariants r_i^2, product r_i | holomorphic two-form | 16-fold degree, A1/D4 contraction data | check triple local equation |
| Y_7 -> seven K3s | delete one root, resolve | seven distinct two-form lines | quadratic extension and contracted curves | c_k(sigma_j)=(-1)^(k!=j) |
| Y_7 -> five-root surface | quotient by two independent flips | function-field invariants | two signs, all two-forms | subgroup order is four, not two |
| E_7 -> Fano-labelled quotient | choose a rank-three quotient / simplex subcode | labelled inertia and quotient field | noncanonical kernel choice | nonzero subcode weights all four |

**Board after the calculation.** (1) Root-cover algebra now controls genus
by weight. (2) K3 symplectic actions give fixed-point and Picard constraints.
(3) Deletion controls individual canonical summands, but loses branch
extension data. (4) Fano structure requires an added quotient, not just
seven labels. (5) Resolution changes quotient models and cannot be dropped.
(6) Hodge correspondences strengthen the connection, as proved separately
in `HODGE_DECOMPOSITION.md`.

## 8. Reproduction

Run `python3 check_deck_characters.py`. The finite universe is all sign
classes and even character words for m=5,...,9, all seven one-flip and
21 two-flip tests, and every six-line deletion of the displayed hourglass.
It compares genus to an independent complete-intersection coefficient
calculation, verifies the canonical determinant twist, enumerates the
fixed-point count by actual branch incidences, and checks the conic/Gale
identities in `KUMMER_GALE_BRIDGE.md`. The script is a finite exact audit,
not a replacement for the sheaf and resolution proofs above.
