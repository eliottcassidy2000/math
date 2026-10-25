# Seven canonical forms, seven K3 quotients, and the different Fano seven

2026-10-25. **PROVED character and code identities with explicit proofs;
CITED classical duality and surface tools; independently audited.**
The general abelian-cover and Campedelli constructions are classical. The
contribution here is a precise transport of those structures to the inherited
root surfaces, with the canonical linearization and the two meanings of
"seven" separated. No surface-isomorphism, rational-point, F1-descent, or RH
claim follows.

The eigensheaf, determinant-twist and K3-transport calculations overlap concurrent [CHARACTER_SEXTICS, C1–C4](../six_seven_geometry_20261025/characters/CHARACTER_SEXTICS.md) and Session60. The amplitude/phase and abstract exterior-cube views are additional probes within these inspected packets. See [INTEGRATION.md](INTEGRATION.md).

## 1. Inherit the actual surfaces

Use Session59's audited
[root-surface theorem](../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md).
Over an algebraically closed characteristic-zero field, take n>=4 distinct
linear forms L_i on k^3 spanning its dual, with at most three branch lines
concurrent. The surface X is the coordinate-square inverse image of their
parameter plane in P^(n-1), and pi:X->P^2 has degree 2^(n-1). X has only
A1 nodes; its resolution Y is crepant, with K_Y pulled back from O_X(n-6).

Closest proved mechanism: that exact complete-intersection and node theorem.
Canonical hostile: a triple branch intersection changes the singularity of
a double-plane quotient even though the universal root cover has A1 nodes.
Corrected near miss: an underlying line-bundle isomorphism need not carry
the intended group linearization. Least-used sidecar: the determinant
character in adjunction.

Board: projective sign deck; even character weights; canonical residue;
six-line double planes; Fano/Hamming quotient; norm weights and phases.
Anchor: explain p_g=1 at six roots and p_g=7 at seven. Niche: exact quotient
maps rather than shared numerical invariants. Wildcard: the different
six-dimensional simplex and six relative phases of a seven-form packet.

## 2. The even subsets are the complete deck-character carrier

The projective sign deck is

    G=F_2^n/<(1,...,1)>.

Its character group is the even-weight subspace of F_2^n. An even subset S
defines chi_S(epsilon)=product_(i in S) epsilon_i. Simultaneous sign acts
trivially exactly because |S| is even.

**Theorem 1 (exact eigensheaves).** There is a G-equivariant decomposition

    pi_*O_X = direct_sum_(S even) O_P2(-|S|/2),             (1)

whose indicated summand has character chi_S.

**Proof.** On U_j={L_j!=0}, set r_i=x_i/x_j for i!=j. The finite free
coordinate algebra has the square-free monomial basis in the n-1 variables
r_i, with r_i^2=L_i/L_j. Even subsets S are in bijection with these monomials:
include j exactly when the chosen subset away from j has odd cardinality.
For |S|=2r, its local generator is

    e_(S,j)=product_(i in S) x_i / x_j^(2r).

Its character is chi_S. On U_j intersect U_k,

    e_(S,k)=(L_j/L_k)^r e_(S,j),

the transition law of O(-r). These bases cover P^2 and prove (1). QED.

The multiplication has not vanished: locally

    e_(S,j)e_(T,j)=e_(S triangle T,j) product_(i in S intersect T)(L_i/L_j).

Thus the characters alone remember a quotient, while these branch factors
retain its algebra and ramification.

Finite-flat duality gives

    pi_*omega_X = direct_sum_(S even) O_P2(|S|/2-3),       (2)

with the same characters (every character has order two). One may derive
this by applying Hom(-,omega_P2) to (1). The finite-flat Hom description is
the classical duality input, explicitly stated in
[Stacks, Lemma 49.15.2](https://stacks.math.columbia.edu/tag/0DWM).
Crepancy identifies H^0(Y,K_Y) with H^0(X,omega_X).

Consequently the dimension of the chi_S canonical eigenspace is

    h_S = binom(|S|/2-1,2) if |S|>=6, and zero otherwise, (3)

and

    p_g(Y)=sum_(r>=3) binom(n,2r) binom(r-1,2).             (4)

All sums are finite. The first decompositions are

| n | Contribution from six-line characters | From eight-line characters | p_g |
|---:|---:|---:|---:|
| 6 | 1 | 0 | 1 |
| 7 | 7 | 0 | 7 |
| 8 | 28 | 3 | 31 |
| 9 | 84 | 27 | 111 |

This reproduces the inherited Hilbert-series dimensions by a different
description and specifies exactly which quotient contributes each form.

## 3. Why the seven coordinate sections have complementary characters

The underlying adjunction isomorphism omega_X=O_X(n-6) requires a determinant
twist as a G-linearized statement. The invariant quadratic equations do not
contribute signs, whereas the ambient residue volume form transforms by
product_i epsilon_i. Thus a canonical section represented by a root monomial
x^alpha has character

    S=[n] triangle {i:alpha_i is odd}.                     (5)

It is even because sum alpha_i=n-6. The total sign of the global scalar
epsilon_i=-1 is therefore one, as required in projective space.

There is an independent complete basis check. The affine root ring is free
over k[s_0,s_1,s_2], with deg(s_j)=2, on all square-free root monomials x_A.
This follows directly by adjoining each x_i through its monic equation
x_i^2=L_i(s). Its degree n-6 piece consists of x_A times base polynomials
of degree (n-6-|A|)/2. Taking S=[n]\A recovers exactly (3).

For n=7, the seven sections represented by x_i have characters [7]\{i},
not singleton characters. In particular flipping root i fixes the canonical
line represented by x_i and negates the other six lines. Its invariant
canonical dimension is one, agreeing with deletion of root i and the
six-root K3 quotient. Omitting the determinant twist gives the wrong answer.

At n=8, distinct-root quadratics give 28 six-line characters. The eight
root squares satisfy five linear relations, leaving a three-dimensional
eight-line character space. At n=9, distinct-root cubics give 84 six-line
characters; the nine remaining character groups each have dimension three.

## 4. Each six-line character is an actual K3 double-plane quotient

For a nonempty even subset S with |S|=2r, quotient by ker(chi_S). The resulting
normal double plane Z_S has algebra O plus O(-r) and equation

    w_S^2=product_(i in S)L_i.                            (6)

Its map from X has degree 2^(n-2). This is a character quotient, different
from deleting root coordinates. For n=6, both X's resolution and the
all-six-line double-plane resolution are K3, but their map has degree 16;
no isomorphism is asserted. For n=7, the seven six-line double-plane quotients
each receive a degree-32 map from X.

For completeness the relevant resolution facts can be checked directly.
An ordinary branch crossing gives w^2=uv, an A1 node. At a triple point the
equation is w^2=uv(u+v). Blow up that base point and include the exceptional
line in the new branch divisor: if B is the original branch divisor, then
B'=strict(B)+E=pullback(B)-2E. Its half-class is rH-E. Since the blown-up
base canonical class is -3H+E, their sum remains (r-3)H. The three new branch
crossings are A1 nodes and their resolutions are crepant. This is the
standard D4 resolution of the original triple-point double cover.

Blow up all selected triple points in this manner. The resulting double
cover of the blown-up plane has O-pushforward

    O plus O(-rH+sum E).

Its H^1 vanishes: adding each E to O(-rH) has quotient O_E(-1), with no
H^0 or H^1, and plane line bundles have H^1=0. Resolving the A1 nodes does
not change H^1. When r=3 the canonical bundle is trivial, so the resulting
smooth surface is K3. More generally its canonical space has dimension
binom(r-1,2), as (3) predicts. The distinct-lines, characteristic-zero,
and no-fourfold-concurrence hypotheses are retained throughout.

Pullback of two-forms along the generically finite quotient identifies
H^0(K_(resolved Z_S)) with the chi_S eigenspace in H^0(K_Y). Resolution maps
may require a common resolution for this comparison; no unwarranted finite
morphism between the smooth resolutions is claimed. The double-plane's
trivial character contributes no two-form because H^0(P^2,K_P2)=0.

Thus at seven roots the canonical space is the direct sum of the pullbacks
of seven one-dimensional K3 canonical spaces. This is a statement about
forms and quotient maps, not a decomposition of the surface into seven
pieces, and not an identification of full integral cohomology lattices.

## 5. The other seven: a Fano quotient kills these seven forms

The seven canonical characters S_i=[7]\{i} span a six-dimensional binary
space. Their sole nonempty dependency is sum_i S_i=0: for a selected set I,
their sum is 1_I if |I| is even and 1_[7]\I if |I| is odd. In particular
there is no three-term dependency. These seven vectors are not the seven
nonzero vectors of F_2^3 and do not carry Fano-line incidence.

There is nevertheless an exact Fano transport. Bijection-label the seven
branch lines by the nonzero vectors alpha_i of F_2^3 and map

    G -> F_2^3,   [v] -> sum_i v_i alpha_i.

This is well-defined because the sum of all seven labels is zero. Its
dual is the binary simplex [7,3,4] code: every nonzero character records the
four labels on which a nonzero linear functional is one. It follows from
(2) that the quotient has no holomorphic two-forms: every nontrivial
eigensheaf contributes O(2-3)=O(-1), not a six-line character.

The orthogonal code is the Hamming [7,4,3] code. It contains the all-one
word, and the kernel subgroup is Hamming/<1>, of order eight. Each of its
seven nontrivial elements is represented by a Fano triple or its complement.
The seven canonical characters restricted to this kernel are precisely its
seven distinct nontrivial characters: their trace is 7 at the identity and
-1 at every other element. This is the reduced regular representation.

At a point lying on branch lines indexed by I, the inertia is generated by
the root sign flips in I. Single and double intersections have independent
nonzero alpha labels. A geometric triple meets the kernel exactly when its
three labels sum to zero. If no geometric triple is Fano-dependent, the
kernel acts freely on X and its node resolution. The quotient then has
K^2=16/8=2 and p_g=q=0. The companion Fano packet checks actual hourglass
labelings and their incidence; those counts are not duplicated here.

The failure is equally concrete. At a dependent geometric triple the
kernel element flips all three coordinates of the A1 cone. Writing that
cone as C^2/{+/-1}, the flip lifts to scalar multiplication by i. The
quotient singularity is 1/4(1,1), and the lifted action fixes the exceptional
P^1 pointwise. Therefore free-action conclusions cannot survive such a
triple. The character calculation p_g=0 does not by itself imply freeness.

**Classical attribution.** Kharlamov–Kulikov, *Surfaces with DIF≠DEF real
structures*, §1.1 Theorem 1.1 and Corollary 1.2 describe the diagonal
four-quadric universal canonical cover of a classical Campedelli surface
and its full sign group. Sections 1.3 and 1.5 describe seven-line covers and
the forbidden dependent-triple condition. The general carrier is recovered
literature, not a new construction here.
[Primary text](https://webdoc.sub.gwdg.de/ebook/serien/e/mpi_mathematik/2005/76.pdf).

## 6. A literal simplex, sphere, and circle-phase transition

For a complex seven-root surface, its seven canonical eigenspaces are
orthogonal under the positive Hodge Hermitian pairing. Indeed, if S!=T,
choose a deck element with chi_S chi_T=-1; invariance of the pairing makes
the cross pairing equal to its negative. Normalize a basis eta_1,...,eta_7.

For a nonzero canonical form omega=sum c_i eta_i, define

    p_i=|c_i|^2 / sum_j |c_j|^2.

The projective form space P(H^0(K_Y))=CP^6 maps onto the six-dimensional
probability simplex. Above its interior the fiber is six relative circle
phases, (S^1)^7/S^1. On a face with m positive coordinates it is (S^1)^(m-1).
Unit coefficient vectors form S^13; the square-root weights lie in its
positive real S^6 section. These are exact but different quotients.

The root deck supplies discrete sign changes in these phases. The full
continuous phase torus here acts on the coefficient vector space; it has
not been shown to arise from automorphisms of the surface. Discarding the
phases loses the projective canonical form and its periods, while retaining
the p_i preserves only normalized norm weights. This is the promised shape
transition with its precise lost coordinate, not a new arithmetic theorem.

## 7. Seven invariant alternating tensors after restriction

Let V=H^0(K_Y) for a seven-root surface. A basis tensor
eta_i wedge eta_j wedge eta_k in the exterior cube has character
S_i+S_j+S_k. Section 5 proves none is zero, so

    dim (exterior^3 V)^G = 0.

This is the abstract exterior power of the seven-dimensional space of
sections, not a wedge product of differential forms on Y; that geometric
product vanishes in degrees above the dimension of the complex surface.

Restrict to a Fano kernel H of order eight. The seven characters become
all nonzero elements of H-dual. A triple now sums to zero exactly when it
is a Fano line, giving

    dim (exterior^3 V)^H = 7.

This is true for every Fano kernel as a representation statement; the
geometric triple test is additionally needed if a free quotient is wanted.
For the labeling 1,...,7, the even representatives of Hamming/<1> are the
simplex [7,3,4] code. This makes the restricted Fano incidence explicit.

The seven-term support is a precise lead toward the alternating tensors
used in octonion/G2 descriptions, after additional real-form and orientation
choices. No such geometric realization is proved here. In particular this
does not claim G2 holonomy for a complex surface. The hostile control is
strong: the full projective sign deck preserves no nonzero alternating
three-form on V at all, despite having the same seven canonical lines.

## 8. Connection ledger, controls, and boundaries

| Source -> target | Exact operation and preserved data | Loss / sidecar |
|---|---|---|
| Root cover -> even code | Projective sign characters; full eigensheaf splitting | Branch multiplication factors are needed for geometry |
| Six-line character -> K3 | Actual double-cover quotient and canonical pullback | Degree16/32 maps need not be isomorphisms |
| Seven coordinate sections -> characters | Adjunction determinant twist | Underlying line-bundle degree alone loses the action |
| Fano labeling -> Campedelli candidate | Simplex/Hamming dual codes and group quotient | Geometric triple incidence decides freeness |
| Canonical form -> six-simplex | Squared coefficient norms | Six relative phases and periods are lost |

The independent standard-library checker compares the eigensheaf formula,
the monomial parity formula, and complete-intersection Hilbert coefficients
for n=4,...,12. It checks every projective even character, every normalized
chart basis, the unique seven-term canonical relation, both binary code
weight distributions, all 35 triples, and the restricted character traces.
It also computes the exterior-cube invariant dimensions 0 and 7 directly.
The all-order claims have the proofs above; computation supplies exact
controls rather than extrapolation.

Run `python3 -B 04-computation/six_seven_gale_cube_20261025/check_character_decomposition.py` from repository root and the same command with
`-O`. `character_results.json` and its optimized counterpart are identical.
There are no inherited enumeration filters and no producer-code imports.

The seven-root surface retains p_g=7 and the inherited finite-torification
obstruction even when it has a p_g=0 quotient. A quotient of a non-torifiable
surface can lose the birational obstruction; the reverse inference is
invalid. Binary characters are F_2 linear algebra, not an F1 descent datum.
No transport to RH, LRC, or unrestricted tournament enumeration is asserted.
