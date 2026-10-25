# Keep the isomorphism: the boundary product becomes exact

2026-10-25. **PROVED** from the quotient-torsor description; finite
groupoid controls in `check_stack.py`. This extension integrates incoming
checkpoints f093b00695 and 2b69da9a76, read during publication. No claim
of historical novelty for the standard quotient-stack identity.

## S1. Complementary quotients over the stack

Let a finite constant abelian group G=H direct_sum K act on a scheme X.
Set R=[X/G]. Then there is a canonical equivalence

    X ≅ [X/H] ×_R [X/K],                             (1)

where the right side is the **2-fibre product**: an object contains an
object of each factor *and an isomorphism between their images in R*.
It is not the product of their coarse spaces over X/G.

**Proof.** The morphism X->R is a G-torsor in the quotient-stack sense,
even at points with stabilizers for the original action on X. Quotienting
this torsor by H or K gives the G/H- and G/K-torsors [X/H] and [X/K].
The homomorphism G -> G/H × G/K is an isomorphism because the subgroups
are complementary. Therefore the product of the two associated torsors
is the original G-torsor. Explicitly, pull back along the fppf cover
X->R: (1) becomes G×X -> (G/H×G/K)×X, the displayed group isomorphism
times X. It descends, proving (1). QED.

The standard torsor definition is [Stacks Project, Section 95.15](https://stacks.math.columbia.edu/tag/04UV); the actual complement argument
above is included. Freeness of H and K on the scheme X is not needed for
(1). If both actions are free, their quotient stacks are the ordinary
quotient schemes, so (1) has the particularly useful form

    X ≅ Y_H ×_[X/G] Y_K,   Y_H=X/H, Y_K=X/K.          (2)

This statement is valid with the original action and markings fixed.
It makes no arbitrary choice of a character basis or branch square root.

## S2. Application to the norm-seven boundary atlas

The incoming [ROOT_STACK proposition](../session63_integral_observers_root_arrangements_20261025/twists/ROOT_STACK.md)
identifies [X/G] with the simultaneous root stack of the seven branch
divisors, including their sections. On L_j!=0 it is the quotient of the
six-root algebra by mu_2^6; at r vanishing branch coordinates it retains
stabilizer mu_2^r. The atlas at a triple remains nodal, not smooth.

For the displayed norm-seven pair, H and K act freely: the branch inertia
injects into both quotient decks at every actual stratum. Consequently
(2) is an exact boundary-compatible reconstruction of X. Replacing the
base stack by its coarse plane gives instead

    X -> Y_H ×_(P^2) Y_K = Y.                        (3)

BOUNDARY B1–B2 prove that (3) is the normalization and compute its
conductor. Thus the algebraic conductor and the stack repair describe
the same comparison map. The old Laurent recovery is precisely its
restriction to the branch complement, where the root stack is the scheme
U itself.

Constant twists do not disappear from the chosen atlases and maps. The
root stack of the section d_i L_i is equivalent to that of L_i by the
specified isomorphism of line bundles with section (multiplication by
1/d_i); no square root of d_i is chosen. The twisted root covers give
different atlases and quotient maps to this common root stack. Keep
these maps, rather than just the abstract isomorphism class of the stack.

## S3. The missing r bits are an isomorphism torsor

At a geometric base point on r branch lines, an ordinary product point
gives two objects xi_H,xi_K of the root stack. The extra datum required
by the 2-fibre product is

    Isom_R(xi_H,xi_K).                                (4)

Over an algebraic closure this set has 2^r elements, a torsor under the
inertia mu_2^r. This is exactly the missing-sheet factor in BOUNDARY (10).
Since Y_H and Y_K are schemes, their objects have no automorphisms; the
2-fibre product is represented by X, not by a new stack with unwanted
stabilizers.

Over k, (4) may have no k-point despite both xi_H and xi_K existing over
k. In the explicit r=1 hostile its field of definition is Q(sqrt(2));
at the displayed triple it is Q(sqrt(2),sqrt(3),sqrt(5)). The residual
squareclass map of BOUNDARY (11) is an arithmetic coordinate description
of this very isomorphism torsor. Thus recording only the inertia *order*
does not repair rational lifting; the actual torsor and marked objects
are needed.

The incoming twist-mass theorem gives a numerical average over such
forms, but an average does not select an isomorphism for a fixed pair.
Its fixed-family extension rule also matters: over even-degree finite
field extensions quadratic twist labels collapse, whereas fresh twists
over the extension form a different family. Our fixed Q examples keep
the actual field extensions, not a fresh-twist average.

## S4. A concrete decoder on finite action groupoids

For a geometric orbit G/I, use additive notation. A 2-product object can
be represented by (x,g), where x is a point of G/I and the second point
is x+g. A change of representatives by (h,k) in H×K gives

    (x,g) -> (x+h,g+h+k)                             (5)

in the binary case. Decompose g=h(g)+k(g). Its recovered original point is

    decode(x,g)=x+h(g) in G/I.                       (6)

Equation (5) leaves (6) invariant; every orbit has the unique representative
(decode(x,g),0). This is a literal decoder, not merely a count equality.
The common source and its group action are retained.

`check_stack.py` checks this decoder for every complementary pair and
every inertia subgroup in binary dimensions 1,2,3, including cases where
H or K is not free (where one must use quotient stacks). For the actual
norm-seven pair it checks all branch strata and verifies that forgetting
the isomorphism has fibre size 2^r, while the decoder recovers every
original point. A noncomplementary hostile retains residual stabilizers
and fails the group isomorphism: complementarity is not optional.

## S5. What incoming Fourier and representation work changes

The incoming [integral Fourier criterion F2](../session63_integral_observers_root_arrangements_20261025/fourier-effectivity/FOURIER_EFFECTIVITY.md)
states that an abstract nonnegative integer histogram requires both
integral Fourier-lattice membership and positivity. Its full labelled
table is needed for arbitrary pinned continuation queries; correlations
alone forget marked phase. This agrees with (4): a mass or covariance
cannot replace the specific isomorphism torsor.

The incoming [specialization theorem](../session63_integral_observers_root_arrangements_20261025/specialization/RESULT.md)
separates invariant-module base change from multiplication. In its C3
example the invariant modules specialize perfectly, but x wedge y=3z.
Here too, the character components alone specialize cleanly while the
branch factors in multiplication determine the missing conductor. These
are exact parallel failures of forgetting multiplication, not an asserted
isomorphism between an exterior algebra and a Kummer algebra.

The [synchronized-action refinement](../session63_integral_observers_root_arrangements_20261025/cancellation/RESULT.md)
also distinguishes diagonal from independent group actions even when
their orbit ranks agree. The shared datum in (1) is explicit: the same
object of [X/G], with an isomorphism. Taking an independent product over
the coarse plane drops precisely that datum.

This closes the compatibility issue at the representation level: use
the full root-stack object and its comparison isomorphism, or equivalently
use the normalized algebra with conductor and residual classes. It does
not make nontrivial arithmetic torsors trivial, and it supplies no new
F1 descent or Collatz termination theorem.
