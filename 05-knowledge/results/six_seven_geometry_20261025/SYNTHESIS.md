# A common action behind several sixes and sevens

2026-10-25. **PROVED** claims below refer to the scoped proofs and independent
audits linked from the README. Finite enumerations are **FINITE-EXACT**;
classical surface constructions are **CITED / RECOVERED**. No novelty claim.

## 1. The numerical threshold has two compatible explanations

For n distinct spanning branch lines on P2, with no four concurrent, the
projective square-root surface is a complete intersection of n-3 quadrics
in P^(n-1). Its minimal resolution is crepant, so

    K=(n-6)H,       H²=2^(n-3).

Six roots balance the canonical degree; seven roots make the root coordinates
canonical sections. Equivalently, projective sign characters are even subsets
S, with canonical eigenspace `H0(P2,O(|S|/2-3))`. The first nonzero forms
therefore occur on six-subsets. At n=6 there is one; at n=7 there are seven.
This also predicts `p_g(8)=28+3=31` and `p_g(9)=84+27=111`.

These are recovered in the concurrent [Session60 deck packet](../session60_six_seven_geometry_20261025/deck-k3/DECK_CHARACTERS_AND_K3_QUOTIENTS.md)
and independently replayed [here](characters/CHARACTER_SEXTICS.md). Seven
canonical characters span a rank-six space for the **full** deck group.
They become the seven nonzero rank-three characters only after choosing
the code subgroup described below. The choice is real data.

## 2. Deleting a root exposes a hidden cube of curves

Fix a deleted line contained in t triple intersections. On the six-root K3,
the omitted square root branches over a divisor of class2H and arithmetic
genus17. Each triple creates a repeated restriction of two branch forms,
and each repeated restriction creates an independent sign choice.

| t | Branch components | Genus per component | Intersection graph | Actual smooth quotient of seven-root Y |
|---:|---:|---:|---|---|
| 0 | 1 | 17 | one vertex | K3 |
| 1 | 2 | 5 | eight intersections | Bl8(K3) |
| 2 | 4 | 1 | square with edge multiplicity4 | Bl16(K3) |
| 3 | 8 | 0 | cube with edge multiplicity2 | Bl24(K3) |

For the actual hourglass, central deletion has t=3 and every peripheral
deletion has t=1. The cube is an intersection calculation on eight explicit
conics, not a resemblance in a drawing. Its matrix is `-2I+2A(Q3)`.

The seven nonconstant Walsh modes are indexed by the seven nonempty supports
of three bits. The three singletons have eigenvalue0 and are actual
hyperplane divisor relations; the three pairs have eigenvalue-4; the full
support has eigenvalue-8. The omitted empty support has eigenvalue4.
Thus the component classes span rank5, not7 or8. In the inherited monoid P2
convention, precisely the three singleton supports are rank-zero points:
the same labels here select divisor relations. The representation and its
marked directions make this connection precise.

[Proof and independent252-deletion checks](root_deletion/DELETION_AND_CUBE.md).

## 3. One group joins the cube, Fano geometry and canonical sections

Label the seven roots by `v_i in F2^3 minus {0}`. Evaluation words
`(a dot v_i)_i` form the simplex code S of dimension3 and weight4. It embeds
in the projective deck group G of dimension6. The Hamming syndrome gives

    0 -> S -> G -> F2^3 -> 0.

This S is also the kernel of the earlier repository Fano point-to-line
sign-gauge map. The arithmetic gauge kernel has acquired an actual surface
action. On the seven canonical sections its representation is the regular
representation of F2^3 with the trivial character removed.

A nonzero code element fixes a root point exactly when an actual concurrent
triple is one of the chosen Fano lines. For the hourglass, write c for the
central label and d1,d2,d3 for the three opposite-pair label differences.
Then `d1+d2+d3=c`. Its action on the eight central conics is the translation
map `a -> (a dot d1,a dot d2,a dot d3)`. Consequently

    no bad concurrent triple
       iff d1,d2,d3 are a basis
       iff the action on the eight conic components is regular
       iff S acts freely on the seven-root surface and its resolution.

This is the strongest unification in the packet. The eight-component regular
representation and its seven nontrivial characters now belong to one group.
The intersection eigenvalues still require the three marked cube directions;
we do not identify a section with a curve or an integral lattice vector.

Exactly16 of the30 Fano structures are free for the hourglass. The quotient
is a smooth minimal Campedelli surface with `p_g=q=0`, `K²=2`: the group has
order8, so the covering surface's `chi=8,K²=16` divide by8. The remaining
12 markings have one bad triple, and2 have three. Their quotient resolutions
have K²0 and-4 respectively because each bad triple contributes two
`1/4(1,1)` singularities. The final incoming
[Session60 theorem](../session60_six_seven_geometry_20261025/seven-lines/RESULT.md)
classifies their minimal models: properly elliptic for one bad triple;
Enriques after contracting four (-1)-curves for three. Its bicanonical descent
argument also shows that the16 free choices give exactly two complex
Campedelli isomorphism classes, with D12 orbit sizes4 and12. The census and
classification are recovered concurrent work; the regular-conic-action
criterion is our additional common mechanism.

The full independent finite universe contains5596 partial Steiner systems
on seven labelled points; each avoids at least8 of the30 Fano structures.
This proves a general existence deduction for every seven-line arrangement
in the stated scope. It counts quotient maps, not isomorphism classes.
[Proof, full census, and common-action lemma](fano_codes/RESULT.md).

## 4. The failed seven-support reflection has two positive repairs

First, resolve quadratic Cremona inversion by blowing up the three coordinate
points of P2. The resulting degree-six del Pezzo has six boundary components
forming a hexagon. Delete its six torus-fixed points: the remaining six
rank-one strata and dense torus have seven literal A1 points. Blowdown
identifies their labels with the seven nonempty supports; the geometric
involution exchanges complementary proper supports and fixes full support.
The order failure is now paid for by exceptional directions: singleton labels
have changed rank from zero to one. The full blowup has13 strata.

Second, on the hourglass itself, `(E,x,z)->(E,-x,-z)` induces the same
3+1+3 permutation after an adapted labelling. It exchanges Fano markings
rather than preserving one Fano plane. The old four-line incidence failure
is compatible with a genuine geometric automorphism on this richer carrier.

There is an exact six-element bridge as well: the peripheral hourglass
slopes `±(1,0),±(0,1),±(1,1)` are the primitive rays of the del Pezzo fan.
Their dihedral symmetry and central inversion agree. Affine line offsets
and concurrency are not fan incidence: opposite toric boundary curves are
disjoint, while the corresponding projective branch lines meet.

The six full coordinate flags are also the six transitive tournaments on
three labels; adjacency of torus-fixed flags is a single adjacent swap,
equivalently a single arc reversal. The two cyclic tournaments are outside
this carrier. Divisibility and the earlier ten-divisor complement laws keep
their existing flags and twists; no arbitrary tournament predicate transfers.

## 5. Quotienting unlocks a symmetry that individual roots forbid

The central six-root K3 has a degree16 map to the double-plane K3 with
branch `(E²-x²)(E²-z²)(E²-(x+z)²)`. On that double plane an explicit
coordinate change turns a paired-line Cremona involution into

    [X:Y:Z] -> [YZ:XZ:-XY],
    F(c(X,Y,Z))=(XYZ)²F(X,Y,Z).

After resolving the three chosen crossings and the twelve remaining nodes,
one lift has eight fixed points and preserves the two-form; the other is
free and gives an Enriques quotient. The four base fixed points are checked
to miss the branch divisor. This is an actual hourglass-derived instance
of a classical Cremona/K3 construction.

The same transformation cannot lift to the six-root or seven-root field.
It sends one defining radicand ratio to

    X(Y-2Z) / (Y(X+2Z)),

which has odd valuation on X=0. Every original root squareclass has even
valuation there. The total product cancels the new coordinate factors to
a square, while the individual-root extension does not. This identifies
the exact information a quotient forgets in exchange for a new symmetry.

The degree-six del Pezzo also has seven anticanonical sections (the center
and six vertices of a lattice hexagon). A smooth double cover branched in
`|-2K|` gives a K3 polarization of square12 and genus7, with eight sections.
This is a further concrete six/seven mechanism, distinct from the seven
canonical sections on the general-type surface. Their graded rings prevent
an accidental identification: degrees6 versus16, quadratic relations9 versus4.
[Explicit models, hostiles and proof](support_reflections/RESULT.md).

## 6. What remains open, and the next decisive probes

- **Recovered closure:** the16 free maps have two complex isomorphism classes;
  the nonfree minimal models are properly elliptic or Enriques. These use
  additional bicanonical, pencil and contraction arguments in the incoming
  theorem, not K² alone. Abstract isomorphism classification of the properly
  elliptic quotients remains **OPEN**.
- **OPEN:** identify integral saturation and arithmetic/Galois descent of
  the branch-cube lattice. The current result is a rational span and an
  algebraically closed component configuration.
- **Recovered closure:** the incoming
  [exact-rank theorem](../session60_six_seven_geometry_20261025/outer-picard/EXACT_RANK_COROLLARY.md)
  proves outer K3 rank19 and hourglass rank58, combining an explicit
  nontorsion elliptic section with two-prime discriminant bounds. Thus the
  seven deletions form exactly two complex isomorphism classes, of ranks20
  and19. Our explicit curves now sit inside a known-rank lattice; its integral
  saturation, rather than the earlier numerical rank choice, is the next gate.
- **Unchanged OPEN problems:** distinct positive integer-square hourglasses,
  full magic squares of squares, and repository LRC(14).

The productive research move was to retain an initially numerical lead long
enough to test three changes of carrier: resolve the reflection's directions,
restrict the deck group by a code, and keep the branch divisor under deletion.
Each produced a map and a loss ledger. The stopped routes now have exact
witnesses: forgotten branch nodes, the canonical determinant twist, failed
Fano normalization, and the odd new-axis valuation.
