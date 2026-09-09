# Independent acceptance of the geometric reconstruction strengthening

**Status: PASS, independent analytic referee.** The
[geometric reconstruction theorem](continuing14_20260908_geometric_fibre_reconstruction.md)
is accepted for every supplier degree q>=3. In particular the two frozen
degree-fifteen examples are nonisomorphic even after algebraic closure of
C(T). The narrower original theorem was correct and is strengthened, not
retracted. This referee is report-only; it is not another finite engine.

The key repair to the original proof is necessary: degrees of closed
punctures all become one over an algebraic closure. The new proof recovers
the constant finite punctures by a valuation invariant of their unlabelled
configuration, rather than treating their original arithmetic labels as
preserved without proof.

Fix a valuation extension at T=infinity, choose t with T=t^(-N), where
N=2q+1, and scale by w=tz. Every moving root of Q(z)=T has valuation -1.
The scaled polynomial reduces to a_N w^N-1 with distinct nonzero roots.
Its unit discriminant and the nonnegative valuations of root differences
force those differences to be units, so all scaled reductions are distinct.
Thus the q finite constant punctures form one scale and the N moving
punctures plus infinity form the other, with no omitted internal collision.

For an unordered quadruple, degeneration of its cross ratio into
{0,1,infinity} is independent of the ordering. A quadruple degenerates
exactly when it has two points in each scale. This can be checked directly
using the producer's formula: the numerator has valuation -1 and the
denominator -2; including infinity gives the same limit. Every other
distribution has four distinct reductions in either z or tz, so cannot
degenerate. This proves both directions of the criterion.

Consequently a fixed finite puncture lies in
(q-1)*binom(2q+2,2) such quadruples, whereas an outer puncture lies in
(2q+1)*binom(q,2). Their ratio is 2(q+1)/q. The two counts differ for
all q>=2, and are 720 and315 in the explicit pair. An isomorphism over
the algebraic closure preserves the actual cross-ratio values, hence
these counts for the SAME chosen valuation. No extra assumption that
the isomorphism respects the valuation or a source model is required.

For q>=3, three constant finite points therefore have three constant
finite images. Their unique Mobius map has coefficients in C. Infinity
must map to a constant outer puncture, of which infinity is the only one.
The map is affine. Evaluating Q_2 o phi-Q_1 at a moving root proves the
polynomial identity because that root is transcendental over C. Conversely
the affine identity sends both derivative-root and moving-root punctures
correctly. This independently establishes the full geometric iff.

The proof pays its boundaries. Equal-sized clusters can be exchanged, and
two constant fixed points do not determine a constant Mobius map. The
q=2 example is an obstruction to this proof step, not a counterexample
to geometric reconstruction in that supplier class. No total-fibration
classification or arbitrary target-parameter change is inferred.

Accepted frozen supplemental pins:

    report 8d41e167994cd06840df49dfcfea0fa4fdee5eb3ded4ab743e43fcb61275a0d2
    source 731299dc13540c4861af3b77d134fa9d4a401dbb47f5e1e1bf56d326d633d661
    output b78deb01eb9dd43a95ea4483fc20e86d60f6df7ba7aa1f974107d35353a95f5a
    cert   07bde5bc5bf523d16040aa13cb97db8fb39fad43a1db640d7c91bda4ee3133a2

The supplemental producer's 18,579 always-active gates test all 9,255
declared quadruples and 27,765 pairings at q=2,3,7, discriminants and
hostiles. Those finite controls supplement the separate all-degree
valuation proof; they do not supply its quantifiers. Normal and optimized
reproduction is required before the packet is filed.
