# Independent review of the wide-family bound and triangle filtration

**VERIFIED by the return lane against the displayed proofs and exact
formulas.** No change to either theorem's declared scope is required.

For [WIDE_REWRITE_SIZE](../fibres/WIDE_REWRITE_SIZE.md), the all-suffix
source formula is
(2^(j+1)−1)(Q(z)+1)2^t+1−5·2^(j−2). The bounds
H(z)−1<=Q(z)<=2H(z)−4 follow from the positive append contributions
and their reverse-complement form. Substituting the lower bound gives
the displayed positive residual after subtracting2^(2j+t), already
positive at the boundary j3,t1. The strong-word length bound limits the
number of eligible pointed z to2^j−1. Summing all permitted t before
summing j gives the convergent weighted geometric sum; no extra log X
factor remains. O(sqrt X) therefore bounds the entire stated family,
including every eligible z and the one additional affine mirrored exit.
It does not bound compositions with unrestricted other constructions.

For [TRIANGLE_FILTRATION](../fibres/TRIANGLE_FILTRATION.md), the corrected
path action is valid without any consecutiveness assumption. Rotating
the labels of a triangle untouched by the new vertices preserves all
used old-old edges because it is an automorphism of the complete old
source. Every position adjacent to a new vertex belongs to a touched
triangle or singleton, hence remains unchanged. The touched set is
invariant. Independent triangle rotations commute and act freely, so
each orbit has the asserted power-of3 size. This pays the universal
quantifier over arbitrary old-new and new-new orientations.

The C3-pair H41 constructor and sink extension give an upper gcd exactly
matching that lower factor when internal gaps vanish. The odd final
triangle uses H5; extra sinks preserve the coprime witness. For arbitrary
gaps, taking min(k,nu) matched pairs and then as many singleton triangle
extensions as remain gives the upper exponent
max(r−k−min(k,nu),0), including the saturated boundary. No equality for
general gap vectors is inferred from the seeded probes. Thus the
correction repairs the proof mechanism while preserving the stated
theorem at that checkpoint. The separate, subsequently independently
reviewed [all-k proof](../arithmetic/TRIANGLE_ALL_K_PROOF.md) now closes
that matching question. This lane also reviewed its contact grouping,
local cover counts, singleton and terminal budgets, and odd-run pairing
argument, and found no gap. The new joint enumeration and complete
within-family birth/arc rules in [TRIANGLE_METAGRAPH](../fibres/TRIANGLE_METAGRAPH.md)
were independently checked by their SCC interval-absorption proofs.
