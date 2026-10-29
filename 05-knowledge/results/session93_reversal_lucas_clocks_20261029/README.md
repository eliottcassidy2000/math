# Lucas reflection, full reversal holonomy, and exact port clocks

2026-10-29. Source snapshot `59cb22d0e93b6c364c4ec5d2954d9fe14ed468d8`.
The session proves the recursive centralizer-generation claim, classifies the
prime-Lucas translation/reflection closure with its exact exception, and
connects both to black/blue presentations and arithmetic observation kernels.
Classical classification and source theorems are credited; no external-priority
or prime-infinitude claim is made.

## 1. The prime-degree gap is strict, and its exception is degree eleven

The [reflection-gap theorem](lucas_groups/LUCAS_GAP_THEOREM.md) states:
if a transitive prime-degree group G≤S_p contains an involution with f fixed
points satisfying **1<f<sqrt(p)**, then G contains A_p, except for the
PSL2(11) and M11 actions at degree11.

The lower inequality cannot include equality: dihedral groups at every odd
prime degree have involutions fixing exactly one point. The proof uses
Jones's CFSG-dependent classification of primitive groups containing a full
cycle, then excludes affine and projective/semilinear alternatives by explicit
fixed-point arguments. Mathieu fixed-point data are cited separately.
This is not presented as a classification-free proof.

For cyclic no-11 words of odd length n=2m+1, reflection fixes F_(m+2) words,
and the exact identity

    L_(2m+1) − F_(m+2)² = 2F_m F_(m−1)

puts the reflection strictly in the gap for n>=5. If p=L_n is prime, the
inherited golden reader identifies the words with F_p. Arithmetic translation
and transported word reflection therefore generate:

- PSL2(11) at n=5;
- S_p when n≡11 mod12;
- A_p for all other odd n>=5 with prime L_n, apart from n=5.

This closes the prime-Lucas classification part of the earlier
`lucas_mirror_gap_spaces_20261028/SYNTHESIS.md` Open C. It does not establish
infinitely many Lucas primes or classify composite-size carriers and restricted
translation menus.

Explicit short words independently certify A29, S199 and A521 using classical
Jordan small-cycle theory rather than CFSG. At degree199, one translated
reflection followed by a specified power isolates a23-cycle. The notes contain
all permutations, exponents and cycle checks.

## 2. A larger symmetry group can force a larger decoding frame

These alternating/symmetric groups cannot themselves preserve a tournament
orientation on the state labels: they are transitive on ordered distinct
pairs. A decoder must retain extra state.

If a finite G-carrier decodes equivariantly to tournaments on the same named
point set, each decoder-state stabilizer has odd order. Hence every decoder
orbit has size divisible by the full2-part of |G|. The concrete lower bounds
are2^24 for A29,2^194 for S199, and2^517 for A521.

The contract is full closure under the entire unguarded group. These are not
bounds on a single certificate or on a carrier permitting only guarded native
operations, and attainment is not proved at those larger degrees.

At the exceptional degree11 the exact minimum is **12**, attained by the
Paley11 orbit. Odd stabilizers and simplicity of PSL2(11) rule out a smaller
orbit; direct enumeration confirms660 group elements,12 images and a55-element
stabilizer. This explains why the older twelve-frame repair is genuinely
exceptional rather than a template that stays small at every prime period.

## 3. Triangle errors generate the full recursive centralizer

Let W_k be the iterated C3 tournament on N=3^k vertices, with digitwise-negation
anti-involution alpha. Write A_k=Aut(W_k) and C_k=C_(A_k)(alpha).

[Full centralizer generation](holonomy/FULL_CENTRALIZER_GENERATION.md) proves

    C_k ≅ product_(j=0)^(k-1) A_j,
    |C_k| = 3^((3^k−2k−1)/4).

The square-root transport between anti-involutions has triangle errors in C_k.
For every k, **exactly k(k−1)/2 individual based triangle errors suffice to
generate all of C_k, and no smaller generating family can do so**.

The proof gives explicit paired subtree rotations. A top rotation and a local
zero-prefix rotation produce each required triangle; the centralizer recursion
supplies the remaining factors. Abelianization proves the sharp lower bound.
This is actual generation, not merely normal generation or numerical evidence
from nine vertices.

The result concerns the natural pair-local square-root convention. Retaining a
global base mark still supplies flat transports, as proved in Session90. There
is no contradiction: a different convention carries different frame data.

## 4. Exactly what port information these loops can erase

The [port/LCA theorem](ports/PORT_ORBITS_AND_LCA.md) gives a complete tuple
invariant for C_k. A port's first nonzero ternary digit determines its block
and side. Negating the suffix on the second side folds each paired block into
a common A_j action. Within that joint folded list, longest-common-prefix
depths and directed cyclic arcs completely determine the ordered-port orbit.
A recursive construction supplies the actual transporter, not just its order.

In particular,

    one-port orbit count = 2k+1,
    ordered two-port orbit count = 8k²+1.

The minimum number of named anchors killing every centralizer ambiguity is

    (3^(k−1)−1)/2,       k>=1.

At27 vertices the group has order243, needs three group generators, but needs
four named port anchors for a pointwise base. These are distinct costs. Naming
the fixed root does not remove this ambiguity.

Individual port orbits and original cross-block arcs are insufficient: at
order9 the pairs(3,6) and(3,8) agree on both, but their folded relative alignment
differs. This is a precise positive LCA repair for this recursive family, not a
universal theorem that arbitrary LCA summaries reconstruct hidden relations.

## 5. Blue classes are flat even when named paths are not

The [black/blue bridge](metagraph/BLUE_QUOTIENT_FLATNESS.md) fixes the types.
For an anti-involution s, let X_s be the Hamilton paths reversed by s and let
C_s be its centralizer. The blue presentations are X_s/C_s. Any transporter
between two reversal marks induces the same map on these quotients, so the
blue-presentation transport composes **strictly and independently of choices**.

A nonidentity triangle error nevertheless moves every named full path, because
automorphisms act freely on ordered Hamilton paths. It can change the labelled
native partner while preserving its isomorphism class and blue presentation.
A shared port retains relative path position: quotienting a path and a port
together is not equivalent to quotienting each marginal separately.

Exact controls:

- Paley7 has9 fixed paths,3 blue presentations and trivial triangle holonomy,
  despite its centralizer having order3.
- W_2 has33 fixed paths,11 blue presentations and full centralizer holonomy3.

Thus flatness is not purity of the blue line, and full holonomy is not purity
of the black line. A marking loop is not an adjacency-changing native move.
It does, however, transport an entire black-core certificate when all path
and port guards are relabelled together. The unique fixed pole stays the same
named target. This reuses the inherited prescribed-black-core theorem without
claiming to prove general black-core connectivity from holonomy.

## 6. Trace recurrences retain less than a chosen arithmetic clock

The [trace-clock theorem](arithmetic/TRACE_CLOCK_CARRIERS.md) works concretely
in O=Z[phi] with the norm-one unit a=phi². For x=u+v phi, two adjacent readings
are

    Tr(x)=2u+v,       Tr(ax)=3u+4v.

Their determinant is5. Over O/(m), traces sampled every d steps have exact
whole-history kernel size gcd(m,5F_(2d)). Two samples already carry everything
that the entire selected-clock history can reveal. Several strides replace d
by their gcd. Coprime clocks can repair decimation loss, but not ramification
at5.

The distinction is explicit:

- Modulo3, stride2 loses information that stride1 retains.
- Modulo33, one norm-one trace can have four preimages, from independent local
  conjugate choices rather than one global branch bit.
- Modulo5, five distinct nilpotent phases have the same entire trace history;
  no additional trace stride recovers them.

The same observations satisfy an exact quadratic conic identity. Away from5,
stride-one norm-one states biject with t²+z²−3tz=−5; at5 this degenerates,
matching the loss. The clock is specified throughout. Simultaneously conjugating
both the source and clock preserves every trace, so an unnamed conjugate cannot
be selected by this recurrence alone.

Finite-field Frobenius gives order bounds, but actual order must be tested.
At the inert prime47 the norm-one group has48 elements while a has order16.
The note supplies a certified order-reduction algorithm after factoring the
bound, not a theorem about prime density or infinitely many primitive values.

An older complementary source is
`THM-3312-exceptional-quadratic-trace-norm-and-cofactor-antidescent.md`: its
quadratic trace/norm passport descends while a cofactor branch does not.
The common map is genuine quadratic trace/norm; no identification of those
Jacobian cofactors with tournament ports, or Jacobian-conjecture consequence,
is asserted. The unspecified recent surd family was not identified, so the
explicit golden model is not silently claimed as a replay of that family.

## Verification and limits

[Independent audit](audit/INDEPENDENT_AUDIT.md) checks the classification
assumptions, projective descent, Jordan words, recursive triangle identities,
blue quotients, literal port orbits, minimum bases and trace-conic identities.
The generation checker uses SymPy's exact permutation-group routines for the
larger finite controls. Other scripts use the standard library. All are run
normally and with `python3 -O`; explicit checks remain active in both modes.

The main closures are precise: prime-Lucas full-group classification, full
recursive triangle generation, and exact port-observation laws. Prime
infinitude/distribution, arbitrary-tournament centralizer generation, general
black-core connectivity and the repository's major arithmetic conjectures are
not solved by these results.
