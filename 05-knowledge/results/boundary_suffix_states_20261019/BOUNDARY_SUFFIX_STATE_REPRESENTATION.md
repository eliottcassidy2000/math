# Finite suffix states retain exterior profile weights

19 October 2026. Status: PROVED and root-reviewed; normal and optimized exact controls agree. This refines the bounded-chain mixture theorem. It gives a bounded interface for actual profile multiplicities, not a classification of which multiplicities arbitrary passive exteriors can realize.

## 1. What the convex hull forgets

The [bounded-chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md) replaces arbitrary exteriors, after conditioning on their extension orders, by finite two-chain atoms. Its exact identity preserves all core-word counts. Taking the unrestricted atom convex hull loses two kinds of information:

1. The multiplicities of lower and upper profiles factor, before extension-count normalization.
2. Profiles on one side come from orders of one fixed active exterior interface. Their weights cannot be chosen independently.

This note retains both coordinates. The finite state representation is exact for each given ambient poset, even when the passive part has arbitrarily many vertices.

## 2. Active and passive lower exterior

Use the same full critical-root-core decomposition P=Q union A union B and maximum incomparability degree k. The roots are oriented critically, so D(x) is a subset of D(y) and U(y) is a subset of U(x); this is used for the one-sided active-size estimate below. Let

    T={a in A : a is incomparable with at least one q in Q},
    A0=A minus T.

Every member of A0 lies below every core vertex. No t in T can lie below a member of A0, because t<a0<q for every q would contradict t's membership in T. Thus all active/passive relations have the form a0<t or incomparability; never t<a0.

Let delta(q)=k-d_Q(q) as before, and put r=max_(q in Q) delta(q). Then

    |T| <= sum_(q not roots) s_A(q)
         <= sum_(q not roots) delta(q).

For a full critical-root core, the only core vertices possibly incomparable with A are the lower-exclusive and common vertices, at most k-1 of them. Each has delta<=k-1, so |T|<=(k-1)^2. The corresponding bound holds upstairs. The combined active size is at most sum_(q not roots)delta(q)<=2(k-1)^2. These estimates count incidences conservatively and are not asserted sharp.

The active poset T has width at most k: an antichain of k+1 active vertices would give each member k incomparable active neighbors and at least one incomparable core neighbor, violating k-thinness. Most importantly, since no active vertex can lie below a passive one, the active maximal vertices of A minus an active final filter S are exactly the maximal vertices of T minus S. They depend only on the bounded active poset.

## 3. Exact tied suffix counts

Let S be an upper filter of T, and define

    H(S)=e(A minus S).

Then S is also an upper filter of A. For each particular linear extension w of S, the number of extensions of A ending in precisely the word w is H(S): concatenate any extension of A minus S with w. Thus different legal suffix orders with the same set S have exactly equal multiplicities, regardless of passive attachments.

Let M(S) be the maximal elements of T minus S. For each a in M(S), put S+a=S union {a}. Define

    J(S)=H(S)-sum_(a in M(S)) H(S+a).                         (1)

If A minus S is nonempty, J(S) counts its extensions whose last vertex is passive. If it is empty, J(S)=1. In particular J(S)>=0. This identity is the exact partition by the last vertex; it does not assume uniform probabilities on the maximal active choices.

Interpret a path as successively deleting last active letters. Starting at S=empty, one may either stop with multiplicity J(S), or delete a in M(S) and continue. Different paths reaching the same filter share the same H(S) and J(S). In normalized form the probability of choosing a next, conditional on a specified suffix already removed, is H(S+a)/H(S); its stopping probability is J(S)/H(S). Products along a path telescope.

For two simultaneously available active vertices a,b, this entails the exact diamond identity

    [H(S+a)/H(S)] [H(S+a+b)/H(S+a)]
      = [H(S+b)/H(S)] [H(S+a+b)/H(S+b)].                     (2)

It is a statement about an actual common suffix set, not an assumption of independence.

## 4. The profile is determined after at most r steps

For an exterior extension sigma, examine its final consecutive run of active vertices and truncate that run at length r. Read it backwards as a deletion path. If it stops before depth r, the preceding vertex is passive, or the entire extension has been exhausted. A passive vertex is below every core vertex and therefore resets every terminal gap to zero.

For each q, alpha_sigma(q) is the length of the trailing block of letters incomparable with q. This block has length at most s_A(q)<=delta(q)<=r. Thus the truncated active suffix determines the entire profile exactly. If the path reaches depth r, no earlier information is needed even when all r retained letters are incomparable with a particular q, since its actual gap cannot exceed r.

THEOREM 1 (exact finite suffix-state representation). Build the deletion tree of legal active suffix words to depth r. A path reaching filter S at depth less than r contributes J(S) to its terminal state. A path at depth r contributes H(S). Map each terminal word to its profile by counting the trailing incomparable block for every q. Summing these terminal weights over words giving alpha yields exactly

    m_A(alpha)=#{sigma in LE(A): alpha_sigma=alpha}.           (3)

All H and J are integers, terminal weights are nonnegative, and their sum, including path multiplicity, is H(empty)=e(A). The analogous upper construction uses initial active runs, or equivalently applies this theorem to the dual.

For r=0 the sole terminal state has weight e(A) and profile zero. If the active interface has fewer than r vertices, exhausting it is included in J(S), so no artificial padding or passive ordering is introduced.

This is exact compression of the passive influence into finitely many extension counts H(S). It does not yet characterize the allowed values of those counts.

## 5. A bounded number of coefficients

At every node there are at most k available active maximal elements, so the deletion-word tree has at most 1+k+...+k^r nodes. The shared-set representation is substantially smaller.

By the classical chain-partition theorem for a width-at-most-k poset ([Dilworth, Theorem1.1, Annals of Mathematics51 (1950),161–166](https://mathweb.ucsd.edu/~asuk/dilworth.pdf), primary theorem checked19 October2026), partition T into at most k chains. A filter S meets each chain in a final segment and is determined by the k segment lengths. Therefore the number of distinct filters of size h is at most binomial(h+k-1,k-1), and

    #{S upper filter of T : |S|<=r} <= binomial(k+r,k).        (4)

The elementary counting deduction is independent of the choice of chain partition. With r<=k-1, there are at most binomial(2k-1,k) H-coefficients. For k=7 this is at most 1716 per side. Word paths still carry profile labels and cannot all be identified merely because they share a set; only their H weights are tied.

## 6. Lower capacity bounds on every active branch

Let a belong to M(S). It is maximal in A minus S. Deleting a maps every extension of A minus S to one of A minus (S+a). For any fixed image, the legal reinsertion slots for a lie after its last predecessor. There are at most d_(A minus S)(a)+1 such slots. Appending a is always allowed, so the map is surjective. Hence

    H(S+a) <= H(S) <= [d_(A minus S)(a)+1] H(S+a).            (5)

Let t_Q(a)=|I(a) intersect Q|, which is at least one, and let i_S(a)=|I(a) intersect S|. All B vertices are comparable with a, so the ambient degree bound gives

    d_(A minus S)(a) <= k-t_Q(a)-i_S(a).

Consequently the following finite-interface inequality is valid:

    H(S) <= [k-t_Q(a)-i_S(a)+1] H(S+a).                      (6)

The bracket is a positive integer and is at most k. For the depth-r representation, these branch inequalities are required only at |S|<r, so no coefficient beyond depth r is needed. Thus every legal active transition has probability at least 1/k, with the sharper displayed denominator when more incidences are already known. Along a specified length-h active suffix, its probability is at least the product of these reciprocal bounds, and in particular at least k^(-h).

This is a LOWER branch-capacity constraint, complementing the upper increment caps from packet147. Arbitrary nonnegative profile mixtures need not satisfy it. The stopping weights J(S) may still vanish; no positive stopping probability is claimed.

When |M(S)|=k and all transition denominators equal k, their sum is already at least one. Equation(1) then forces zero stopping probability and every transition probability exactly1/k. More generally the sum of reciprocal bounds cannot exceed one; a proposed active interface violating that necessary condition is excluded without counting its passive part.

## 7. Exact rank-one profile multiplicities

Define m_A(alpha) and m_B(beta) by the finite state representations. For each core word w, let T_(alpha,beta)(w) be its lift count in the corresponding two-chain atom. Then

    C_P(w)=sum_(alpha,beta) m_A(alpha)m_B(beta) T_(alpha,beta)(w). (7)

The multiplicity array before normalization is therefore rank one, over nonnegative integers. If pi_(alpha,beta) denotes the aggregated normalized mixture weight of that atom, then

    pi_(alpha,beta)/e(C_(alpha,beta))
       = m_A(alpha)m_B(beta)/e(P).                           (8)

Its two-by-two minors vanish. The normalized mixture weights pi themselves are generally not rank one; dividing by the atom extension total is essential.

Actual positive supports form a rectangle: supp(m_A) cross supp(m_B). In fact their coordinatewise maxima RECOVER the actual exterior allocations exactly:

    max_(alpha in supp(m_A)) alpha(q)=s_A(q),
    max_(beta in supp(m_B)) beta(q)=s_B(q).                   (9)

To prove the first equality, let I=I_A(q). It is an upper filter of A, and its complement is the strict A-downset of q. Concatenate an extension of A minus I with an extension of I. Its terminal q-gap is exactly |I|=s_A(q), while no gap can be larger. The upper proof is dual. Different q may attain their maxima in different profiles; no simultaneous-maximizing profile is asserted.

The support maxima consequently satisfy max alpha(q)+max beta(q)<=delta(q). Conversely, for a product support, imposing alpha(q)+beta(q)<=delta(q) on every positive pair is equivalent to this coordinatewise allocation bound, since independently maximizing profiles can be combined.

There is also an exact and quantitative maximal-gap mass formula. With s=|I|<=r,

    Prob[alpha(q)=s]=e(A minus I)e(I)/e(A)
                   =H(I)e(I)/H(empty)
                   >= e(I) k^(-s).                        (10)

Indeed alpha(q)=s occurs precisely when the final s vertices are I. Each of its e(I) suffix orders has the same multiplicity H(I); the branch-capacity bound applies along every such suffix. This constrains the mass at a support maximum, not just its existence. For s=0 the formula gives1.

The shared-set ties, flow identities and branch lower bounds constrain each factor m_A and m_B before forming this product. Factorization of the coarser increment array U_u V_v alone does not retain them.

## 8. A sharp passive-bias control

Take a four-antichain core Q={x,y,q,r}, with roots x,y. The active lower exterior consists of incomparable a,b, with a below x,y,r and incomparable with q, and b below x,y,q and incomparable with r. Take a passive chain c1<c2<c3, all below Q and below a, and all incomparable with b. There is no upper exterior.

This nine-vertex poset is five-thin: b has exactly five incomparable neighbors (a,r,c1,c2,c3), and every other vertex has fewer. The lower exterior has five extensions. Its coefficients are

    H(empty)=5, H({a})=4, H({b})=1, H({a,b})=1,
    J(empty)=0, J({a})=3, J({b})=0, J({a,b})=1.

The two legal complete active suffix words (a,b) and (b,a) each have multiplicity1, exactly as the shared-set rule demands. The remaining three extensions stop after suffix(a). The profile “q has gap1, r has gap0” has multiplicity4; the opposite profile has multiplicity1. Thus passive attachments can create bias, but not arbitrary weights. The first transition selecting b attains the lower bound1/k=1/5 exactly; the next transition after deleting a attains the sharper bound1/4.

This example rules out both assuming uniform active choices and allowing an unconstrained profile distribution concentrated on one of these two profiles. The independent standard-library controls verify these values, a depth-one truncation variant, a non-chain lower interface, its upper dual, and the empty-active case.

## 9. Remaining realizability question

For a fixed bounded active poset and fixed incidence with Q, the inequalities above are necessary for its H array. They are not asserted sufficient for realization by a finite passive poset. In particular, integrality, shared extension counts of overlapping deletions, and any additional inequalities from the passive attachments remain available constraints.

The next useful mathematical target is to characterize or separate this finite H-array region, rather than treating arbitrary atom mixing as actual exterior realizability. A global seven-thin argument still needs coverage of active interfaces and core morphologies; this bounded representation alone does not prove balance or Gold.

## 10. Exact controls and reproduction

Run from this directory:

    python verify_suffix_state_representation.py
    python -O verify_suffix_state_representation.py

Normal and optimized receipts match on209 explicit guards. The verifier uses only the Python standard library and checks selected exact fixtures, not a census. Its SHA-256 is 4dab33a999ce590d946b334517c7703a05c3dd78c60ed08580429d0bdf1fc481. The general representation and capacity claims are established by the proofs above; the controls also reject confusing normalized atom weights with the rank-one multiplicity matrix.
