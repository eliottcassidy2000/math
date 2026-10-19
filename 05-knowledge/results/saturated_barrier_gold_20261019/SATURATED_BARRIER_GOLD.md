# Saturated barriers, a new induced Gold obstruction, and two free ports

19 October 2026. Status: PROVED by elementary structural reductions and a portable exact finite certificate, independent review PASS. No general seven-thin, 1/3–2/3, or Peczarski Gold theorem is claimed. No historical novelty is asserted.

## 1. The new scope

Let Q be the twelve-vertex poset with strict predecessor masks

    (0,0,0,0,7,13,14,12,206,5,685,535).

This is the second no-good-pair seed in the [height-three witness note](../height_three_seven_thin_no_good_pair_20261019.md). Its unweighted half pair was already known. What follows concerns arbitrary ambient attachments, whose restriction to Q need not be uniform.

Fix the same two questions throughout:

    first 7<9; then 4<5, on either first branch.             (1)

For a poset P, write E=e(P). If a first-answer branch has b extensions and a complete answer path inside it has c extensions, its Gold margin is E-b-c.

**Theorem A (induced-pattern exclusion).** Every finite seven-thin poset containing an induced copy of Q has all four margins of (1) strictly greater than E/9. More precisely, its smallest normalized margin is at least

    50695/453252 > 1/9.                                     (2)

The constant in (2) is sharp for this fixed tree within this induced-pattern family. The first question also satisfies the sharp bounds

    66465/151084 <= Pr(7<9) <= 84619/151084.                 (3)

This requires no designated critical roots, no full-root-core assumption, and no bound on the number of exterior vertices. In particular, an inclusion-minimal seven-thin counterexample to weak Gold cannot contain Q as an induced subposet.

There is an independent strengthening beyond bounded incomparability degree. Put

    S={1,4,5,6,7,9}.

**Theorem B (two unbounded attachment ports).** Suppose a finite poset P contains Q induced, every vertex outside Q is comparable with every vertex of S, and each of 0,3,8,11 has at most one incomparable neighbor outside Q. Then the tree (1) is strict Gold, and every margin is strictly greater than

    (386/6245) E > (3/50) E.                                (4)

There is no bound on the number of exterior incomparable neighbors of 2 or 10, nor on exterior sizes or internal degrees. The constant 386/6245 is the sharp infimum for this fixed tree in this larger family. Section 8 gives actual finite posets approaching it.

Thus the transfer does more than optimize the already closed ten-core tree menu: it closes a different induced-core morphology and also an unbounded-degree attachment family.

## 2. Inheritance and the structural change

The [bounded-chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md) preserves every core-word lift count after conditioning on exterior extension orders. The [three-tree cover](../three_tree_gold_cover_20261019/THREE_TREE_GOLD_COVER.md) uses that reduction for a particular ten-vertex full critical-root core. The [suffix-state](../boundary_suffix_states_20261019/BOUNDARY_SUFFIX_STATE_REPRESENTATION.md) and [attachment-correlation](../boundary_attachment_correlation_20261019/PASSIVE_ATTACHMENT_HULL_AND_FKG.md) notes explain information lost by arbitrary convexification.

Here a single tree is valid throughout the relevant cone. Consequently shared suffix counts, rank-one multiplicities, and nonlinear correlation inequalities are not required for this transfer. This is a stronger uniform certificate, not a claim that those realizability restrictions disappear in general.

The new structural object is a connected, dominating set of degree-saturated core vertices. It forces the two-sided exterior decomposition from induced containment alone. The second operation replaces an unbounded terminal-gap coordinate by a positive combination of a bounded-cap model and a vertex-deletion model. These are exact maps of extension-event counts.

## 3. Saturated-barrier lemma

Let Q be an induced subposet of a k-thin P. Let S be a nonempty subset of Q satisfying:

1. every s in S has exactly k incomparable neighbors already inside Q;
2. the incomparability graph induced on S is connected;
3. every q in Q minus S is incomparable with at least one member of S.

Then the exterior vertices partition into A and B such that

    A<S<B,       A<B,

and every A/Q relation is a<q or incomparability, while every Q/B relation is q<b or incomparability.

**Proof.** Saturation leaves no outside incomparable neighbor for any s in S. An outside vertex z is therefore comparable with all S. Its orientation cannot change across an incomparability edge s||t: s<z<t would imply s<t, and the reversed orientation gives the other contradiction. Connectivity makes the orientation uniform on S. Put z in A if z<S, and in B if S<z. Any s in S gives A<s<B.

If a in A lay above q in Q minus S, choose s in S incomparable with q. Then q<a<s, a contradiction. Thus a<q or a||q. The upper statement is dual. ∎

The saturation hypothesis can be replaced by the explicit assumption that every outside vertex is comparable with every S, as in Theorem B. Connectedness and domination remain sufficient for the same conclusion.

Both structural conditions matter. For two two-element antichains in ordinal sum, all four core vertices are degree-one saturated, but their incomparability graph is disconnected. A new vertex can lie between the components. For a core r<{s,t} with s||t, the connected saturated set {s,t} does not dominate r; a new vertex can lie above r and below s,t. These are controls on the decomposition lemma, not counterexamples to Gold.

For the displayed Q, the incomparability degrees are

    (6,7,3,6,7,7,7,7,6,7,5,6).

Its saturated set is exactly S. The incomparability graph on S is connected, and S dominates Q. All six remaining vertices are boundary ports:

    lower ports: 0,2,3;       upper ports: 8,10,11.

The former are minimal and each lies below some S vertex; the latter are maximal and each lies above some S vertex. Therefore A is below all three upper ports, and B is above all three lower ports.

## 4. Chain compression with a barrier instead of roots

No critical-root argument is needed after Section 3. Fix exterior extension words sigma of A and tau of B. Define

    alpha(q) = number of A letters following q's last A predecessor,
    beta(q)  = number of B letters preceding q's first B successor.

As usual, alpha is antitone along Q, beta is isotone, and both vanish on S. Replace A by a lower chain of length max alpha and B by an upper chain of length max beta, with

    c_i<q iff i<=max(alpha)-alpha(q),
    q<d_j iff j>beta(q),

and put the whole lower chain below the whole upper chain.

For any core extension word w, let pre_S(w) be the word before its first S vertex, and let suf_S(w) be the word after its last S vertex. Define G(t_1,...,t_m) to count weakly decreasing nonnegative integers z_1>=...>=z_m with z_i<=t_i, with G(empty)=1. The conditioned lift count is exactly

    G(alpha(pre_S(w))) G(beta(reverse(suf_S(w)))).           (5)

Indeed A precedes all S and can interleave only with pre_S(w); the last-predecessor constraints give its terminal gaps. The initial A letters before every retained gap have one forced position once sigma is fixed. Delete that initial segment. The dual deletion removes a forced final B segment. Each remaining gap assignment determines exactly one merge. A<B separates the lower and upper choices.

Thus, for every core-word event C simultaneously,

    N_P(C) = sum_(sigma,tau) N_(alpha_sigma,beta_tau)(C).     (6)

In particular, all four cells of (1), their total, and their Gold margins transfer linearly. Normalized mixture coefficients are the atom extension totals divided by E, with multiplicities from exterior-order pairs. They are not uniform over atoms.

When P is k-thin, the usual bounds alpha(q)+beta(q)<=k-d_Q(q) and support-size<=k hold. The replacement chains remain k-thin, by the same degree argument as in the earlier compression theorem. The full-root bound |Q|<=2k is not used or claimed for this generalized barrier setting.

## 5. The 240-atom certificate for Theorem A

For this Q, every allowed conditioned lower profile is exactly

    alpha=(a,0,b,c,0,0,0,0,0,0,0,0),
    a,c in {0,1}, b in {0,1,2,3,4}.

Every allowed conditioned upper profile is exactly

    beta=(0,0,0,0,0,0,0,0,f,0,g,h),
    f,h in {0,1}, g in {0,1,2}.

These are 20 lower and 12 upper profiles. Their supports are disjoint, so all 240 pairs obey the shared budgets. Every pair is an actual seven-thin chain model containing Q, of order 12 through 18. Generic bounded monotone profile enumeration independently confirms this explicit rectangle.

Order the four cells as TT,TF,FT,FF. Direct integer ideal dynamic programming gives the following exact minimum normalized margins over the entire 240-atom family:

    TT: 78113/453252,
    TF: 25363/180473,
    FT: 50695/453252,
    FF: 4768/24065.                                         (7)

Every number is greater than 1/9. Equation (6) therefore proves Theorem A for arbitrary finite exteriors. The receipt records profiles attaining all four minima and the endpoints of (3). In particular the FT minimum is realized by

    alpha=(1,0,4,0,0,0,0,0,0,0,0,0),
    beta =(0,0,0,0,0,0,0,0,1,0,2,0).

This atom has E=453252 and cells (102269,97126,148700,105157). Its FT margin is 50695, so the uniform fixed-tree constant in (2) is sharp.

Nothing is inferred merely from each atom having some Gold strategy. The SAME tree and the SAME four linear inequalities hold on all atoms and hence on their entire cone.

## 6. A general affine-tail and deletion criterion

The remaining argument is reusable for other cores. Consider one lower boundary port u, with gap cap b, while every other lower cap is at most an integer r. For any ordered boundary word w and b>=r,

    G(alpha_b(w)) = G(alpha_r(w))
      + (b-r) * 1_{w begins with u} * G(alpha_r(w without u)). (8)

If u does not appear first, a preceding coordinate bounds its gap by r, so increasing b changes nothing. If u is first, each newly permitted first gap r+1,...,b leaves the same unrestricted bounded suffix gap family. This proves (8) directly, including r=0 and the empty word.

When u is minimal in Q, core extensions beginning with u are in bijection with extensions of Q minus u, by deleting that first letter. The coefficient of b-r in (8) is therefore the lift law of the corresponding actual chain model of Q minus u. The other boundary side is unchanged. The upper version deletes a maximal vertex v when its reversed suffix begins with v, equivalently when v is last in the core word.

**Finite free-port criterion.** In an anchored two-sided core decomposition, suppose u and v are distinct boundary vertices outside the nonempty anchor set, with u minimal and v maximal in Q. Suppose all lower gap caps except that of u are bounded by r, and all upper caps except that of v are bounded by s. For comparisons avoiding u,v, it suffices to check:

- lower free cap 0,...,r, or delete u;
- upper free cap 0,...,s, or delete v;
- every allowed bounded choice of the other profile coordinates.

For free cap b>=r, use the r-cap model with coefficient 1 and the deletion model with coefficient b-r. Do the same on the upper side and multiply the two expansions. All coefficients are nonnegative. Every simultaneous cell-count vector is consequently in the cone of those finitely many base/deletion models.

The bounded profile choices may be over-enumerated, provided each proposed generator is a valid poset and the target inequality is checked on it. No assertion of arbitrary-mixture realizability is needed. This is a sufficient transfer criterion; it does not assert the existence of a common tree for other cores.

The geometric content is exact: an unbounded integer attachment ray contributes a deletion-face generator after a finite breakpoint. It is not an analogy between unrelated dynamical systems. The preserved quantities are all comparison-event counts on retained vertices; the discarded free-port rank is restored by putting a deleted minimal vertex first, or a deleted maximal vertex last, when interpreting the coefficient.

## 7. The 144-generator certificate for Theorem B

Section 3 still gives the same A/B decomposition under Theorem B's explicit outside-comparability assumption. Now

    a,c,f,h in {0,1},       b,g arbitrary nonnegative integers.

There are two free ports, u=2 and v=10. Apply Section 6 with r=s=1. For each of the 16 binary choices (a,c,f,h), there are three lower states (cap 0, cap 1, delete 2) and three upper states (cap 0, cap 1, delete 10). This gives 144 actual generator posets, on 10 through 14 vertices.

The exact normalized margin minima over these generators are

    TT: 771/6245,
    TF: 3231/24464,
    FT: 386/6245,
    FF: 3116/16323.                                         (9)

Each generator count is obtained directly by ideal dynamic programming. As a separate check, every deletion generator is reconstructed by taking the appropriate first or mixed finite difference of full-core models at free cap values 0,1,2. The two computations agree on all four cells.

Every finite free-cap model has a positive base term retaining both free ports, with each free cap 0 or 1. That base belongs to Section 5's bounded family and has all margins strictly greater than (386/6245)E. The other terms have margins at least that multiple by (9). Therefore every finite free-cap model has the strict inequality (4), and summing (6) proves Theorem B for arbitrary exterior posets.

The unbounded theorem uses an exact positive cone, not a numerical extrapolation from large caps. The affine-tail formula is the all-cap proof.

## 8. Sharpness as a genuine unbounded attachment limit

For n>=1 take the actual two-chain model with

    alpha=(1,0,n,0,0,0,0,0,0,0,0,0),
    beta =(0,0,0,0,0,0,0,0,1,0,n,0).

It satisfies Theorem B. It has 12+2n vertices; its core degrees at 2 and 10 grow with n, so this is genuinely outside a fixed-thinness family. The same affine expansion gives

    E(n)    = 12490 n^2 + 78784 n + 122888,
    m_FT(n) =   772 n^2 +  8213 n +  18243.

Consequently

    lim_(n to infinity) m_FT(n)/E(n) = 386/6245.

This proves sharpness of the uniform fixed-tree infimum in Theorem B. Every finite member remains strictly above it. It does not prove optimality among all possible comparison trees.

## 9. Reproduction, scope, and next question

Run in this directory:

    python verify_saturated_barrier_gold.py
    python -O verify_saturated_barrier_gold.py

The standard-library outputs are byte-identical. The verifier reconstructs the induced core, saturated connected dominating set, complete bounded profile universe, and all 240 bounded atoms. Their four cell counts are independently reconstructed from all 61,180 core extensions aggregated into 896 boundary-word/cell fibers. It also reconstructs all 144 base/deletion generators, checks them against full-core finite differences, tests the affine-gap identity at general thresholds, and verifies selected unbounded sharpness models and the two structural hostile examples. The receipt and manifest freeze the exact counts and hashes.

The 240-atom certificate is an all-exterior theorem for one new induced pattern, not a census of seven-thin posets. The general saturated-barrier and free-port lemmas may be used without this particular seed, but still require a common linear comparison certificate. The harder ten-point refinement Q'=(0,0,0,6,0,19,2,82,86,215) is not settled here: its 11,769-atom adaptive-cover investigation remains exploratory, with no claimed Gold obstruction or closure.

The useful next structural target is another saturated dominating pattern whose unrestricted attachment cone has a common strategy, or a core with several unbounded ports where the affine deletion expansion becomes a higher-dimensional positive polynomial decomposition. Neither is inferred from the present two-port result.
