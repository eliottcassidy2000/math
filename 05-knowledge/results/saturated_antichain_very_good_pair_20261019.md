# One-slack antichains force ambient very good pairs

2026-10-19. PROVED structural restriction; root-reviewed and frozen. An all-size structural sufficient condition, not a full morphology exclusion or a solution of seven-thin balance/Gold. No historical novelty claim is made.

## 1. Ambient theorem

Let P be a finite poset, S an antichain of m >= 3 vertices, and C a disjoint set of t vertices each incomparable with every vertex of S. Suppose each s in S has incomparability degree at most m+t. Then S contains a subset of size at least ceil(m/2) whose every pair is very good. In particular P has a balanced pair and weak Gold, using the good-pair consequences cited below.

Here a very good pair is an incomparable pair which, in P or its dual, has equal strict downsets and both directed upset differences are chains. Singletons and empty sets count as chains. The conclusion is about the full ambient P, not a probability computed on an induced subposet.

### Proof

The m-1 other members of S and the t members of C consume m+t-1 incomparable neighbors of each s. Consequently each s has at most one incomparable neighbor outside S union C.

Partition the remaining vertices into D, those below at least one member of S; U, those above at least one member of S; and W, those incomparable with every member of S. No vertex belongs to both D and U: s < z < s' would contradict the antichain condition. Every z in D is either below or incomparable with each s; dually every z in U is above or incomparable with each s.

If W is nonempty, the one-slack bound forces W to have exactly one vertex and forces every member of D below every s and every member of U above every s. All members of S then have identical downsets and upsets, so every pair is an exact twin pair.

Otherwise W is empty. Write D_s = D minus the strict downset of s, and U_s = U minus the strict upset of s. These are precisely its remaining incomparable neighbors, so |D_s|+|U_s| <= 1. Let d be the number of s with D_s nonempty and u the number with U_s nonempty. Their supports are disjoint, hence d+u <= m.

The m-d vertices with D_s empty have identical downsets D. Each has upset U with at most one point removed, so every pair of their upsets has directed differences of size at most one. They form a very-good-pair clique. Dually the m-u vertices with U_s empty form such a clique. The larger clique has size at least ceil(m/2), which is at least two. This proves the theorem. No order or cardinality assumptions on D, U or C were used.

## 2. Seven-thin consequences

In a seven-thin poset, the theorem applies whenever an antichain S with m >= 3 has t common incomparable vertices and m+t >= 7. (If m+t > 8 the degree bound already makes the hypotheses impossible.) Therefore any hypothetical seven-thin counterexample to either ordinary balance or weak Gold satisfies all of the following necessary restrictions:

- Its width is at most six: any antichain of seven supplies m=7,t=0.
- No six-element antichain has even one common incomparable vertex.
- No five-element antichain has two common incomparable vertices.
- No four-element antichain has three common incomparable vertices.
- No three-element antichain has four common incomparable vertices.

Equivalently, the common incomparability neighborhood of any antichain S of size 3 through 6 has size at most 6-|S|. In particular, the induced poset on I(v) has width at most five for every vertex v.

These statements are necessary conditions for counterexamples, not assertions that satisfying them yields a counterexample. They do not require a critical root, a majority order, or a canonical choice.

## 3. Application to large-exclusive root cores

Use the complete critical-pair core notation Q={x,y} union L union T union Z, where L=D(y) minus D(x), T=U(x) minus U(y), and Z=I(x) intersect I(y). Every member of L is incomparable with x, and every member of T is incomparable with y.

If |L|=6 and L is an antichain, take S=L and C={x}. The theorem supplies at least three members of L whose every pair is an ambient very good pair. The dual statement holds for an antichain T of size six. Seven-thinness also forces Z empty in these cases, since d_inc(x)=1+|L|+|Z| <= 7.

Thus each remaining c=0 branch (l,u) with u=6, and also l=6 when present, loses all cores whose six-element exclusive side is an antichain. In particular a residual (6,6,0) core must have a comparability on each of its two exclusive sides. The restriction covers all ambient attachments and arbitrary orders on the opposite side. It is not limited to height-two posets: exterior D and U can have arbitrarily large height.

This is a proper subfamily restriction, not a full exclusion of (2,6,0), (3,6,0), (4,6,0), (5,6,0), or (6,6,0). It neither assumes nor duplicates the pending (2,2,4), (3,3,1), or canonical (2,4,1) scans. More generally any six-element antichain inside a root incomparability neighborhood is forbidden, even when it is not the whole exclusive side.

## 4. Why an ambient proof is needed

Adding exterior vertices can change the induced-core extension distribution. An exact control below uses a c=0 (6,6) core with antichains L and T and every relation L<T, L<y, x<T. Adding a vertex a below x,y,T and all but one member of L preserves seven-thinness and the complete root core, but changes the probability of the exceptional L vertex preceding another L vertex from one half to 193/343. The core has 25,920,000 extensions, with 12,960,000 placing the selected exceptional vertex first; the ambient poset has 29,635,200 extensions, with 16,675,200 placing it first. Therefore core balance is not being transferred by an unjustified uniform-weight argument. The theorem instead constructs the very good pair directly in P.

## 5. Dependencies and controls

The primary source was checked directly on 2026-10-19: Definition 5 on page 3 gives the equal-downset and two chain-difference conditions, while Definition 1 and Theorem 2 on page 2 give good pairs and the balanced-pair consequence. The definition and balanced-pair consequence of very good pairs are attributed to Zaguia, Definitions 1,5 and Theorem 2, author-deposited paper https://arxiv.org/pdf/1610.00809v3. The weak Gold consequence uses the already established [good-pair Gold synthesis](good_pair_gold_synthesis_20261019.md), in its convention allowing determined or repeated comparisons. The height-two result [packet112](seven_thin_height_two_very_good_gold_20261019.md) is contextual only and is not needed for this proof. The argument does not use canonical packet126 or any morphology exclusion as a premise.

The accompanying standard-library control checks all permutation-orbit representatives of assignments of no defect, one lower defect, or one upper defect to antichains of sizes 3 through 8, using three labeled possible defects on each side. It verifies the promised clique bound by explicit set differences. These are abstract defect configurations, not a census of posets; realizability is unnecessary for this conservative check. A separate exact ideal-DP calculation verifies the ambient exterior-weight example and a higher-height fixture. Guards remain enabled under optimization. Normal and optimized runs agree on all 86,186 guards. There are 6,399 permutation-orbit representatives, representing 6,725,544 labeled assignments across these six antichain sizes. The finite checks are regression controls; the elementary proof above establishes the all-size theorem.
