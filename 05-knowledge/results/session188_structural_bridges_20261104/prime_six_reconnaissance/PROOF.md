# Every six-block quotient converse is visible in a fixed ambient

## 1. Finite endpoint certificate

For a tournament Q, let S_i(Q) and E_i(Q) count directed Hamilton paths starting and ending at its physical quotient label i. Put p_i=S_i+E_i modulo two.

**Finite six-name theorem.** Every tournament on six labeled vertices has p nonzero.

This statement is certified exhaustively, rather than extrapolated from a sample. Encode the15 unordered pairs in lexicographic order, with bit one meaning smaller label points to larger. The table covers every integer code0 through32767 exactly once. It retains the complete integer S and E vectors and their parity mask. A forward subset/last-vertex recursion computes E; applying the same recursion to the converse computes S. Initial singleton states have value one. Every Hamilton path has a unique penultimate state, so this recurrence is exact and exhaustive. Counts are at most5! per endpoint and fit safely in the stated integer types.

The complete parity-support histogram is:

-13,920 tournaments have exactly two odd coordinates;
-12,480 have exactly four;
-6,368 have all six.

In particular no code has p=0. The table also records a proper nonsingleton module witness whenever one exists, found by testing every subset. Exactly10,320 codes are prime. Their support histogram is720,5,760,3,840 respectively. Fifteen explicit relabeling orbits partition those prime codes: fourteen of size720 and one of size240. These latter census facts are additional controls; the baseline theorem uses the complete labeled table and does not depend on isomorphism reduction.

`endpoint_table.tsv` has no header. Its15 integer columns are code, first proper-module mask (zero if prime), p-mask, S0,...,S5, E0,...,E5. The rows are ordered by code. The verifier regenerates and compares the entire table directly, not only a digest. An independent audit must validate both completeness and the recurrence before this finite certificate is used.

## 2. Arbitrary nonempty children and fixed physical ambient

Let Q be ANY six-name tournament, with arbitrary nonempty child tournaments A_i on fixed disjoint physical sets. Put T=Q[A_i] and U=Q^op[A_i]. Corresponding target interiors may also be replaced while preserving their complete source ordered-cover profiles. Choose any index with p_i=1, for example the first set bit in the retained table.

The already proved endpoint-parity transfer (`../parity_profile_transfer/PROOF.md`, Sections1 and3) gives a single child-uniform physical query:

    a=all names of A_i,     c=empty,     x->y.

Its full Hamilton response difference is2 modulo4. The proof does not require any child to be nontrivial and does not assume equal child sizes. It uses the exact uniform-block start identity and the fact that modulo two only one child run survives. It is a statement about the complete physical response, including every cubic and quartic contribution.

Substitute T or U at the root of ANY fixed ambient tournament R. Give each port the original root's incidence at every exterior name, on both sides, leaving all exterior arcs fixed. The augmented core is then a common module. Its complete ordered-cover contrast is(2,0,0,...) modulo4, because all same-order higher cover residues are universal modulo4. Substitution multiplies the first coordinate by H(R), which is odd. Thus the same physical query separates in every fixed ambient, still by2 modulo4.

This proves an unbounded-child-order theorem from the finite six-name certificate. It is not unrestricted reconstruction of arbitrary twelve-name or larger tournaments, and no physical child name is collapsed in an actual query.

## 3. Source-relative prime-root families

For the following consequence only, assume Q is prime. Consider simultaneous XOR reversals of ORIGINAL modules contained in the core, with the whole core selected oddly after parity reduction. Every proper original module then lies within a child. The inherited hereditary monotonicity theorem (`../module_decomposition/PROOF.md`) implies that target child cover profiles are coordinatewise at least the source profiles; whole-child converse preserves them.

Use the endpoint query z from Section2 and its core converse-complement Jz=(1-c,1-a). Outside the core both queries keep the same root incidences, not complementary exterior masks. Complete augmented-profile duality (`../monotone_dual_queries/PROOF.md`) gives opposite exact baseline ambient contrasts. Select the member whose source-minus-baseline difference is negative. Increasing target child profiles cannot cancel this strict negative difference, since each child is uniform to both ports and exact run substitution has nonnegative coefficients.

Consequently the fixed two-query battery {z,Jz} detects every such source-relative original-module family. For a fixed source and ambient, one baseline-sign-selected member works for all those targets. No laminarity of the descendant selections is assumed. This is not pairwise injectivity between arbitrary family outputs, nor does it cover operations changing exterior names or the ambient quotient. The prime condition is needed here, not for the baseline six-block theorem.

## 4. Scope boundary

The endpoint theorem is finite in quotient order: it covers every six-name quotient, not every even order. The prime eight-name tournament with outgoing rows

    (238,4,240,118,195,146,34,74)

has p=0. The verifier retains its full integer endpoint vectors and a direct module test. This does not assert that its complete response is blind; other methods already distinguish it. It only prevents an unjustified general-even-order extension of the six-name endpoint certificate.

The existing highest-changed-node restoration theorem can use this new finite result after independent review. This packet itself does not claim to solve arbitrary simultaneous reversal noncancellation.

## 5. Executable controls and full named data

In addition to the complete labeled table, the verifier uses exact physical counts on singleton, all-pair and unequal mixed children, with arbitrary-fixed-ambient controls and profile-preserving child orientation changes. It retains an entire4096-query physical response bank for a six-name prime example, both orientations, together with all484 Boolean coefficient coordinates of bidegree at most(2,2). Thus cubic and quartic terms are retained explicitly rather than replaced by endpoint data in the physical controls.

There are also source-relative root-family controls, including a crossing pair of original intervals inside a transitive child. Those finite examples support the implementation; the family proof is the inherited monotonicity/duality argument above.

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. Both modes regenerate the complete table and all physical data. No package installation is needed. The exact frozen publication allowlist is in README.md.
