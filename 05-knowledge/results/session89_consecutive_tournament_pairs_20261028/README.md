# Consecutive tournament pairs: coordinates, symmetry, and reusable computation

Session 89, 2026-10-28. The starting observation is exact:

    n² = binom(n,2) + binom(n+1,2).

This session turns that bit correspondence into explicit maps, distinguishes the different isomorphism actions, implements paired enumeration, and connects the construction to self-converse tournaments and a smallest closed arithmetic-decision carrier. Classical mechanisms and prior repository results are credited; external novelty is not asserted merely because a construction was missing from an earlier packet.

## 1. The labelled bijection works; the quotient must be specified

Split an n×n bit matrix into its upper triangle, lower triangle and diagonal. One triangle orients T_n; the other triangle and the n diagonal bits orient a second tournament U_(n+1), with its new vertex distinguished. The inverse recovers every bit.

There are three different uses of “consecutive pair”:

1. **Independent shared-label pair:** T_n and rooted U_(n+1) are independently chosen, but the first n labels are shared.
2. **Nested extension:** U restricted to those n vertices must equal T. This has only one independent tournament layer.
3. **Overlapping stream:** a middle tournament occurs in two consecutive pairs. Its canonical representative, group and pointed data can be reused, but the two relative alignments must still be joined correctly.

Ordinary matrix conjugation permutes rows and columns together. Natural pair isomorphism relabels the two tournaments together. These actions are not the same under the fixed coordinate split:

| n | Natural aligned-pair classes | Ordinary matrix-conjugacy classes |
|---:|---:|---:|
|0|1|1|
|1|2|2|
|2|8|10|
|3|88|104|
|4|2,752|3,044|

The difference cannot be removed by finding a more ingenious S_n-equivariant bijection. The zero matrix is fixed by every permutation, whereas an even-order subgroup cannot fix a tournament. A movable orientation frame repairs a larger, explicitly anchored object; it does not make the unanchored quotient disappear.

Read [the codec, affine action and transpose law](encoding/ENCODING_AND_ACTIONS.md).

## 2. Exactly when an orientation anchor exists

For a faithful permutation subgroup H≤S_n, an H-equivariant matrix/pair codec exists **if and only if |H| is odd**. An invariant tournament supplies the reference orientation: on each anchor arc i→j, read the two directed matrix entries relative to that arc. The obstruction is the symmetric matrix pairs00/11, which need a tie-breaking orientation.

The anchor can be found from supplied generators by parity propagation on unordered edges. A compact solution takes O((s+1)n²) elementary visits for s generators; enumerating all anchors still costs its exponential output size. A contradiction supplies a generator word reversing an edge. Checking only the generators' orders fails: two3-cycles can generate A4, which has no invariant tournament.

A bare named vertex leaves S_(n−1), so it is insufficient for n≥3. A class plus a vertex in an **already oriented tournament** is more useful, but it leaves that vertex's actual stabilizer. For P7, one root still leaves a group of order3. Treating it as a unique labelling would erase alignment data.

Containment H≤Aut(A) is enough for the codec; equality is not guaranteed. The regular C3×C3 action on nine points has16 invariant anchors, all with full automorphism group of order81. Thus the declared subgroup remains part of the cache contract.

Read [the anchor theorem and certificates](anchors/anchor_theorem.md).

## 3. Enumerate shared alignments once, then reuse residual mask orbits

Choose representatives of two n-tournaments with groups G and H. Their relative core alignments are double cosets

    G \ S_n / H.

After choosing an alignment g, the new-root mask is classified under the **concrete** residual group

    K = G ∩ gHg^−1.

The implementation caches complete mask and stabilizer recipes by K, not just by its order. It then returns the full paired representative, named rooted-child class, and complete pair stabilizer.

The complete paired counts through n=5 are

    1, 2, 8, 88, 2752, 279968.

At n=5, 8,784 core alignments require only five distinct residual-group recipes. Seven fresh-setup-inclusive benchmark trials gave median times:

| Method | Median seconds |
|---|---:|
| Joint-candidate reference |0.546617|
| Shared-core reuse |0.263014|
| Rooted-child-first double cosets |1.644683|

The tested shared-core route is2.08× faster; the rooted-child-first route is exact but slower. All three agree on every output, not only the total. The n=5 stabilizer histogram is279,448 trivial,512 of order3, and8 of order5. Independent literal enumeration of all2^(n²) objects covers n≤4; n=5 is a complete three-method bank comparison, not a second33-million-object labelled replay.

Nested extensions are extracted separately, giving1,2,4,12,48,296. For an overlapping stream, after putting the middle object M into a common frame, compatible joins use

    A \ Aut(M) / B,

with the corresponding intersection stabilizers and actual transports. Multiplying class counts or forgetting the middle frame does not implement this join.

Read [the generation proof and output contract](generation/ALIGNED_GENERATION.md) and [reproduction guide](generation/README.md).

## 4. Share scalar-count work, not the already-summed weights

For r≥1 independent tournament layers and s binary vertex attributes, the inherited channel is

    F_n(r,s) = sum_(odd partitions λ of n) 2^(r e(λ)+s length(λ))/z_λ.

A000568 is F_n(1,0); independent aligned pairs use F_n(2,1); nested extensions use F_n(1,1). Ordinary matrix conjugacy includes even partitions too. The r=0 case cannot silently retain the tournament odd-cycle filter.

The fused implementation computes A000568 and paired counts in one response-state traversal, with separate exact weight accumulators. Against the strongest repeated Session87-based computations, seven-run median CPU gains were1.41×,1.28×,1.19× at bounds80,120,160. Peak Python memory rose1.64×,1.76×,1.80× respectively. No smaller state quotient or new count formula is claimed.

Squaring the first accumulated channel is wrong. Cores(3,7) and(5,5) reach the same terminal response state but have different exponents and denominators; their exact weighted moment defect is1536/7. Shared traversal preserves the missing second-moment information by carrying it independently.

Read [the fused counter and measured tradeoffs](counting/README.md).

There is also an exact all-order low-bit consequence. For odd n,

    v2(F_n(r,s)) = r(n−1)/2+s.

For even n≥4, F_n(r,s) is divisible by2^(rn/2+2s), and the next normalized bit is φ(n)/2 modulo2. Thus the bound is exact precisely when φ(n)/2 is odd. The n=0 and n=2 exceptions are handled explicitly. The proof isolates the lowest-valuation cycle types, rather than extrapolating from a table. For ordinary tournaments this determines every odd-order valuation; for paired counts it gives v2=n at odd n. These are cheap checks on computed counts, not certificates of every output bit. See [the proof and independent checks](counting/LOW_BITS.md).

## 5. A structural self-converse lift, not merely equal counts

There is an explicit isomorphism-class bijection between aligned consecutive pairs at n and self-converse tournaments on2n+1 vertices.

Build vertices0,p_i,q_i. A designated anti-involution fixes0 and swaps p_i,q_i. Orient p_i→q_i to normalize each pair; one core tournament determines the positive layer, its converse determines the negative layer, and the discrepancy graph and root bits determine the remaining arcs.

Every odd self-converse tournament admits such an anti-involution. All its anti-involutions are conjugate under its automorphism group, so forgetting the chosen anti-involution causes no additional identification of pair classes. This closes the precise structural gap left open by the earlier Session22 numerical identity. The singleton is handled separately.

But symmetry groups are not preserved wholesale:

    Aut(pair) ≅ C_Aut(W)(τ).

The number of anti-involutions is the index of that centralizer. Already the two n=1 pairs have trivial groups, while their lifts are the transitive triangle and C3, with automorphism orders1 and3. A class bijection is weaker than a groupoid equivalence.

Complete independent tournament-bank checks through order seven give self-converse image counts1,2,8,88. The map gives a concrete reduction to odd self-converse tournament isomorphism, not a polynomial-time isomorphism algorithm; a simpler ordinal-sum reduction is also recorded to avoid claiming a new complexity equivalence.

Read [the lift, decoder and centralizer theorem](encoding/ODD_SELF_CONVERSE_LIFT.md) and [the exact prior-status comparison](history/SELF_CONVERSE_PRIOR_STATUS.md).

## 6. A smallest closed carrier for a declared matrix/arithmetic decision

The concurrent nine-prime bank forms this3×3 array:

    19  59   2
     7  83  41
    89 139  13

Its square-σ subsets are exactly the32 binary cut matrices

    C_3 = { M : M_ij=a_i+b_j over F2 }.

Four anchored rectangle parities decide membership. Ordinary independent row and column permutations preserve the code. Natural paired-tournament relabelings do not: one transposition sends the zero matrix to the prime subset{59,7}, whose σ is480, not a square. Only the identity among the six natural S3 actions preserves the whole code.

The follow-through gives an all-size, operation-relative minimality theorem. Let C_n be the same rectangle code, allow every **named cell toggle** and natural paired relabeling, and observe membership in C_n. For n≥3 the exact future-decision equivalence is

    M ~ M' iff M−M' belongs to D_n,
    D_n = { M : M_ij=a_i+a_j+c }.

Thus the smallest carrier has n²−n bits. A canonical representative clears the first row and stores the remaining(n−1)×n rectangle. The quotient has affine update rules; D_n is invariant under the linear parts, not necessarily under the full affine maps. Exceptions: n=1 needs zero bits, n=2 needs one.

At n=3 the current yes/no decision is one output bit. Supporting every future toggle needs the four-bit rectangle syndrome; adding natural relabeling needs exactly six bits. All512 matrices fall into64 future profiles, eight matrices each. The remaining2×3 carrier has a constant-row decision kernel; as three F2² columns this is the same diagonal-kernel query model recovered in Session87. The original prime-coordinate lifts are not identified with the F64 field, its multiplication, or its Frobenius action.

This is a proved smallest closed **decision** carrier for the stated menu. It neither identifies tournament classes nor settles the unrestricted attachment-carrier problem. For n>3 the theorem is about the combinatorial cut code; no unproved general prime-bank realization is assumed.

Read [the bank and action obstruction](history/concurrent/BANK_MATRIX_BOUNDARY.md) and [the all-size carrier proof and APIs](history/concurrent/FUTURE_BANK_CARRIER.md).

## 7. Historical repairs and honest open-problem connections

The March square-codec idea was genuine, but two scripts implemented a lossy upper coordinate. The old map forgot m−1 bits, so its image had dimension m²−m+1. Its own side-three display filled seven of nine cells. We corrected the narrow coordinate expressions in both live scripts and retained warnings that their old saved rank/determinant data were not full-matrix censuses. Historical output was not rewritten, and unrelated explorations were not revalidated.

The [history and complexity map](history/HISTORY_AND_COMPLEXITY_BRIDGES.md) recovers seven earlier mechanisms and separates several legitimate targets:

- Ordinary binary matrix conjugacy carries graph isomorphism; the induced action on the pair codec is not natural paired-tournament isomorphism.
- Simultaneous GL intertwiners are a relaxation of permutation intertwiners. An invertible linear map need not contain a valid relabeling; the structured module-isomorphism problem must not be confused with arbitrary matrix-space nonsingularity/PIT.
- Dyadic monomial substitution can preserve formal coefficients while making degree exponential. It is not a general PIT derandomization.
- Decoding a zero-one permanent preserves that labelled counting problem, not tournament Hamilton-path counts or natural pair class functions.
- General tournament reconstruction from vertex-deleted cards has counterexamples; it cannot be used as an unproved shortcut to the separate simple-graph reconstruction conjecture.

These are explicit reductions, obstructions and testable next targets, not claimed solutions of major open problems. Classical formulas and literature are cited in the linked notes.

Read [the repaired all-size square codec](history/REPAIRED_SQUARE_CODEC.md) and [live-code regression instructions](history/repairs/README.md).

## Verification and use

New mathematical checkers use active exception checks and are replayed under normal and optimized Python. They are standard-library code except the AST-isolated regression of the two existing legacy NumPy codecs, which requires NumPy. No raw result receipt is needed as a committed input. Historical baselines are loaded from existing repository files with exact hashes, not vendored copies.

The proof, producer finite universe and independent-auditor universe are stated separately. Timings are matched local observations with setup included where specified; they do not imply an asymptotic improvement. Final publication checks verify exact file contents and preserve all unrelated repository work.

The common lesson is constructive: choose a canonical object, retain the residual action and transports, and reuse exactly the data the next pair requires. A coordinate bijection, a class bijection, a symmetry-group equivalence, and a future-decision quotient are four distinct achievements.
