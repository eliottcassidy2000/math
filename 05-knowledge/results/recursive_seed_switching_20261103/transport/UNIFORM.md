# A uniform native implementation of the inherited shared-clock +4 family

**PROVED:** an explicit native algorithm for every integer clock t≥1.
**FINITE-EXACT:** saved guards and independent counts at t=1,2,3,4,8.
This is one uniform map between two families, not a returning family, a +2
map, a two-digit induction, or unrestricted Hamilton-path spectrum coverage.

## Inherited family and the remaining implementation step

Session182, **`session182_physical_returns_20261104/induction/PROOF.md`,
sections2 and5**, proves the following counts, terminal-port strongness
criterion and all-strong completed adjacency schedule. Let A=TT_t, t≥1,
and let M(A;B,W) have A→B, W→r, and r→A∪(B\W). Start with B=TT3 on
0→1→2 and W={0,1,2}. The target is B'=TT5 on 0→1→2→3→4 and
W'={0,2,4}. Then

    H(source) = 7·2^t−5,
    H(target) = 21·2^t−11 = 3H(source)+4.

Both full tournaments are strong. The actual root defects E_r(2)−H are
6 and12. These are inherited statements, not new numerical realizations.
Session182 explicitly leaves native Hamilton-complement guards uncertified
for this family. The increment here is a uniform algorithm constructing them.

## The native program for every clock

Retain every old physical name, write the root as39, and use any distinct
prefix names disjoint from {0,1,2,3,4,39}. Perform these steps:

1. Birth4 as a terminal suffix vertex and an in-port. It beats exactly39
   among the old vertices; all other old vertices beat4.
2. Birth3 immediately before4 in the suffix, as a nonport. It beats exactly4
   among the existing vertices; all other vertices, including39, beat3.
3. Reverse the single old pair 1→39, retaining every other old incidence.

After each of the first two births the old induced graph is unchanged. The
suffix orders and port sets are respectively

    (0,1,2; {0,1,2}),
    (0,1,2,4; {0,1,2,4}),
    (0,1,2,3,4; {0,1,2,4}),
    (0,1,2,3,4; {0,2,4}).

The terminal suffix vertex is a port at each stage. The inherited exact
terminal-SCC criterion therefore proves full strongness at every completed
adjacency stage. The intermediate birth-stage counts and root defects are
15·2^t−13 with defect14, and23·2^t−19 with defect20. They also follow
directly from the inherited transitive-port recurrence. No H or strongness
invariant is asserted inside the following temporary-helper gadget.

For step3 use the direct pair-block algorithm in **`transport/FOLLOWUP.md`**.
On the post-birth tournament delete1 and39 temporarily for the purpose of
constructing a word, and build a Hamilton word Q on the remaining vertices
by first-out-neighbor insertion. Before the first Q vertex beaten by either
endpoint, insert the block(1,39) or(39,1), choosing its final endpoint to
beat that Q vertex. Every earlier Q vertex dominates both endpoints. If
there is no such Q vertex, append either block. Thus the full word is legal
in the current tournament or its single-pair reversal, and the selected
pair is consecutive.

Apply the inherited one-helper arc-flip gadget from **METAGRAPH section5,
`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`** in that orientation.
If necessary traverse its two involutions in reverse order. It births one
helper, performs two legally guarded off-path complements, and deletes only
that helper. The exact endpoint reverses1→39 and changes no other pair.
This constructs the required guards for every t≥1 in O((t+6)²) adjacency
tests. It uses exactly two permanent births, one temporary helper birth and
deletion, and two native complements, with peak order t+7. All t+4 source
names survive. The native gadget and family are credited to their original
proofs; the contribution is their explicit uniform composition.

In this family the construction simplifies to two explicit words. Let h be
the helper, beating exactly1, and let A denote the ordered prefix word:

    p = (1,39,A,0,2,3,4,h),
    q = (1,h,39,A,0,2,3,4).

The word p is legal before the first complement:1→39→A→0→2→3→4→h.
After that complement the previously unprotected edges h→1 and39→h have
reversed, so1→h→39; all remaining edges of q were protected by p. Hence q
is legal too. The protected-edge sets differ exactly on
{1–39,4–h,1–h,h–39}. Two off-path complements reverse precisely their
symmetric difference. Deleting h therefore leaves just the required reversal
1→39. This is a direct all-clock guard proof, with O(t) word-output length;
the generic direct compiler reproduces these same words in the finite controls.

## Exact controls and transported named data

`uniform.py` constructs all four completed stages directly from named
adjacency at clocks1,2,3,4,8. It imports the direct compiler only for its two
guard words. The separate elementary interpreter from
`zero_ready_component_returns_20261103/component/audit_native.py` reconstructs
the helper birth and both literal off-path complements independently, then
checks that helper deletion gives exactly the target matrix.

For all20 completed graphs, an exact small-feedback path-bin count and a
separate endpoint subset DP agree with the four displayed H formulas.
Reachability independently checks strongness. Each graph retains the root's
eight path-cover/endpoint fields, its actual defect and current golden,
F49 andF343 observers. Source and final target also retain the ordered
pair(39,1), both five-field endpoint responses, P2 and ordered W2. Name1
persists physically while its root-port membership changes; these current
responses are recomputed, not transported by a field automorphism. No full
joint17 observer or arbitrary-operation closure is asserted for this family.

The first four numerical pairs are the inherited9→31,23→73,51→157 and
107→325; clock8 gives1787→5365. Both legal-guard positive controls and
illegal/repeated-name hostile controls are checked with explicit exceptions.

Reproduce the stdout receipts with:

    python3 transport/uniform.py > transport/uniform.json
    python3 -O transport/uniform.py > transport/uniform.optimized.json

The normal and optimized receipts must be byte-identical. The finite checks
support the implementation; the preceding argument supplies its all-clock
quantifier. Iterating from the TT5 output under both digits still requires
new transformations and a returning representation class.
