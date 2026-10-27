# Session75 Phase2 independent audit: black cut counterfamilies

**Verdict: ACCEPTED in the scopes below.** Audit date 2026-10-27.
No shared repository edits or remote actions were performed by this auditor.
The source research pin is `96e35642110104cc547516a63f37eebcf930f707`.
The final producer milestone's manifest SHA256 is
`5bc89f4d9673325311088321d984811a2ddbf123be6739e81ff75b314f5d28a3`.
All 43 entries of that manifest were checked against the frozen source bytes.

## 1. Accepted conclusions and limits

1. **All-order refutation of A.** For every order n >= 10, the specified
   even/odd rooted-converse family has lambda = 1, minimum semidegree delta = 2,
   and black distance to nonstrong tournaments exactly 2. The initial source
   is strong and black-active. Each minimum-cut decreasing cone has exactly
   nine Hamilton words, all blue in the full anti-automorphism sense.
2. **Minimum order 10 is established**, not merely the smallest constructed
   example. The complete ordinary-class universes through order 9 are
   independently authenticated, and every strong black-active source has
   literal geodesics attaining both lambda and delta. In particular all 353
   strong order-7 classes are included, not only regular classes.
3. **Separate survivor supplement accepted.** Replacing the original second
   path by the supplement's Q_core gives a two-black route to one named source
   plus a genuinely strong deleted core at every member of both families.
   Thus both target distances equal 2 on these families. This is a separate
   proof and route; it is NOT a property of the original nonstrong route.
4. **Unique-unbalanced-minimum-shore repair accepted.** A strong tournament
   with a unique minimum outgoing shore of size different from n/2 has black
   nonstrong distance lambda.
5. **Global B remains OPEN.** Neither the finite classes nor these infinite
   families prove the core-distance equality for every strong black-active
   tournament. No all-order black connectivity statement is accepted here.

The order-4 fixed-shore obstruction is a hostile to an arbitrary-shore method,
not a global counterexample to A. The audit explicitly finds a black exit at
another shore there. At order 10 there are 99 blue Hamilton words in total;
only nine of these lie in the minimum-cut cone. Those counts must not be
conflated.

## 2. Independent structural audit

### Carrier, current, and literal gauge

For a directed Hamilton word P, reversing its path arcs changes each cut by
`1_A(last P) - 1_A(first P)`: every transition into/out of A cancels except the
endpoint imbalance. In particular, only the first and last named degrees
change, by -1 and +1. Therefore lambda is a lower bound for nonstrongness,
and delta is a lower bound for reaching a source or sink, hence for the
singleton-plus-strong-core target.

The native nonpath-complement operation K satisfies K_P(T) = op(R_P(T)).
For an R-route T_j and words P_j, retain the SAME named initial source and use
native states op^j(T_j), with word P_j reversed exactly when j is odd.
This gives the literal K-route, not a converse-merged or score-state route.
Both target predicates and blackness transport. Our checker additionally
constructs every K result by complementing nonpath unordered pairs and checks
actual anti equations at both endpoints. If the producer constructed a finite
witness on the source's converse, this audit FIRST reverses every word to
restore the original named source before applying the alternating gauge.

### Rooted orbit and the complete exit cones

Deleting the root r of U = C3[T_a,C3,{r}] leaves T_a -> C3. Thus every
root-ending Hamilton word is the increasing A word, one of three cyclic B
words, then r. Rotations of B fix the rest and transitively permute exactly
these three words. This is a rooted word-orbit statement, not purity of U.

The even join has one U-to-V bridge. Any Hamilton word beginning in U and
ending in V must use exactly that bridge, no reverse crossing, and exhaust
U before V. There are exactly 3 times 3 such words. The reversed V block is
another root-ending U word; the rooted automorphism between those words,
and its inverse on the opposite block, give the full positional
anti-automorphism. Cross arcs are uniform except the paired roots.

Every nonempty root-avoiding subset of U has at least two outgoing internal
arcs: two B vertices already beat r; one B vertex beats r and its absent
successor; no B vertices leaves every A vertex beating all three B vertices.
A proper root-containing subset has an internal outgoing arc plus the bridge.
Two proper nonempty block intersections contribute an internal outgoing arc
each. One-block V shores have ample cross arcs. Complement after the
block-exchanging anti-automorphism covers the remaining whole-block cases.
Hence the even minimum shore is uniquely U.

In the odd case all direct V-to-U arcs remain; the bridge is U-root -> w ->
V-root. The same cut partition excludes every shore except U and U+w.
Both decreasing cones are forced to be alpha_U, w, beta_V, so coincide with
nine words. Reversal fixes w and satisfies its exceptional-root equations.
Any one-step nonstrong target would require a decreasing initial unit cut;
all its possible words have now been exhausted and are blue. Thus no black
one-step exit exists in either parity.

### Upper-bound routes and terminal types

The stated first word is a Hamilton cycle cut inside V. Its endpoint degrees
exclude an anti-automorphism, so it is black. Both endpoints lie outside all
initial minimum shores, leaving those cut values at 1. Every other cut was
at least 2 and changes by at most 1; the intermediate is strong.

The even Q is directed using the reversed internal alpha paths and newly
reversed cross edge. If its positional reversal were an anti-automorphism,
the two full induced blocks would be converse-isomorphic. They are not:
the discrepancy is one arc flipped between equal-degree B vertices, changing
the sum of internal squared degrees by exactly 2. The odd Q instead has
endpoint degree sum n-4, not n-1, which already excludes blue. Both reach a
zero U cut. These are nonstrong targets, with no core claim attached.

For the SEPARATE Q_core, both words end at the original indegree-2 vertex
v = a+2. The endpoint current makes it a source. The second word is directed
and black by the verified endpoint-degree obstruction. The two internal U
reversals cancel, restoring strong U. In V minus v, the triangle
r -> a -> a+1 -> r is strong, both B vertices beat every A vertex, and every
A vertex beats r. Thus that entire block is strong. In even order the one
opposing cross arc U-root -> (a+1)_V joins the strong blocks. In odd order
U-root -> w -> (a+1)_V gives a directed cycle of the blocks and w. Therefore
the ACTUAL deleted core is strong. This verifies the arbitrary-a proof,
not just the finite controls.

### Unique unbalanced shore repair

Choose a minimum-shore-decreasing Hamilton word at every strong step (cut a
Hamilton cycle at an arc entering the shore). After j steps its cut is
lambda-j, while every other cut is at least lambda+1-j. Hence it remains
uniquely minimum. An anti-automorphism would take that shore to the complement
of its image, another minimum outgoing shore. Uniqueness would force its
size to be n/2. The unbalanced hypothesis therefore forbids ALL
anti-automorphisms at every intermediate source; all selected paths are
black. The lower bound and constructed lambda-step route agree.

## 3. Independent finite evidence

All code under this audit root is independently authored, with no producer
module imports. Shared imports are only between this audit's own scripts.
Family construction uses a pairwise adjacency-matrix predicate, rather than
the producer's bit-flip constructor. The final source/word/native endpoints
are then explicitly compared with the frozen producer receipts.

`audit_families.py` checks a = 1..25 and a = 100, both parities. It verifies
rooted rotations and root-ending path completeness through a=8, exact minimum
shores, every listed full anti equation, both literal operations/colors,
strong intermediate and nonstrong target. Complete Hamilton-word and direct
cut exhaustions independently give:

| n | Hamilton words | Black words | Proper cuts | Minimum shores |
|---|---:|---:|---:|---|
| 10 | 657 | 558 | 1,022 | 992 |
| 11 | 1,215 | 1,116 | 2,046 | 992, 2016 |
| 12 | 2,601 | 2,394 | 4,094 | 4032 |
| 13 | 4,995 | 4,788 | 8,190 | 4032, 8128 |

Every one of these black one-step neighbors is strong. At larger controls,
minimum shores are computed exactly, not sampled: for each size k,
`min c(A) = sum of k smallest outdegrees - k(k-1)/2`. All cutoff ties are
enumerated, so the complete minimizing family is retained.

`audit_survivors.py` independently checks the separate core route on the
same 52 controls and verifies every structural core arc described above.
It binds all 40 producer survivor cases to literal named K endpoints. It
also replays the parent's distinct n10 route and n11 route. The parent's n11
labels differ: the family-to-parent permutation is
`(0,1,2,3,4,6,7,8,9,10,5)`. This permutation is explicitly verified on all arcs;
the parent's source is not silently identified with the family source.

### Ordinary source universe and provenance

`SOURCE_TABLE_BINDINGS.json` records exact byte matches, SHA256s, and git paths
at pin 96e3564 for all six source tables. No inherited class cardinality is
used alone as a completeness argument.

`audit_universe.cpp` is an independent canonicalizer: start with degrees,
iteratively refine using outgoing counts to color cells, then enumerate
within the resulting cells and take the MAXIMAL arc word. Canonical-key
multiplicity gives the automorphism order. It checks the table's order,
rejects isomorphic duplicates, and sums n!/Aut. Refinement colors are intrinsic,
so every automorphism preserves the cells and every isomorphism is covered.
The disjoint orbit masses equal the entire labelled tournament universe:

| n | Ordinary classes | Labelled orbit mass | Strong classes |
|---|---:|---:|---:|
| 4 | 4 | 64 | 1 |
| 5 | 12 | 1,024 | 6 |
| 6 | 56 | 32,768 | 35 |
| 7 | 456 | 2,097,152 | 353 |
| 8 | 6,880 | 268,435,456 | 6,008 |
| 9 | 191,536 | 68,719,476,736 | 178,133 |

This differs from the producer's raw score-block canonicalization; at n9
it checks 5,737,246 refined labelings versus the producer's 27,282,165.
These are ordinary classes, never converse-merged classes.

`audit_finite.py` binds every compact witness to its exact source index/code,
checks strongness, recalculates lambda=delta, and replays every eligible named
core route of that length. It checks 296,985 R edges and 296,985 literal K
edges. At n5 the one strong pure-blue exclusion is exhaustively verified
(15 words); C3 is independently verified as the order-3 boundary (three words).
Thus all 184,535 strong black-active ordinary classes in orders4..9 are
covered. The original nonstrong and core lower bounds prove exact distances,
which establishes the global A counterexample's minimum order10.

## 4. Reproduction, acceptance, and stopping point

Run `bash run_verification.sh` in this directory. It checks the frozen local
inventory; compiles fresh C++ normal/optimized configurations; reruns all
independent checks using normal Python and Python -O; and compares each
result with the paired receipts. Every essential reject is an explicit error,
not an optimization-sensitive assert. `audit_controls.py` confirms rejection
of malformed tournaments, non-Hamilton words, reversed nondirected words,
a blue move, and a nonstrong target wrongly presented as a strong core.

Paired family, survivor, finite, universe, and hostile-control receipts are
byte-identical. No mathematical blocker remains in the accepted scopes.
The infinite proofs, finite minimum-order gate, and global-B open question
are separate. Reserved theorem namespaces were not used as dependencies or
promoted by this audit. This is the coherent verified milestone and stopping
point for the requested independent Phase2 audit.
