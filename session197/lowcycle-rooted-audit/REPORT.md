# Independent rooted low-cycle audit

2026-11-05. **PASS for the all-order laws, start-one grammar, bridge
equivalence, bounded image, and feedback-vertex family, after one local
statement correction.**

Audited ROOTED.md hash:
`ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`.
Verifier hash:
`7c2ae10b5e77e07e607e040fa438dba6ffe78b92a61f1623f444815d69a9c303`.
Compiler/search hash:
`f770a8e7650ae1102c877dbb7d4d52f4533317cad60cf911102c1d13e08e1bb1`.

## Proof findings

The one-seam start count is `d*s`: a path beginning at the new root must
immediately use the sole reverse cross arc, then exhaust the old piece,
then traverse the gadget minus its root. The cycle count is `s*E`, with
rotation counted once by starting after that reverse arc. The missing
forward cross arc can obstruct the other boundary only when both pieces
are singletons. The stated total-order-at-least-three exception is exactly
right. The inherited three-coordinate laws also have the correct singleton
delta term and no missing empty-slot contributions.

For start-one completeness, deleting the root forces a path to begin in
the initial SCC. All start responses within that strong component are
positive. Their sum times the later-component product can be one only
when there is exactly one outneighbor port, its start response is one,
and every later SCC is a singleton. Strongness forces the last tail vertex
to point to the root, so the resulting gadget has ending response one.
The recursive component is strictly smaller. The converse grammar is
valid directly from the start and cycle product laws, including the
one- and two-vertex initial outputs.

The bridge equivalence is correct for strong C=1 tournaments. Deleting a
strong bridge b->a and taking the vertices reachable from b gives the
required two-side decomposition: otherwise a replacement b-to-a path
would restore strongness without the arc. The unique-ending-path induction
is valid on a nontrivial terminal SCC. Its chosen outdegree-one vertex
loses to the earlier SCCs and to the external root, so no hidden external
outgoing arc is introduced. The singleton terminal-component case is
handled separately. Finally, every Hamiltonian path from an outdegree-one
vertex closes into a cycle, so its start response equals C.

The original lemma wording omitted **nontrivial** before “rooted strong
tournament”; the singleton was a counterexample to that isolated statement.
The producer corrected it before signoff. Every actual application already
separated the singleton case, so the main equivalence was unaffected.

The integer gadget compiler is bijective because the final zero bit gives
the highest nonzero bit of L-1. The parity formula for the starting response
also holds at the empty-base singleton. The even and odd recurrences have
the correct ranges; in particular the excluded odd recurrence at ell=1
would create an inadmissible trailing-one word.

The bounded search is genuinely arbitrary-depth and arbitrary-gadget.
`h'>=L*R` excludes all larger L, and `R'>=h+R>R` gives a finite increasing
coordinate. A preterminal R>X cannot lead to a final h<=X. A terminal
state with R>X must, however, have its h recorded first. The implementation
does so. A concrete edge test is X=5, where h=5 is a needed terminal value
with R above the cutoff. No vector saturation or unproved stabilization
is used.

The triangle family has t disjoint induced triangle blocks, forcing at
least t feedback deletions. Removing the newly marked vertex of each
block deletes every reverse seam and leaves a transitive ordered join,
giving the matching upper bound. The order, unique-cycle claim, matrix,
and characteristic polynomial all agree. Its recurrence is valid from
t=3 with the listed initial values.

## Independent computations

`audit.py` counts each one-apex gadget by its two increasing base runs,
independently of the source's recursive compiler. All 5,001 profiles agree,
with 255 also checked by direct adjacency subset counting.

A separate heap search processes states in increasing R rather than the
source's breadth-first order. It reproduces every value, state count,
transition count, and minimum-discovery-depth statistic at bounds
1,2,3,5,47,501,5001. In particular:

- Bound 501: 2,727 states, 3,900 transitions, discovery depth 9.
- Bound 5001: 81,585 states, 117,980 transitions, discovery depth 13.

All listed missing odd values agree exactly.

Independent graph checks include 841 gluings with literal cycle counts,
235 rooted recursive decompositions, 234 instances of the corrected
nontrivial unique-ending-path lemma, and all 89 C=1 classes through order
eight for the bridge/outdegree-one/start-one equivalence. The complete
class input is from this auditor's independent primary census, not the
producer's generated tournament files. All strong classes through order
seven and all C=1 classes of order eight are covered by the decomposition
checks.

`family.py` reconstructs four triangle-family graphs through order 13,
independently counts H and C, checks the explicit feedback set, and directly
rejects every smaller candidate feedback set. It confirms H=5,35,229,1483,
C=1, and feedback vertex numbers 1,2,3,4.

Both source verifier modes pass with identical receipts. The separate
order-nine exploratory census in Section 7 was not independently rerun
in this audit; it is explicitly finite evidence and supports none of the
all-order conclusions signed off here. The global C=1 bridge question and
the bridge-free C=2 classification remain open as stated.

Results: `AUDIT.json`, `FAMILY.json`, `closure-*.json`, and associated logs.
Run `python audit.py` and `python family.py` to reproduce. The former uses
the sibling independent H67 census's `classes*.txt` files. No author source
was changed by the auditor.
