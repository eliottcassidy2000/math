# The finite two-polarity failure is not an arbitrary-reader theorem

**PROVED:** the end-only exclusion below and explicit infinite counterexamples
to a corresponding start-only exclusion. The all20-target statement in
[`reader/ENTRY_DEFECTS.json`](../../reader/ENTRY_DEFECTS.json) remains
**FINITE-EXACT** and is not enlarged here.

Write Z_s(B,w) for TT_s→B plus an outer root r beaten only by w∈B, s≥1.
Assume B is nontrivial strong. The inherited
[`two_entry_funnels_20261103/READY.md`, §4](../../../two_entry_funnels_20261103/READY.md)
classifies the end defects E_v(2)−H after the two-birth reset. It does not
simultaneously exclude the converse start defects S_v(2)−H.

## Precise end-only extension

If |B|≥4, then no vertex of Z_s(B,w) has end defect1 or2, for every s≥1.
The same conclusion holds for B=C3 when s≥2. At B=C3,s=1 the inherited
order-five H9 exception has end-ready vertices, so that boundary matters.

For s=1 the graph is exactly the standard reset output of B: first birth
r beaten only by w, then birth the prefix vertex beaten only by r.
READY§4 therefore applies when |B|≥4.

For s≥2 let the ordered prefix be a1,…,as. The outer root has end defect0.
Every old B vertex other than w has indegree at least s+2≥4, and w has
indegree at least s+1≥3. Prefix ai has indegree i. The inherited bound
b_v≥d_v(d_v−1) excludes every ai with i≥3. For a2 the stronger two-inport
bound gives b_a2≥2H(B)≥6, because its residual tournament after deleting
the two inports r,a1 is TT_(s−2)→B. Finally a1 has sole inport r, so

    b_a1 = start_r(Z_s−a1) = H(B) ≥3.

The last equality holds because a Hamilton path starting r must traverse
the remaining nonempty prefix first and then B. This proves the end-only
claim. It says nothing about small start defects at old vertices.

## A strong order-four seed defeats the two-polarity generalization

Take B with outgoing rows[8,1,3,6] on names0,1,2,3, and mark w=0. Vertex1
beats only0. For every s≥1 it still has unique out-neighbor0 in Z_s(B,0).
The converse of the single-inport defect identity gives

    S_1(Z_s,2)−H(Z_s) = end_0(Z_s−1).

An actual Hamilton path ending0 in Z_s−1 must start at r: its only possible
predecessor would be0, which is already the final vertex. It must then use
the entire transitive prefix before entering B−1. The latter is the cycle
0→3→2→0, with exactly one Hamilton path ending0. Therefore

    S_1(Z_s,2)−H(Z_s) = 1        for every s≥1.

For s=1 a literal source is rows[40,1,3,6,15,30] on
[0,1,2,3,a1,r]. Independent subset DP gives H=19 and S_1(2)=20.
Thus a converse/re-mark supplies a genuine+2 reset entry even though the
end-only theorem excludes end-ready vertices in that same graph.

The empty-reader boundary is not the explanation. Starting from the same
marked B, append any all-zero reader of length l≥2, where a zero adds a
vertex beating only the retained mark0. Call the resulting strong component
C_l and its final reader vertex z_l. This vertex still beats only0 in
Z_s(C_l,0), hence

    S_z_l(2)−H = end_0(Z_s−z_l) = end_0(C_(l−1)).

The first equality is the unique-outport identity; the second is the same
forced outer-root/prefix argument. The exact reader recurrence is e'=q
under a zero and q'=q. Here q=H(B−0)=1, since B−0 is transitive. Consequently
end_0(C_(l−1))=1 for every l≥2. Arbitrarily long nonempty readers and every
positive prefix clock therefore also admit a start-ready1 vertex.

**Scoped selector corollary.** More generally, fix any nontrivial strong
seed (B,w) and q=H(B−w). Let an arbitrary binary reader word end in00, and
let z be its final vertex. In Z_s of the resulting component, for every s≥1,

    S_z(2)−H = q.

Indeed, deleting w turns every reader birth, whether0 or1, into a global
sink birth. Thus q is invariant under the whole reader. The penultimate0
sets the marked endpoint count e_w to q, and the final0 gives z the sole
out-neighbor w. The forced outer-root/prefix argument above then gives
sigma_z=e_w=q. Ending in0 makes the component strong: w reaches every old
reader vertex through a different seed vertex, and every old vertex other
than w beats the newest vertex, which beats w. In particular q=1 is an
exact sufficient selector for a start-ready1 port after any reader ending00.
More precisely, for **any nonempty reader word ending0**, let e_seed be the
number of seed Hamilton paths ending at w. The final reader vertex z has

    sigma_z = S_z(2)−H =
        e_seed,  if the whole reader word is0;
        q,       if the word ends00;
        0,       if the word ends10.

The unique-outport/forced-prefix proof gives sigma_z=e_w of the component
before the last birth. For the empty preceding word this is e_seed; a
preceding0 sets e_w=q; a preceding1 appends a global sink, setting e_w=0.
These are all possibilities. Consequently the **terminal vertex** supplies
sigma1 iff the word is0 with e_seed=1, or it ends00 with q=1. It supplies
sigma2 iff the whole word is0 and e_seed=2, since the tournament deletion
count q is odd. No assertion is made about the other ports by this exact
terminal-vertex classification.

The penultimate reader vertex specifically cannot supply sigma1 or2 either.
Its exact out-neighbors are w and the final reader vertex. Its outdegree2
already excludes sigma1. If it had sigma2, applying READY§2 to the converse
would require a transitive residual and one of its five two-port words.
Every old prefix vertex beats both outports, so it gives a00 letter in that
converse word. The five accepted words contain no00. Since the prefix is
nonempty, sigma2 is impossible. This argument does not exclude other marks.

Actual bank0 (rows[4,1,2], mark0) provides two scoped positive entries at
prefix clock2: reader11100 has H469 and last/penultimate start defects1/168;
reader01100 has H475 and corresponding defects1/182. These are alternative
reader targets, distinct from the minimum-order selected targets whose
all-port defects were recorded in `ENTRY_DEFECTS.json`.

The first missing implication is precisely the attempt to transfer an
end-defect classification to start defects in the same actual tournament.
Converse swaps those predicates but also changes the source architecture.
Additional hypotheses on the actual seed and reader, or explicit all-port
readouts as supplied for the20 targets, are required. No broader reader
classification or new search is claimed here.

## Portable literal controls

`check_entry_scope.py` independently computes subset path endpoints on ten
explicit matrices: the H19 source above, the C3/one-prefix H9 boundary, and
the same strong four-vertex seed with reader00 and prefix clock2. The last
case checks start defect1 at the newest reader vertex. Two further controls
use the mixed words10100 at mark0 (q=1) and100 at mark1 (q=3), both at prefix
clock1, checking the generalized selector's exact value. The final two
controls reproduce the explicit H469/H475 alternative entries and their
last/penultimate defects. Three terminal-law controls use the lone0 at
marks0 and1 (seed endpoint counts1 and2), and10 at mark0 (defect0).
No search is used.
Reproduce from repository root:

    python 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/check_entry_scope.py
    python -O 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/check_entry_scope.py --output 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/entry_scope_optimized.json

The normal and optimized receipts must be byte-identical.
