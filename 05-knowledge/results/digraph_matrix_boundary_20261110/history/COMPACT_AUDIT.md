# Independent review of the oriented compilers

2026-11-10. **PROVED constructions; FINITE-EXACT independent controls.**
This note audits the root session's constructions. Vertex splitting is a
standard gadget pattern; no external priority claim is made. The universal
oriented spectrum is already proved in Session235 PB5. The improvement here
is the explicit logarithmic-order strong witness for every target at least3,
obtained from the already proved binary semicomplete construction.

## 1. Retained inheritance and corrections

The owner’s proposed condition was already refuted and sharply repaired in
[Session234, DG1–DG5](../../session234_digraphs_20261109/PROOFS.md).
The complete ordinary-seven and ordinary-twenty-one semicomplete families
are proved there at all orders. Loops do not enter the count.

[Session235, PB1–PB7](../../session235_arc_pair_models_20261109/theory/PROOFS.md)
already distinguishes all three interpretations of “11 means empty,” proves
complement–reverse parity, gives all nonnegative oriented counts, and gives
unbounded strong oriented examples at each fixed total at least2. Its
[core transport](../../session235_arc_pair_models_20261109/census/CORE_TRANSPORT.md)
already proves ordinary H7 maps to decoded H1 and ordinary H21 maps exactly
to {1,3,5,7,15}. Our initial small reread independently reproduced the core
histograms, but contributes no new classification.

The binary semicomplete source is
[THM-5071, statement(4a)](../../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md).
MISTAKE-748 warns that its one-reciprocal exclusions of1 and7 are inherited
from Session234; only the exclusion of4 was new there. Neither that theorem
nor the present compiler proves the tournament spectrum conjecture.

## 2. Generic exact compiler on 2n+1 vertices

Let D be any loopless digraph on n>=1 vertices. Form S(D) on vertices
`s` and `v_in,v_out` for every old vertex v. Its arcs are exactly

    s -> v_in,
    v_in -> v_out,
    u_out -> v_in whenever u -> v is an arc of D.

There are no loops or reciprocal pairs. Input loops must be discarded:
otherwise a loop would create a reciprocal pair on its two split copies.

Every Hamiltonian path starts at s, which has indegree zero. Each v_in has
the unique outgoing neighbor v_out. It cannot be last: v_out's only
incoming neighbor is v_in, so v_out could not already have occurred.
Therefore the path has the form

    s, p1_in, p1_out, ..., pn_in, pn_out,

where p1,...,pn is a Hamiltonian path of D. Compression and expansion are
inverse bijections, so H(S(D))=H(D). One extra source is sufficient; no
terminal vertex is needed. The empty input needs a separate convention.

Without the extra source, the same split graph satisfies the exact boundary

    H(split(D)) = H(D) + n C_Ham(D),

where C_Ham counts directed Hamiltonian cycles up to cyclic rotation.
Paths starting at an in-copy encode ordinary D paths. A path starting at
v_out must end at v_in, since all other in-copies must be followed by their
out-copies; it encodes a D Hamiltonian cycle rooted at v. This explains
precisely why an unforced split is not a count-preserving compiler.

## 3. Strong logarithmic-order compiler

Fix h>=3 and put k=bitlength(h-2). Let S be the subset of {0,...,k-1}
specified by

    h-2 = sum_(i in S) 2^(k-1-i).

The highest bit is present, so 0 belongs to S. On the k+2 vertices
`0,...,k-1,a,b`, take precisely the following arcs:

    i -> j for 0 <= i < j < k;
    i -> a for every i < k;
    a -> b;
    b -> i for i in S.

The graph is oriented. It is strong because every chain vertex reaches
a then b, and b reaches0 then every chain vertex. Compressing a,b into
one vertex w gives the transitive (k+1)-vertex tournament with sink back
arcs at S. THM-5071(4a) gives its path total h-1.

Hamiltonian paths using a->b compress bijectively to those h-1 source paths.
If a->b is not used, a must be last because it has no other outgoing arc;
b must be first because it has no other incoming arc. The chain order is
then forced, giving exactly

    b,0,1,...,k-1,a.

This path exists since0 is in S. Hence the new graph has exactly h paths.
Its order is

    bitlength(h-2)+2 = ceil(log2(h-1))+2.

This is an explicit universal upper bound, not a minimum-order theorem.
At h=7 it gives an optimal order-five witness; at h=21 it uses order seven,
whereas Session235 already has an order-six witness. H0 is an edgeless pair;
H1 is a singleton; the inherited pair of triangles sharing one vertex is a
strong oriented H2 witness. The displayed binary recipe itself has domain
h>=3, and no behavior outside that domain is inferred.

## 4. Independent exact checks

`compact_audit.py` imports no producer code. For all4165 loopless labelled
digraphs of orders1–4, it compares the input's direct permutation count
with an independently implemented subset-DP count on the generic output.
It also verifies every lifted path, the unsourced Hamilton-cycle correction,
orientation, and invariance under adding all input loops.

The strong compiler is independently reconstructed for every h=3,...,512
and selected points surrounding powers of two through4099. For every one,
the auditor checks orientation, strong connectivity, target subset-DP count,
and compressed-source count. At target orders at most8 it additionally
lists all vertex permutations and checks the unique path avoiding a->b.
The normal and optimized JSON receipts must be byte-identical.

From the repository root:

    python3 -B 05-knowledge/results/digraph_matrix_boundary_20261110/history/compact_audit.py
    python3 -B -O 05-knowledge/results/digraph_matrix_boundary_20261110/history/compact_audit.py

The all-order claim follows from the bijections above; the finite checks
audit their implementation and their boundary cases.
