# Independent order-eleven audit: the entire marked endpoint-four layer closes

2026-10-26. **ACCEPTED / FINITE-EXACT:** every directed closed blue marked
Hamiltonian path P with endpoint outdegrees (6,4) in an order-eleven
tournament has three successive directed black K moves to its exact labelled
K_P target. Combining this isolated finite theorem with the sealed n>=12
result gives the same local replacement theorem for every n>=11.

No n>=12 structural theorem is silently applied at eleven. The present
checker constructs and verifies every relevant source system **at n=11**.
It does not claim that eleven is optimal or that the whole black core is
connected. This separate packet leaves the first n12/n14 audit immutable.

## 1. Complete source universe and independent implementation

Use P=(0,...,10), with a=0, b=1, x=2, y=8, c=9, z=10. A raw independent
bit encodes each of the 55 forward unordered-pair arcs. The equations impose
only the P path, z->a, the complete P anti-reversal, the exact three-element
incoming interior-port set I chosen from {2,...,9}, and the two actual branch
bits T(1,9), T(2,9).

The resulting 4*binom(8,3)=224 cells are pairwise disjoint, nonempty, and
exhaustive for the marked source universe. Every cell has affine dimension
14, so the certificate covers exactly 224*2^14=3,670,016 labelled source
tournaments with the marked path fixed. There is no strongness, mixedness,
module, clone, isomorphism, shared-port, or extra-interior filter.

The independent `check.py` reads only its local `found.json`, the producer's
seven proposed literal source/word certificates. It imports no producer
algorithm. Its raw Gaussian solver, source equations, old words, and full
h-family leaf generator are copied from the earlier sealed independent
audit, then applied directly at eleven. The manifest pins that predecessor
implementation. Every source-blue augmented system is built with the actual
preceding K flips. The same 96 small raw equation-system controls compare
Gaussian consistency and the complete decoder with all 64 assignments per
system, including both positive and inconsistent cases.

## 2. Base cells: four singleton failures and no omitted source

There are 212 base cells. All their old three-word ledgers are uniformly
directed and have the correct full symbolic labelled K_P target. In 208
cells all three actual words are universally black. The exact remaining
augmented blue systems are:

| I | Middle | Chord | Old branch | Blue step | Blue dimension |
|---|---|---|---|---:|---:|
|{3,4,7}|closed|x->c|closed_one|3|0|
|{3,5,7}|closed|x->c|closed_one|3|0|
|{4,6,7}|closed|c->x|closed_zero|3|0|
|{5,6,7}|closed|c->x|closed_zero|3|0|

Every dimension-zero system has exactly one literal source. Each entire
dimension-14 original cell is therefore settled by its old ledger outside
that singleton and by the supplied replacement at that singleton. These are
not four sampled representatives of larger unresolved families.

## 3. The previously deferred endpair cells: three more singleton sources

The remaining twelve original cells have closed middle and I={h,8,9},
2<=h<=7, with either x,c direction. The existing literal h-family words are
refined by their actual complementary chord guards into 25 nonempty affine
leaves. Their two-way intersections are checked inconsistent, and their
2^dimension cardinalities sum to the full original source cell. Thus the
refinement is both disjoint and exhaustive.

Every old leaf ledger remains uniformly directed and has the exact symbolic
target at eleven. Twenty-three leaves are uniformly three-black. The only
exceptions are:

| I | Chord | Further guard | Old ledger | Blue loci |
|---|---|---|---|---|
|{6,8,9}|x->c|3->1|h6_backward|Q2: one dimension-zero source|
|{7,8,9}|c->x|1->3|h7_forward|Q2 and Q3: one dimension-zero source each|

The last two singleton loci are **different sources**. The checker forms
their exact union and replays each old word literally; it does not conflate
“Q2 can be blue” with “Q3 is blue on the same source.” Their union has size
two. The endpair exception universe therefore has exactly three distinct
sources, not two and not an unenumerated affine family.

The source partition now closes:

    208 uniformly repaired base cells + 4 patched base cells
    + 12 completely covered endpair cells = all 224 cells.

At the refined level, there are 23 uniform endpair leaves and two patched
leaves, with the latter needing exactly three distinct literal repairs.

## 4. All seven replacements and the actual clone merger/split

The independently decoded blue-locus union equals the supplied seven source
matrices exactly: four base sources and three endpair sources, no duplicates,
omissions, or extras. On every one, the checker verifies the original source
path and anti-reversal, exact incoming and reflected outgoing port labels,
endpoint scores, branch bits, and the relevant leaf's chord guards.

Each of the 21 new words is a complete Hamiltonian permutation, directed at
its actual source, and black by an explicit anti-reversal mismatch. The
receipt records the ordered pair, its actually reflected ordered pair, their
unequal bits, the complete source row matrix at each stage, and the actual
word-endpoint incoming/outgoing ports. Three full K operations, all path-edge
XORs, and complete literal target equality are checked separately.

The clone sidecar is nontrivial in this packet. For source

    (126,156,440,1008,1248,1474,1926,1793,1555,1079,15),

the first two moves produce

    (126,184,314,1008,1472,1488,1926,1541,1667,1079,15).

Here vertices 5 and 4 have identical arcs to every other vertex and 5->4,
so the ordered pair (5,4) is a genuine two-vertex transitive clone block.
There is no twin pair at the original source or after the first move. The
third move splits this block: the exact K_P target has no twin pair. These
facts are asserted directly against the literal matrices, not inferred from
core weights. No quotient is substituted anywhere. Arbitrary clone
configurations in the much larger uniform source cells are likewise retained
because all their original pair variables remain present.

## 5. Hostile boundary, mechanism, and disposition

The known endpoint-four order-seven source
(6,4,8,19,103,79,15) still has pure-blue target
(122,124,121,112,32,64,16). All 5040 target permutations are checked; exactly
nine are directed Hamiltonian paths and none is black. Own-K invariance of
blue/black status precludes a black move entering that target. Thus a genuine
small-order obstruction remains; local ledger closure at eleven is not an
unfiltered all-order connectivity proof.

The mechanism is an exact conditional split: uniform blackness handles the
complement of each old blue locus, and raw Gaussian rank collapses the actual
blue loci to seven singletons. Their verified replacements remove the final
exceptions. The source partition, not the search menu that discovered the
words, supplies the universal quantifier.

Reproduction in this directory:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json

Both executions pass with byte-identical receipts. Gates use explicit
exceptions, not removable assertions. `MANIFEST.json` pins this report, the
checker, literal input, and both outputs. The audit makes no repository edits
or pushes.

**Final disposition:** accept complete order-eleven marked closed endpoint-four
closure and the resulting n>=11 local replacement threshold. Orders below
eleven, threshold optimality, and whole black-core connectivity are outside
this certificate.
