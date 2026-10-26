# Independent lower-boundary audit: the endpoint-four threshold is exactly eight

2026-10-26. **ACCEPTED / FINITE-EXACT:** complete closed marked endpoint-four
replacement at n=8,9,10. Together with the separately sealed n>=11 result,
every such marked blue path at every n>=8 has three successive directed black
K moves to its exact labelled K_P target. The explicit failure at n=7 makes
8 the optimal integer threshold for this **unfiltered local replacement
statement**. Whole black-core connectivity remains OPEN.

This packet is separate from the sealed n12/n14, n11, and bridge audits. It
performs no search below eight. The known order-seven example is checked only
as the sharpness obstruction.

## 1. Full source universe, not a successful word menu

Use the positions of the marked path as labels P=(0,...,n-1), and put
a=0, b=1, x=2, y=n-3, c=n-2, z=n-1. This is a bijective coordinate choice;
the resulting words can be translated back to the original vertex labels.
It is not a quotient or an identification of different tournaments.

The independent checker creates one raw bit for every forward unordered
pair. It imposes only the directed P path, closing arc z->a, full actual
P anti-reversal, an exact three-element incoming interior-port set I from
{2,...,n-2}, and the branch bits T(b,c),T(x,c). All these conditions are
signed equalities between raw bits and/or a constant. There is no strongness,
mixedness, clone, module, isomorphism, disjoint-port, or extra-interior filter.

The signed-component solver is independent of the producer's Gaussian
elimination. Its complete source partitions are:

| n | All source cells | Dimension in every cell | Literal marked sources |
|---|---:|---:|---:|
|8|40|4|640|
|9|80|7|10,240|
|10|140|10|143,360|

Every cell is nonempty. The cells are exhaustive and disjoint because the
source has one exact I and one value for each branch bit. Their count is
4*binom(n-3,3). These are labelled tournaments with P fixed, not isomorphism
classes or Hamiltonian-path counts.

The checker imports only local independent `core.py` word/literal helpers,
copied from the first sealed audit; it imports no producer source-family,
elimination, enumeration, or search code. The sole producer data are thirty
proposed literal source/word certificates in local `found.json`.

## 2. Every old-family completion is either uniformly repaired or explicitly listed

The old base family omits only closed-middle cells I={h,y,c}. Their inherited
literal h-family words are refined by actual complementary chord guards.
Each refined leaf is checked nonempty. Every pair of leaves inside an
original cell is checked inconsistent together, and the sum of their exact
2^dimension cardinalities equals the cardinality of the unsplit source cell.
Thus these are complete, disjoint refinements, not selected examples.

| n | Base cells | Endpair cells | Total refined leaves | Uniform three-black leaves | Leaves needing replacement | Literal sources repaired |
|---|---:|---:|---:|---:|---:|---:|
|8|34|6|47|45|2|16|
|9|72|8|89|84|5|5|
|10|130|10|151|142|9|9|

In each uniform leaf, every word is a full permutation, every join is forced
directed at its actual stage, and adjoining the complete actual anti-reversal
equations is inconsistent. The full symbolic flip set equals K_P, separately
from legality and blackness.

### Order eight: a legality failure requires the whole leaf

The two failed leaves have closed middle, I={4,5,6}, guard 3->1, and either
direction on x versus c. They use the old `h4_backward` words. Vertex5 is
now an incoming port, so their proposed join0->5 is forced false. Both the
first and second old words are illegal; the third formal word is directed.
The formal edge XOR still equals E(P), but that does **not** make the old
sequence a legal detour.

Each entire source leaf has dimension3. The independent signed-component
decoder enumerates all eight assignments in each leaf, producing sixteen
distinct literal sources. No blue-word condition is added. This is crucial:
repairing only a blue locus would not cover a legality failure. Every one of
the sixteen sources has a verified new three-black detour.

### Orders nine and ten: singleton blue loci

At n=9, exactly five old base leaves can have a blue first word; every blue
augmented system has dimension0. All old joins and targets are valid. The
port/branch combinations are:

* I={2,3,4}, open middle, x->c;
* I={2,3,5}, closed middle, c->x;
* I={2,3,5}, open middle, c->x;
* I={3,4,7}, open middle, x->c;
* I={3,5,7}, open middle, c->x.

At n=10, eight base leaves and one h3-forward endpair leaf have a possible
blue word, again with dimension0 in all nine cases. Their complete exact
list and literal sources are retained in the receipt. Every old join and
target is valid. The distinct union has nine sources.

For a legal old leaf, every source outside its exact blue-locus union keeps
the old three-black ledger. Every source inside that union has a new literal
repair. The independently decoded exception union at orders8/9/10 equals
the thirty supplied source matrices exactly, without omission, duplicate,
or extraneous input. This source partition is what supplies the universal
quantifier; no assumption about completeness of a search menu is needed.

## 3. Ninety actual black moves; clone mergers retained exactly

Each literal repair is checked from its original source: full tournament
matrix, P arcs and anti-reversal, endpoint scores (n-5,4), exact incoming and
reflected outgoing ports, both branch bits, and every leaf guard. Each new
word is a directed Hamiltonian permutation in its actual successive source.
All ninety words have explicit anti-reversal mismatch witnesses. The receipt
records the ordered pair, its actual reversed positional image, the unequal
arc bits, every intermediate matrix, and each word's endpoint-port labels.
The final complete rows equal K_P(T), not merely an isomorphic tournament.

Clone blocks are recomputed from their definition at **every** literal
source, intermediate state, and target. The checker exhausts every nonempty
vertex subset, tests whether it is a transitive module, retains exactly the
maximal ones, records their actual internal vertex order, checks that they
partition the vertices, and reconstructs every original arc from those
ordered blocks and their literal quotient incidences. This is a sidecar
validation, not a quotient replacement for any K operation.

For example, the order-eight source

    (14,4,8,18,167,207,159,15)

has maximal nontrivial clone block (7,0). After the first move its nontrivial
block is (4,5). After the second move the block is (7,3,0,1). At the exact
target the nontrivial blocks are (2,3) and(4,5). The literal receipt retains
all singleton blocks as well. Such mergers and splits show concretely why
one must not assume a fixed weighted clone quotient commutes with K.

## 4. Sharpness and exact scope

The known order-seven source and target are independently replayed:

    T=(6,4,8,19,103,79,15),
    K_P(T)=(122,124,121,112,32,64,16).

P=(0,...,6) is directed, closed, and blue, with endpoint degrees (2,4).
All 5040 target permutations are enumerated. Exactly nine are directed
Hamiltonian paths, and all nine are blue. Blue/black status of a word is
unchanged across its own K edge, so no black edge can enter this target.
There is no positive-length black detour to it.

The complete n8/n9/n10 finite certificate plus the separately audited n>=11
theorem therefore gives all n>=8, while the failure at 7 prevents any smaller
integer uniform threshold. This is the sharp threshold for replacing a
**closed blue endpoint-four marked move**. It does not say every blue move
lies in this layer or is closed, and it does not prove the whole black core
is connected.

## 5. Reproduction and final disposition

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json

Both runs pass with byte-identical receipts. All gates are explicit
exceptions rather than removable assertions. The manifest pins the checker,
independent local core, input table, both receipts, and this report.

**Final disposition:** accept the thirty new low-order repairs, complete
marked closed endpoint-four closure at 8/9/10, and the sharp n>=8 consequence
after combining the earlier audited tail. No search or claim below the
known order-seven sharpness obstruction is made. No repository edit or push
was performed by this audit.
