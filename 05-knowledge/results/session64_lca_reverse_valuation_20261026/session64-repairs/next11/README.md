# Complete order-eleven closed endpoint-four replacement

2026-10-26. **FINITE-EXACT:** every directed, closed, blue marked
Hamiltonian path at order11 with endpoint outdegrees `(6,4)` has
three directed black K moves to its exact labelled K target.

Combined with the separately sealed n12/n14 packet and the inherited
n13 and n>=15 theorems, this gives the uniform closed endpoint-four
replacement theorem for **n>=11**. This is a local replacement result,
not full black-core connectivity, and does not cover open marked paths.
No claim of optimality is made.

## 1. Full affine universe

Use `P=(0,...,10)`, with actual incoming ports
`N-(0)={10} union I`, where I is any three-element subset of `{2,...,9}`.
Impose exactly P's directed path arcs, closing arc `10->0`, complete
positional anti-reversal, and both branch bits `T(1,9),T(2,9)`.
There is no clone, module, strongness, disjoint-port or isomorphism filter.

All `4*binomial(8,3)=224` affine source cells are nonempty.

| Partition | Complete cells | Unchanged passing cells/leaves | Failed cells/leaves |
|---|---:|---:|---:|
|Inherited base bank|212 cells|208 cells|4 cells|
|Closed-middle endpair I={h,8,9}, h=2..7|12 cells, split into25 leaves|23 leaves|2 leaves|

The endpair leaves use the literal h2/h3/h4/h>=5 formulas inherited by
the n>=12 theorem, **rechecked at n11 rather than assumed valid there**.
Their guards are genuine source arcs. Each leaf is nonempty; every pair
of leaves in one source cell has inconsistent intersection; and the sum
of their exact powers-of-two cardinalities equals the unsplit cell's
cardinality. Thus the25 leaves form a complete disjoint partition.
Every old ledger in all these cells is universally directed and has
path-edge XOR E(P) and the exact labelled target. Only blackness can fail.

## 2. Exactly seven exceptional tournaments

Each feasible blue locus listed below has affine dimension0: the source
is uniquely forced. There are seven distinct sources, not seven samples.

| I | Middle | x,c chord | Inherited branch | Blue step(s) |
|---|---|---|---|---|
|{3,4,7}|closed|x->c|closed_one|3|
|{3,5,7}|closed|x->c|closed_one|3|
|{4,6,7}|closed|c->x|closed_zero|3|
|{5,6,7}|closed|c->x|closed_zero|3|
|{6,8,9}|closed|x->c|h6_backward,3->1|2|
|{7,8,9}|closed|c->x|h7_forward,1->3|2 and3, on two **different** sources|

Outside these loci the displayed inherited ledger is already all black.
Inside them, `found.json` supplies seven replacement triples. The checker
reconstructs every raw affine locus and verifies that the union of its
full solution sets equals exactly those seven literal source rows.
Each replacement is legal and black at every actual successive source,
and ends at the complete literal `K_P(T)`, not merely an isomorphic copy.
Six selected triples use the cyclic-shift bridge from the parent packet;
the seventh uses a short prefix reorder. This distinction is incidental
to the exact finite theorem: every triple is replayed independently.

## 3. A genuine intermediate clone merger is retained

The source

```
I={7,8,9},
T=(126,156,440,1008,1248,1474,1926,1793,1555,1079,15)
```

has the repair

```
Q1=(1,2,3,4,5,6,7,8,9,10,0),
Q2=(1,5,6,7,2,3,4,8,9,10,0),
Q3=(0,1,5,6,7,2,3,4,8,9,10).
```

T and the first intermediate state are clone-free. After Q2, the exact
ordered maximal transitive clone block is **(5,4)**. Q3 splits it again.
All other vertices remain singleton blocks. Of the28 stored states in
the seven repairs, this is the only nonsingleton clone block.

The certificate retains literal vertices4 and5 throughout, their internal
arc, every external incidence, every flipped pair, and the full block
intersection maps at merger and split. No contraction is applied and no
claim that weighted-core K moves commute with canonicalization is needed.
Every move also stores complete source rows, endpoint neighborhoods and a
specific failed anti-reversal equality proving blackness.

## 4. Verification and limits

```
python check.py > normal.json
python -O check.py > optimized.json
cmp normal.json optimized.json
```

Both runs PASS and are byte-identical. This checker imports the sealed
parent `../check.py` by a checked SHA256; that parent pins the inherited
Gaussian solver and old-bank word source. The new endpair source splits,
complete conditional-union proof, and seven-source replay are local.
The independent auditor uses its separate raw forward-pair system.

The full theorem is at one finite order. The n>=11 consequence combines
it with the separately proved n>=12 range; it is not extrapolation.
Order10 and below, open marked paths, other endpoint layers and the full
black core remain outside this result. The parent order-seven pure-blue
target continues to obstruct an unfiltered all-small-order claim.
