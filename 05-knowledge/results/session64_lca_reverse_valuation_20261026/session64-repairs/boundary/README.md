# The sharp closed endpoint-four threshold is eight

2026-10-26. **FINITE-EXACT:** every directed closed blue marked path with
endpoint outdegrees `(n-5,4)` at n=8 or n=9 has a three-black-move
replacement to the exact labelled K target.

Combining this finite head with the independently checked n10/n11
packets, the sealed n12/n14 repairs, and the inherited n13 and n>=15
results gives the uniform theorem for **every n>=8**. The threshold is
sharp for this unfiltered closed endpoint-four statement: an explicit
order-seven source has a pure-blue target with no black move entering it.
This is not full black-core connectivity and says nothing about replacing
all open blue paths or other endpoint-degree layers.

## 1. Complete finite source universe

Use `P=(0,...,n-1)` as literal labels. Impose its path and closing arcs,
complete actual anti-reversal, exactly three middle incoming ports I at0,
and the two genuine source bits `T(1,n-2),T(2,n-2)`. Every I from
`{2,...,n-2}` is used. No clone, module, strongness, shared-port or
isomorphism filter is imposed.

| n | All source cells | Base cells | Endpair cells | Endpair guard leaves |
|---|---:|---:|---:|---:|
|8|40|34|6|13|
|9|80|72|8|17|

Every source and every guard leaf is nonempty. The guard leaves are
pairwise disjoint within each endpair source cell, and their exact
powers-of-two sizes sum to the whole original cell. Thus no source
completion is omitted.

## 2. Order nine: five singleton blue loci

At n9 all 17 endpair leaves pass unchanged, as do 67 of 72 base cells.
The remaining five base cells have one possible blue first word:

| I | Middle | x,c chord | Old branch |
|---|---|---|---|
|{2,3,4}|open|x->c|open_release|
|{2,3,5}|open|c->x|open_release|
|{2,3,5}|closed|c->x|closed_zero|
|{3,4,7}|open|x->c|open_release|
|{3,5,7}|open|c->x|open_release|

Every corresponding blue locus is an affine singleton. Off those five
sources the old ledger is already three-black; each singleton has a new
literal three-black certificate in `found.json`.

## 3. Order eight: repair entire illegal leaves, not just blue loci

At n8 all 34 base cells pass, and 11 of 13 endpair guard leaves pass.
The two exceptions have

```
I={4,5,6},       6->1,       3->1,
```

with either value of the 2,6 chord. The inherited h4-backward word tries
`0->5` in both its first and second moves, but5 is an actual incoming
port. These are **legality failures**, not symmetry failures. Each full
failed leaf has dimension 3, hence eight completions.

The checker enumerates the **entire original guarded affine leaf** in
these two cases, not a blue sublocus. All 16 literal sources get three-black
repairs. This distinction is essential to the full order-eight result.
Together with the five n9 singletons, exactly 21 exceptional sources are
replayed. The source-union equality is checked in both directions.

## 4. Actual large clone mergers are kept literal

For example, the order-eight source

```
T=(14,4,8,18,167,207,159,15),       I={4,5,6}
```

starts with the ordered clone block `(7,0)`. Its chosen repair has,
successively, the clone blocks `(4,5)`, then the size-four ordered block
`(7,3,0,1)`, and finally two blocks `(2,3)` and `(4,5)` in the target.
Every other block at those states is singleton. These are actual changes
in the literal tournament, not formal weighted-core substitutions.

Each certificate stores all tournament rows and every word, flipped
pair, endpoint neighborhood, blackness witness, ordered clone block,
and block-intersection transport map. Endpoints and original port labels
are never collapsed. A size-four merger may cross the marked path's old
boundaries; quotient path lifting is neither used nor presumed valid.

## 5. Exact order-seven obstruction and sharpness

Retain

```
P=(0,1,2,3,4,5,6),
T=(6,4,8,19,103,79,15),
B=K_P(T)=(122,124,121,112,32,64,16).
```

P is directed, closed and blue at T, with endpoint scores `(2,4)`.
The checker tests all 5040 permutations of B. Exactly nine are directed
Hamiltonian words, and all nine are blue. Blue/black status of a word
is unchanged across its own K edge, so no black K move can enter B.
Consequently there is no positive-length black detour from T to B.
This is a genuine connectivity obstruction for this local edge, unlike
any of the failed menus repaired above. It proves that 8 cannot be lowered
in the unfiltered all-n statement.

## 6. Reproduction and exact status

```
python full_check.py > full-normal.json
python -O full_check.py > full-optimized.json
cmp full-normal.json full-optimized.json
```

Both receipts are PASS and byte-identical. Gates remain active under
optimized Python. The new checker imports the sealed parent Gaussian
engine and old-bank source through a checked SHA256 of `../check.py`;
new guard partitions, whole-illegal-leaf enumeration, literal replay and
the order-seven obstruction check are local. The independent review is
a separate reconstruction, not another invocation of this implementation.

The conclusion is an exact finite head plus already proved higher-order
tail. It is not a sampled tournament census, a conjectured threshold,
or a proof of the unrestricted black core's connectivity.
