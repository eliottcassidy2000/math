# Endpoint Transfer Witnesses S95

The endpoint-transfer rank theorem splits into mechanisms.

The matrix laws are counting:

```text
row sum    = 2^(n-1) parent fiber
column sum = child fiber
```

The rank question is different. It asks whether parent quotient classes remain
distinguishable after reducing all endpoint-extension counts modulo 2.

## Three Mechanisms

The witness script separates three notions:

```text
private odd child  =>  immediate triangular proof
full matching      =>  support-level Hall condition
full GF2 rank      =>  actual parity injectivity
```

For unmerged tournaments, the strongest thing happens through `6->7`:

```text
tournament ranks:        [1, 2, 4, 12, 56]
tournament matchings:    [1, 2, 4, 12, 56]
tournament private rows: [1, 2, 4, 12, 56]
```

Every parent class has a private odd child. If this persists, the row-rank
conjecture has a concrete proof target: construct the private child.

For complement-merged tournaments:

```text
merged ranks:        [1, 2, 3, 10, 34]
merged matchings:    [1, 2, 3, 10, 34]
merged private rows: [1, 2, 3, 9, 28]
```

Private witnesses no longer cover every row by `5->6`, but full rank survives.
Merging hides some triangular pivots without destroying parity information.

For even graphs:

```text
even ranks:          [1, 1, 2, 6, 8]
even matchings:      [1, 1, 2, 7, 8]
even private rows:   [1, 1, 2, 3, 7]
even dependencies:   [0, 1, 1, 1, 8]
```

At `5->6`, even graphs have a full matching but rank `6 < 7`. This is the key
negative example: Hall support is not enough. The rank defect is genuine
mod-2 cancellation.

## What The Even Dependencies Say

Some even-graph parent rows are literally zero modulo 2. Their endpoint
children all appear with even multiplicity. The first examples are small cycle
orbits:

```text
n=3->4: triangle row is zero
n=4->5: 3-edge path/cycle-shape row is zero
n=6->7: eight parent rows are zero or dependent
```

This suggests that high automorphism 2-adics in the graph lens pair endpoint
signatures before the quotient sees them. The tournament lens avoids this
because tournament automorphism groups have odd order.

## New Slogan

The tournament quotient remembers endpoint parity because its symmetries are
odd. The even-graph quotient forgets endpoint parity because its symmetries
contain powers of two.

That may be the cleanest bridge between:

- Redei oddness,
- automorphism tax,
- endpoint recursion,
- and the even/cycle-space rank defect.
