---
id: HYP-1786
status: OPEN
source: endpoint_transfer_witnesses_s95.py
session: codex-2026-05-29
---

# HYP-1786: Endpoint Private Child Witnesses

## Statement

For the unmerged tournament-isomorphism quotient of the fixed-path tiling cube,
every parent class at level `n` has a private odd endpoint-extension child at
level `n+1`.

Equivalently, for every tournament class `C` at level `n`, there exists a child
class `D` at level `n+1` such that:

```text
A_n(C,D) is odd
A_n(C',D) is even for all C' != C.
```

If true, the full GF(2) row rank conjecture for unmerged tournament endpoint
transfer follows immediately by THM-356.

## Evidence

Exact endpoint-transfer witness computation for transitions `2->3` through
`6->7`:

```text
tournament ranks:        [1, 2, 4, 12, 56]
tournament matchings:    [1, 2, 4, 12, 56]
tournament private rows: [1, 2, 4, 12, 56]
```

The number of parent rows with at least one private odd child equals the total
number of parent classes at every tested level.

## Contrast

The complement-merged quotient still has full row rank, but private rows stop
covering all parents:

```text
merged ranks:        [1, 2, 3, 10, 34]
merged private rows: [1, 2, 3, 9, 28]
```

The even-graph quotient loses rank:

```text
even ranks:          [1, 1, 2, 6, 8]
even matchings:      [1, 1, 2, 7, 8]
even private rows:   [1, 1, 2, 3, 7]
```

At `5->6`, even graphs have a full support matching but deficient rank. The
obstruction is algebraic cancellation, not Hall support.

## Possible Proof Route

Find a canonical extension signature for each parent tournament class that
creates a child class whose endpoint-deletion profile has exactly one odd
parent.

Potential invariants to use:

1. child `H`, `F`, and automorphism tax;
2. endpoint-deletion multiset of parent classes;
3. score sequence plus endpoint exposure;
4. a lexicographically extremal signature, such as all endpoint arcs down/up
   after putting the parent in canonical base-path form.

The computational evidence suggests the unmerged quotient may have a hidden
triangular ordering, while complement merging obscures the private witnesses.
