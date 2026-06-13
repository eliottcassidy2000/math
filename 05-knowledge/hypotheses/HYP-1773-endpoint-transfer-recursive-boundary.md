---
id: HYP-1773
status: PARTIALLY_CONFIRMED
source: endpoint_transfer_bucket_recursion_s95.py
session: codex-2026-05-29
---

# HYP-1773: Endpoint Transfer Recursive Boundary

## Exact Theorem

Let `F_n(C)` be the fixed-path tiling fiber of isomorphism class `C` at `n`.
Add a new terminal vertex to the fixed path:

```text
0 -> 1 -> ... -> n-1 -> n.
```

Each parent tiling has `2^(n-1)` endpoint-extension children, because the only
new free arcs are from the new terminal vertex to old vertices `0..n-2`.
Therefore:

```text
Q_{m(n+1)} = Q_{m(n)} x Q_{n-1}
m(n) = C(n-1,2).
```

Let `A_n(C,D)` count endpoint-extension children of parent class `C` landing in
child class `D`. Then:

```text
Σ_D A_n(C,D) = 2^(n-1) F_n(C)
Σ_C A_n(C,D) = F_{n+1}(D)
```

After complement merging, with merged bucket masses `M_n(U)`:

```text
Σ_V B_n(U,V) = 2^(n-1) M_n(U)
Σ_U B_n(U,V) = M_{n+1}(V)
```

These are exact quotient forms of the cube product.

## Parity Boundary

Modulo 2, every parent row has even sum because `2^(n-1)` is even. In the
unmerged transfer, every child column has odd sum because every `F_{n+1}(D)` is
odd. Hence:

```text
XOR(parent parity rows) = all child classes.
```

In the merged transfer:

```text
XOR(parent parity rows) = SC child-node indicator.
```

So self-complementary child nodes can be recovered from endpoint-extension
parity alone.

## New Open Conjecture

The parity transfer matrix has full row rank over `GF(2)` for all `n`.

Measured ranks for transitions `n -> n+1`, `n=2..6`:

```text
unmerged ranks: [1, 2, 4, 12, 56]
parent classes: [1, 2, 4, 12, 56]

merged ranks:   [1, 2, 3, 10, 34]
parent merged:  [1, 2, 3, 10, 34]
```

If true, endpoint insertion loses no parent-level parity information. The
quotient tower is recursively injective over `GF(2)` even though it is massively
many-to-one over integers.

## Evidence

Exact enumeration through child level `n=7`:

```text
class support:       [2, 5, 22, 181, 2338]
merged support:      [2, 5, 19, 112, 1312]
class odd entries:   [2, 4, 16, 124, 1390]
merged odd entries:  [2, 4, 12, 62, 668]
SC boundary weights: [2, 2, 8, 12, 88]
```

All row/column sum checks and parity-boundary checks passed.

## Relation To Older Threads

This clarifies the "fast time scale" from the two-reductions thread: Mode A
vertex insertion is not only a local Hamiltonian-path update, it is a global
transfer operator between quotient levels.

It also links to the bridge-exposure thread. Local insertion of an arbitrary
vertex needs the richer boundary object `Q_T`; endpoint insertion is the
fixed-path slice where those local signatures assemble into a clean quotient
transfer matrix.
