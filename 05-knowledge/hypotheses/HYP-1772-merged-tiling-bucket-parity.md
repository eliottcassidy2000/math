---
id: HYP-1772
status: CONFIRMED
source: merged_tiling_bucket_constraints_s95.py
session: codex-2026-05-29
---

# HYP-1772: Merged Tiling Bucket Parity

## Statement

In the fixed-path tournament-tiling explorer, let

```text
m = C(n-1,2)
F(C) = number of tilings in isomorphism class C
M(U) = merged tiling mass of a complement-merged node U.
```

Then:

```text
Σ_C F(C) = 2^m
F(C) is odd for every unmerged class C
```

After complement merging:

```text
M(U) is odd  iff  U is self-complementary
M(U) = 2 * odd  iff  U is a non-self-complementary complement pair
Σ_U M(U) = 2^m
```

So the power of two is split into exactly `SC_n` odd buckets and
`(V_n-SC_n)/2` buckets with 2-adic valuation exactly 1.

## Proof Skeleton

The previous fixed-path fiber theorem gives:

```text
F(C) = H(C) / |Aut(C)|.
```

Redei's theorem says `H(C)` is odd. A tournament automorphism group has odd
order: if it had even order, Cauchy's theorem would give an element of order
2, but an involutory automorphism must swap some pair of vertices and therefore
send the edge between them to its reverse. That is impossible in a tournament.

Therefore `F(C)` is odd. Complement preserves `H` and `|Aut|`, hence preserves
`F`.

If `C` is self-complementary, its merged node has mass `F(C)`, odd. If `C` is
not self-complementary, the merged node contains `C` and `C^op`, so its mass is
`F(C)+F(C^op)=2F(C)`, exactly twice an odd number.

## Edge Balance

Let `lambda_u` be the number of one-tile cube edges staying inside merged node
`u`, and `tau_uv` the number of cube edges between merged nodes. Then:

```text
2 lambda_u + Σ_{v!=u} tau_uv = m M_u.
```

Consequently:

```text
Σ_{v!=u} tau_uv == m M_u (mod 2).
```

When `m` is odd, precisely the self-complementary merged nodes have odd
weighted cross-incidence. When `m` is even, every merged node has even weighted
cross-incidence.

## Evidence

Exact fixed-path enumeration for `n=3..7`:

```text
SC classes = odd merged buckets: [2, 2, 8, 12, 88]
SC odd mass: [2, 6, 52, 240, 7308]
NSC even mass: [0, 2, 12, 784, 25460]
merged nodes: [2, 3, 10, 34, 272]
odd weighted cross-incidence nodes: [2, 2, 0, 0, 88]
```

No failures were found in:

```text
F(C)=H(C)/|Aut(C)|
F(C) odd
|Aut(C)| odd
merged parity iff SC
NSC merged bucket v2 exactly 1
weighted edge-incidence balance
edge-incidence parity
```

## Meta-Level Reading

The merged explorer graph is not merely distributing arbitrary odd numbers into
buckets. It is distributing odd Hamiltonian-path orbit counts subject to:

1. complement pairing,
2. exact 2-adic bucket type,
3. cube-edge incidence balance,
4. parity boundary conditions on the quotient graph,
5. second-moment inflation from hiding complement chirality.

These constraints should be part of any attempt to determine the merged
structure at all `n`.
