# Tiling Bucket Constraints S95

The useful simplification from the tournament-tiling explorer is that the
fixed-path cube is a Szele slice:

```text
2^m tilings, m = C(n-1,2).
```

Every isomorphism class bucket has size

```text
F(C) = H(C) / |Aut(C)|.
```

Since `H(C)` and `|Aut(C)|` are both odd, every unmerged bucket is odd. This
turns the explorer into a constrained partition of a power of two into odd
Hamiltonian-path orbit counts.

## Merged Bucket Shape

Complement merging has an exact 2-adic signature:

```text
self-complementary node:      odd
non-self-complementary pair:  2 * odd
```

So the merged UI is a partition

```text
2^m = odd + odd + ... + 2*odd + 2*odd + ...
```

The odd buckets are not just a statistic; they are exactly the
self-complementary classes. This gives a cheap detector: if the merged node
mass is odd, it is SC; if it is even, it is NSC and has no higher 2-adic
divisibility.

Measured sequences for `n=3..7`:

```text
V classes:                         [2, 4, 12, 56, 456]
SC classes / odd merged buckets:   [2, 2, 8, 12, 88]
merged nodes:                      [2, 3, 10, 34, 272]
SC odd mass:                       [2, 6, 52, 240, 7308]
NSC even mass:                     [0, 2, 12, 784, 25460]
```

One tiny theorem that falls out: `SC_n` must be even. The total is even, every
NSC bucket is even, and the SC mass is a sum of `SC_n` odd numbers.

## Edge Buckets

The one-tile flip cube adds another set of constraints. For merged node `u`:

```text
2 lambda_u + Σ tau_uv = m M_u.
```

Loops count twice because both endpoints stay in the bucket. Cross edges count
once at each endpoint. This gives the mod-2 boundary condition:

```text
Σ tau_uv == m M_u (mod 2).
```

When `m` is odd, odd weighted cross-incidence marks exactly the SC nodes. When
`m` is even, no node has odd weighted cross-incidence.

Measured:

```text
merged silent loop edges:          [0, 5, 20, 494, 5504]
merged cross edges:                [1, 7, 172, 4626, 240256]
odd cross-edge buckets:            [1, 1, 8, 8, 124]
odd weighted cross-incidence nodes:[2, 2, 0, 0, 88]
```

The n=7 line is especially telling: `m=15` is odd, and the 88 odd-incidence
nodes are exactly the 88 SC buckets.

## Second Moment

Merging hides complement chirality and inflates pair counts. If an NSC pair has
half-fiber `f`, then:

```text
before: f^2 + f^2
after:  (2f)^2
bonus:  2f^2
```

So:

```text
Σ_merged M^2 = Σ_unmerged F^2 + 2 Σ_NSC-pairs F^2.
```

Measured:

```text
Σ merged M²:   [2, 30, 592, 48724, 5547240]
merge bonus:  [0, 2, 52, 20924, 2360944]
```

This is a clean way to measure how much complement merging erases.

## Creative Constraints To Try Next

1. Treat the merged graph as an electrical network with node mass `M_u` and
   conductance `tau_uv`. The equation `2lambda_u + Σtau_uv=mM_u` is a weighted
   degree law, so all unknown edge buckets live in a transportation polytope.

2. Try reconstructing SC nodes from edge parity alone when `m` is odd. If node
   labels are missing but edge weights are known, parity recovers the
   self-complementary subset.

3. Look for inequalities bounding loop mass:

   ```text
   0 <= lambda_u <= C(M_u,2)
   2lambda_u <= mM_u.
   ```

   Sparse loop nodes are "locally rigid" under tile flips; heavy loop nodes are
   internal neutral zones.

4. Use the second moment and edge degree law together. The bucket masses must
   satisfy both:

   ```text
   ΣM_u = 2^m
   ΣM_u^2 = known second moment
   2λ_u + Στ_uv = mM_u
   ```

   This is much stronger than knowing the number of merged nodes.

5. The odd cross-edge-bucket count is not the same as the odd-incidence node
   count. That gap may encode how SC nodes pair through NSC corridors. For
   n=7, only 124 of 2123 cross buckets have odd weight, but they induce exactly
   88 odd-incidence SC vertices.

6. The SC odd mass `S_n` is a missing sequence. It is neither the number of SC
   classes nor the total tiling count; it is the Hamiltonian-path orbit mass
   carried by self-complementary structure. It should be compared against
   grid-symmetric tilings in the explorer.

7. Since every NSC bucket is `2*odd`, any merged bucket divisible by 4 would
   immediately falsify the fixed-path fiber model or the complement-pairing
   convention. This is a useful sanity check for future UI/export code.

The higher-level pattern is that complement merging converts a pure odd
partition of `2^m` into a parity-labeled quotient whose edge incidences still
remember the hidden odd partition.
