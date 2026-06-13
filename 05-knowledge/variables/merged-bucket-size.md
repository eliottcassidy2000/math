# Variable: Merged Bucket Size

**Symbol:** `B_M`, `W_d(M,N)`
**Type:** integer bucket size; integer transport matrix
**Defined in:** THM-345

## Definition

Fix the base Hamiltonian path `n-1 -> n-2 -> ... -> 0`, and let
`Q_m = {0,1}^m` be the staircase tiling cube with `m = C(n-1,2)`.

The merged-bucket map is

```text
pi : Q_m -> G_n/Z_2,
```

where tournament isomorphism classes are identified with their
tournament-complement classes.

For a merged node `M`,

```text
B_M = |pi^{-1}(M)|.
```

For Hamming layer `d`,

```text
W_d(M,N) = #{(t,S): pi(t)=M, |S|=d, pi(t xor S)=N}.
```

This is an ordered transition count.  The corresponding unordered line
count between distinct buckets is `W_d(M,N)` for `M<N`, since symmetry
makes the two directed counts equal.

## Values at Small n

| n | m | total tilings | merged buckets | SC buckets | NS buckets |
|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 2 | 2 | 2 | 0 |
| 4 | 3 | 8 | 3 | 2 | 1 |
| 5 | 6 | 64 | 10 | 8 | 2 |
| 6 | 10 | 1024 | 34 | 12 | 22 |

Bucket histograms through `n=6` are saved in
`05-knowledge/results/merged_bucket_constraints_s5.out`.

## Equations

- Bucket partition:

```text
sum_M B_M = 2^m.
```

- Bucket parity:

```text
B_M is odd       iff M is SC,
B_M == 2 mod 4  iff M is NS.
```

- Hamming-layer row sums:

```text
sum_N W_d(M,N) = B_M * C(m,d).
```

- Transport symmetry:

```text
W_d(M,N) = W_d(N,M).
```

- Diagonal parity:

```text
W_d(M,M) == 0 mod 2.
```

- Cross-outflow parity:

```text
sum_{N != M} W_d(M,N) == B_M * C(m,d) mod 2.
```

By Lucas' theorem, `C(m,d)` is odd iff every 1-bit of `d` is also a
1-bit of `m`.

## Relationships

- Related to `H(T)`: for an unmerged class `C`,
  `tilings(C)=H(C)/|Aut(C)|`.
- Related to `delta_proj`: bucket transport matrices are the quotient
  side of projection-defect computations.
- Related to `G_n/Z_2`: `B_M` is the stationary mass of a merged node for
  the random Hamming-layer walk on tilings.

## Tags

#merged-metagraph #tiling-buckets #hamming-layers #parity #transport
