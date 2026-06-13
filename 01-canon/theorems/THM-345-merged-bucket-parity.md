---
id: THM-345
name: merged-bucket-parity
status: PROVED
date: 2026-05-29
session: kind-pasteur-2026-05-29-S5
scripts:
  - 04-computation/merged_bucket_constraints_s5.py
results:
  - 05-knowledge/results/merged_bucket_constraints_s5.out
---

# THM-345: Merged Bucket Parity and Hamming-Layer Transport Constraints

## Statement

Fix the base Hamiltonian path `n-1 -> n-2 -> ... -> 0`, and let
`Q_m = {0,1}^m` be the staircase tiling cube, where
`m = C(n-1,2)`.

Let

```text
pi : Q_m -> G_n/Z_2
```

send a tiling to its merged tournament isomorphism class, where a class
is merged with its tournament complement.  For a merged node `M`, define
the bucket size

```text
B_M = |pi^{-1}(M)|.
```

Then:

1. If `M` is self-complementary, then `B_M` is odd.
2. If `M` is non-self-complementary, then `B_M == 2 mod 4`.
3. The buckets partition the fixed-base tiling cube:

```text
sum_M B_M = 2^m.
```

For each Hamming layer `d = 1,...,m`, define the ordered transport matrix

```text
W_d(M,N) = #{(t,S): pi(t)=M, |S|=d, pi(t xor S)=N}.
```

Then:

1. `W_d` is symmetric.
2. `row_sum_M(W_d) = B_M * C(m,d)`.
3. `W_d(M,M)` is even for every bucket `M`.
4. The cross-outflow parity is

```text
sum_{N != M} W_d(M,N) == B_M * C(m,d) mod 2.
```

Equivalently, a bucket has odd cross-outflow in layer `d` exactly when
`M` is self-complementary and `C(m,d)` is odd.  By Lucas' theorem, the
second condition is `d & ~m = 0` in binary.

## Proof

For an unmerged tournament isomorphism class `C`, the fixed-base tiling
count is

```text
tilings(C) = H(C) / |Aut(C)|.
```

This is the orbit-stabilizer count for labeled tournaments in `C` whose
fixed base path is Hamiltonian.

By Redei's theorem, `H(C)` is odd.  Also, `|Aut(C)|` is odd.  Indeed, if
`Aut(C)` had even order, Cauchy's theorem would give an automorphism of
order 2.  Such an involution has a transposed pair `{u,v}`; but it would
send the oriented edge `u -> v` to `v -> u`, contradicting preservation
of the tournament orientation.  Therefore `H(C)/|Aut(C)|` is odd.

If the merged bucket `M` is self-complementary, it consists of one
unmerged class, so `B_M` is odd.  If `M` is non-self-complementary, it
consists of a complement pair `C, C^op`.  Complement preserves both
Hamiltonian path count and automorphism group size, so

```text
B_M = tilings(C) + tilings(C^op)
    = 2 * H(C)/|Aut(C)|
    == 2 mod 4.
```

This also proves the partition identity, because every fixed-base tiling
belongs to exactly one merged bucket.

For the Hamming-layer matrix, each tiling has exactly `C(m,d)` layer-`d`
neighbors, proving the row-sum identity.  The map

```text
(t,S) |-> (t xor S, S)
```

is an involution between ordered transitions from `M` to `N` and ordered
transitions from `N` to `M`, so `W_d` is symmetric.  When `M=N`, the same
map pairs all diagonal transitions with no fixed points because `d>0`;
hence `W_d(M,M)` is even.  Subtracting the even diagonal from the row sum
gives the cross-outflow parity formula.

Lucas' theorem gives the binary criterion for `C(m,d)` odd.

## Exact Verification

`04-computation/merged_bucket_constraints_s5.py` independently enumerates
the fixed-base tiling cube for `n=3,4,5,6`, constructs the merged buckets,
and checks every Hamming layer.

Summary:

| n | m | tilings | merged buckets | SC buckets | NS buckets | violations |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 2 | 2 | 2 | 0 | 0 |
| 4 | 3 | 8 | 3 | 2 | 1 | 0 |
| 5 | 6 | 64 | 10 | 8 | 2 | 0 |
| 6 | 10 | 1024 | 34 | 12 | 22 | 0 |

The script checks:

- bucket-size parity;
- exact bucket sum `2^m`;
- symmetry of every `W_d`;
- row sums `B_M*C(m,d)`;
- even diagonal;
- cross-outflow parity.
- aggregate cross-lines by SC-SC / SC-NS / NS-NS type.

## Consequences

1. **The old S202 "merged tiling excess" narration is corrected.**
   Merged buckets still partition the fixed-base tiling cube exactly.
   There is no excess term.

2. **SC status is visible from bucket size parity alone.**
   A merged node is SC iff its bucket size is odd.

3. **Every Hamming layer has a parity lower bound.**
   If `C(m,d)` is odd, each SC bucket has odd cross-outflow, so layer `d`
   contains at least `SC_merged(n)/2` cross-bucket tiling lines.

4. **Complement-tiling lines are constrained.**
   Since `C(m,m)=1`, every SC bucket has odd cross-outflow under the
   all-tiles flip.  Thus complement-tiling cannot be purely internal on
   the SC spine.

5. **The bucket matrices form a reversible transport system.**
   Normalizing rows by `B_M*C(m,d)` gives the random Hamming-layer walk on
   merged buckets, with stationary measure proportional to `B_M`.

6. **The sea appears as transport excess.**
   At `n=6`, NS-NS cross-lines dominate every layer in the S5 table,
   even in Lucas-active layers where the parity obligation comes only
   from SC buckets.  This separates the unavoidable parity constraint
   from the geometric bulk of the merged metagraph.

## Related Objects

- `B_M` and `W_d`: see `05-knowledge/variables/merged-bucket-size.md`.
- General quotient bucket balance: THM-346 and
  `05-knowledge/variables/tiling-bucket-balance.md`.
- Projection-defect work: `04-computation/projection_defect_waggly_layers_s1.py`.
- Merged tiling count precursor: `04-computation/merged_tiling_counts_s202.py`.
- Tiling count identity: `05-knowledge/results/tiling_count_theorem.md`.
