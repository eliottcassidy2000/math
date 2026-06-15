# HYP-2521: Permanent Gaps as Alpha-Lattice Baby-Hodge Holes

**Status:** SYNTHESIS with exact `n<=6` evidence and theorem-backed continuation at
`n=7,8` (codex-2026-06-15).

## Claim

The permanent gaps `H=7` and `H=21` are best understood as missing low-level
`alpha`-lattice points of the interacting conflict graph `Ω(T)`, not as scalar
accidents.

For `n<=8`, `alpha_3=0`, so

```text
H = 1 + 2*alpha_1 + 4*alpha_2.
```

Here:

- `alpha_1` is the total directed odd-cycle count, the first mean-field shadow;
- `alpha_2` is the first interaction correction, counting compatible
  vertex-disjoint odd-cycle pairs.

The free board / mean-field envelope at fixed `alpha_1=m` is huge:

```text
0 <= alpha_2 <= C(m,2).
```

But the realized tournament lattice is sparse. The forbidden values occur when
the line `alpha_1 + 2*alpha_2 = (H-1)/2` misses the realized lattice entirely.

## Exact evidence at n<=6

`04-computation/alpha_lattice_mean_field_gap_codex.py` exhaustively computes the
realized fibers `(alpha_1, alpha_2)` and the finer spectral fibers
`(c3,c5) -> alpha_2`.

At `n=6`, the candidate alpha-points are:

```text
H=7:  alpha_1 + 2*alpha_2 = 3
      candidates (3,0), (1,1)

H=21: alpha_1 + 2*alpha_2 = 10
      candidates (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)
```

All of them are absent from the realized `n=6` lattice.

This is sharper than "H is missing": the *entire low-level occupation pattern*
needed to build `H=7` or `H=21` is absent.

## Continuation at n=7,8 from existing theorems

- `THM-201`: when `alpha_1=3` at `n=7`, the only realized `Ω` shape is `P3`,
  not `K3`, so `alpha_2=1`, not `0`. Thus `H=7` is exactly the missing point
  `(alpha_1,alpha_2)=(3,0)`.
- `THM-202`: `Ω(T)=P4` is impossible, blocking one of the `H=21` graph shapes.
- `THM-079`: the full `H=21` line is skipped by the realized `(alpha_1,alpha_2)`
  pairs through `n=8`. The exact jump table there is the higher-`n` continuation
  of the `n=6` lattice picture.

So the "baby-Hodge hole" interpretation is not confined to `n=6`; it matches the
known structural proofs at `n=7,8`.

## Why this matters

This reframes the OCF as:

```text
spectral moments / board counts = mean-field reference,
alpha_2 and higher alpha_j      = interaction corrections,
permanent gaps                  = missing occupation vectors.
```

The free board does not know the intersection closure laws of `Ω(T)`. The real
tournament does. The gaps measure that deviation.

In particular:

- `"Why is 7 forbidden?"` becomes:
  the free low-level occupation vector `(3,0)` is not realizable by any true
  conflict graph `Ω(T)`.
- `"Why is 21 forbidden?"` becomes:
  the whole affine line `alpha_1 + 2*alpha_2 = 10` misses the realized
  low-level lattice.

## Tournament Analysis / Assumption Challenge

The natural vertices are not runners or cycle lengths. They are the occupation
vectors:

```text
(alpha_1, alpha_2), or more finely (c3,c5,alpha_2).
```

The quotient preserves the OCF scalar `H = 1 + 2*alpha_1 + 4*alpha_2` and
destroys the underlying overlap geometry. The hidden tournament structure lives
exactly in that destroyed overlap layer.

Challenged assumption: forbidden `H` values need not first appear as missing
scalars. They can appear earlier and more cleanly as missing low-dimensional
occupation vectors.

## Open route

Promote this from synthesis to proof program:

1. Describe the realized `(alpha_1, alpha_2)` region for general `n<=8`.
2. Prove the jump law on each `alpha_1` fiber directly from overlap closure in `Ω`.
3. Lift from the low-level alpha-lattice to the general baby-Hodge moment image
   `(c3,c5,...,H)`, with `alpha_2` the first interaction correction and
   `alpha_3` the next.
