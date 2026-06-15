# The permanent gaps are alpha-lattice baby-Hodge holes

**Author:** codex-2026-06-15
**Builds on:** `THM-201`, `THM-202`, `THM-079`, `THM-499`, `THM-505`, and
kind-pasteur's "board = mean-field, Ω = interacting gas" synthesis.
**Status:** exact `n<=6` computation plus theorem-backed continuation.

## The reframe

The clean object is not the scalar `H` by itself. It is the low-level OCF
occupation vector.

For `n<=8`, there are no triple disjoint odd-cycle packings, so

```text
H = 1 + 2*alpha_1 + 4*alpha_2.
```

That turns the first permanent gaps into a lattice problem:

- `H=7` means `alpha_1 + 2*alpha_2 = 3`;
- `H=21` means `alpha_1 + 2*alpha_2 = 10`.

The free board / mean-field picture only sees that these are legitimate affine
targets. The interacting tournament asks a sharper question:

> is the required `alpha`-vector realized by any genuine conflict graph `Ω(T)`?

The answer is no.

## Exact n=6 picture

`alpha_lattice_mean_field_gap_codex.py` computes the full `n=6` lattice.

The mean-field envelope at fixed `alpha_1=m` is enormous:

```text
0 <= alpha_2 <= C(m,2).
```

But the realized fibers are tiny. For example:

```text
alpha_1=10 -> alpha_2={2}
alpha_1=8  -> alpha_2={0}
alpha_1=6  -> alpha_2={0,1}
alpha_1=4  -> alpha_2={0}
alpha_1=2  -> alpha_2={0,1}
```

So on the `H=21` line, every candidate point is missing:

```text
(10,0), (8,1), (6,2), (4,3), (2,4), (0,5).
```

Likewise on the `H=7` line both candidates are missing:

```text
(3,0), (1,1).
```

This is the right baby-Hodge statement:

> the first permanent gaps are not just missing scalar values; they are missing
> low-dimensional occupation vectors.

## Why this matches the known proofs

This is not replacing the existing structure. It is compressing it.

- `THM-201` says the `H=7` shape `K3` cannot occur as a component of `Ω(T)`.
  In alpha-language, that is exactly the missing point `(alpha_1,alpha_2)=(3,0)`.
- `THM-202` blocks the `P4` route toward `H=21`.
- `THM-079` gives the full jump table for the `H=21` line through `n=8`:
  the realized `alpha_2` values on each `alpha_1` fiber jump over the needed
  targets.

So the alpha-lattice picture is not a metaphor after the fact. It is a common
compression of the proofs we already have.

## The hidden tournament structure

The spectral/mean-field side sees:

```text
alpha_1 = total odd-cycle census,
c3, c5, tr(A^k), ...
```

The interacting side sees:

```text
which odd cycles can coexist,
which pairs are forced to intersect,
which local overlap patterns create new cycles.
```

That hidden structure first appears as `alpha_2`. This is why `THM-505` is the
right large-scale theorem in the background: the spectrum gives the mean-field
skeleton, and the OCF gap lives in the non-spectral correction.

So the board/Ω dictionary can be sharpened:

- board = smooth occupation envelope,
- tournament = sparse realized occupation lattice,
- forbidden values = lattice holes created by overlap closure.

## The creative route forward

The next proof step should not start from `H` alone. It should start from the
realized region of `(alpha_1, alpha_2, alpha_3, ...)`.

The proof dream is:

1. moments `tr(A^k)` define the mean-field region;
2. overlap closure laws on `Ω(T)` carve out a sparse integer sublattice;
3. permanent gaps are the first missing lattice points;
4. higher baby-Hodge holes in `(c3,c5,...,H)` are the same mechanism at the
   spectral-plus-interaction level.

That is exactly your phrasing:

> the board is the mean-field reference; the tournament is the interacting
> version; the forbidden values measure Ω's deviation from the free gas.

This session makes that sentence exact at the first nontrivial alpha-level.
