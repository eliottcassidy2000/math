# Double Round Robin as Vertex Doubling

**Session:** codex-2026-06-01-S501c
**Computation:** `04-computation/double_round_robin_blowup_s501.py`
**Canon:** THM-378

## The reframing

A tournament is a one-game round robin: every unordered pair has exactly one
directed result.  Doubling the vertices does not just make a larger tournament.
It turns every old pair into a `2 x 2` fiber block.

The important doubled blocks are not lex blocks.  Lex blowup makes the old
winner sweep the block `4-0`, so it magnifies the old hierarchy.  The SC-style
doubled block is a double round robin:

```text
old winner wins one perfect matching,
old loser wins the complementary perfect matching.
```

So the quotient pair record is always `2-2`.  The old result has not
disappeared; it has moved from score to matching parity.

## Voltage bits

For each old edge `{u,v}`, choose a voltage bit `sigma_uv`.  If `u -> v` in
the base tournament, orient

```text
u_r -> v_{r+sigma_uv}
v_r -> u_{r+sigma_uv+1}
```

for `r=0,1`.  The all-zero voltage is exactly the SC blowup from the old notes:
same sheet follows `T`, opposite sheet follows `T^op`.

Changing the names of the two sheets over a vertex is a gauge transform:

```text
sigma_uv -> sigma_uv + tau_u + tau_v.
```

The real data are triangle parities

```text
p_abc = sigma_ab + sigma_bc + sigma_ac.
```

After choosing a root, the independent parities are `p_0ij`; there are
`binom(n-1,2)` of them.

## Why this matters

This puts the old SC blowup in the same language as Tournament Analysis.
There is pairwise data on each old edge, a binary switch choosing which fiber
matching carries the old winner, and a gauge action that removes fake sheet
labels.  The leftover object is a cube of dimension `binom(n-1,2)`, exactly the
fixed-base tiling dimension.

This is the precise sense in which vertex doubling is a double round robin:
the quotient scores are all tied, but the hidden voltage parities still affect
path structure.  The S501c script shows this explicitly.  For cyclic `n=5`, all
64 gauge classes have the same near-regular doubled score sequence and the
same quotient team record `(8,8,8,8,8)`, yet their Hamiltonian path counts
range from `14937` to `15565`.

## Consequences

1. The SC blowup is the trivial-voltage point in a larger family of signed
   double-round-robin lifts.
2. Score collapse is not information loss; it is relocation from quotient
   score to voltage parity.
3. Lex blowup and SC/double-round-robin blowup are opposite extremes:
   lex stores the old result as dominance, while SC stores it as parity.
4. The equality of dimensions with the tiling cube suggests that fixed-base
   tournament coordinates can be read as hidden schedule-twist parities in a
   doubled world.

## Next checks

- Compare the voltage parity vector with the fixed base-path tiling bits under
  explicit labelled examples.
- Compute how `H(D_sigma(T))` depends on triangle parities for all labelled
  `n=5` bases, not just transitive and cyclic bases.
- Test whether pressure-DAG peel layers survive or simplify under
  double-round-robin voltage lifts.
