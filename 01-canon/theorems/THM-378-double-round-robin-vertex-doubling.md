---
id: THM-378
name: double-round-robin-vertex-doubling
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S501c
depends_on:
  - THM-371
---

# THM-378: Vertex doubling is a double-round-robin voltage lift

## Statement

Let `T` be a tournament on a finite vertex set `V`, and let

```text
sigma : {{u,v} : u != v} -> F_2
```

be a voltage bit on every unordered pair.  Define a tournament
`D_sigma(T)` on `V x F_2` as follows.

First put the internal fiber arc

```text
(v,0) -> (v,1)
```

for every `v in V`.  For every base arc `u -> v` in `T` and every
`r in F_2`, put

```text
(u,r) -> (v,r + sigma_uv),
(v,r) -> (u,r + sigma_uv + 1).
```

Then:

1. `D_sigma(T)` is a tournament on `2|V|` vertices.
2. For every base pair `{u,v}`, the `2 x 2` block between the doubled fibers
   is a double-round-robin block: the winner in `T` wins one perfect matching
   of the `K_{2,2}` fiber block, and the loser wins the complementary perfect
   matching.  Thus the projected pair record is always `2-2`.
3. The clone score sequence is universal.  If `|V|=n`, every `(v,0)` has
   out-degree `n`, and every `(v,1)` has out-degree `n-1`, independently of
   `T` and `sigma`.
4. The all-zero voltage `sigma=0` is the repo's SC blowup: same-sheet arcs
   follow `T`, cross-sheet arcs follow `T^op`.
5. Relabelling sheets by a function `tau: V -> F_2` sends

```text
sigma_uv |-> sigma_uv + tau_u + tau_v.
```

The gauge-invariant voltage data are triangle parities

```text
p_abc = sigma_ab + sigma_bc + sigma_ac.
```

With a chosen root `0`, the parities `p_0ij` for `1 <= i < j <= n-1` are
independent and complete.  Hence signed double-round-robin vertex doublings
over a fixed labelled base tournament have

```text
2^binom(n-1,2)
```

sheet-flip gauge classes.

## Proof

For a fixed base pair `{u,v}`, assume first that `u -> v` in `T`.  The rule
orients the two cells

```text
(u,r) -> (v,r + sigma_uv),  r=0,1,
```

in the direction of the base winner.  These two cells form one perfect
matching of the `2 x 2` fiber block.  The second rule orients the remaining
two cells

```text
(v,r) -> (u,r + sigma_uv + 1),  r=0,1,
```

in the opposite direction; these cells form the complementary perfect
matching.  Therefore every cross-fiber pair receives exactly one orientation.
The internal arcs orient the only pairs inside a fiber.  Thus every unordered
pair in `V x F_2` receives exactly one arc, so `D_sigma(T)` is a tournament.

The same block computation proves the double-round-robin statement: each
base pair splits into two matching games won according to `T` and two matching
games won according to `T^op`.  The quotient record is therefore `2-2`, while
the original base relation is retained by which matching carries the old
winner.

For scores, fix `v`.  Vertex `(v,0)` beats `(v,1)` internally.  Against each
other base vertex `w`, exactly one of the two clone cells incident to `(v,0)`
is an out-neighbor of `(v,0)`, regardless of whether `v -> w` or `w -> v` and
regardless of `sigma_vw`.  Therefore

```text
d^+((v,0)) = 1 + (n-1) = n.
```

The same argument without the internal win gives

```text
d^+((v,1)) = n-1.
```

If `sigma=0`, the winning matching for a base arc `u -> v` is
`u_0 -> v_0` and `u_1 -> v_1`, while the complementary matching is
`v_0 -> u_1` and `v_1 -> u_0`.  This is exactly the SC blowup rule.

Finally, replace every label `(v,r)` by `(v,r + tau_v)`.  For an arc between
fibers `u` and `v`, sheet difference changes by `tau_u + tau_v`; hence the
voltage changes as claimed.  Around any triangle `{a,b,c}`, the added terms
are

```text
(tau_a+tau_b) + (tau_b+tau_c) + (tau_a+tau_c) = 0 in F_2,
```

so triangle parity is invariant.

Conversely, choose a root `0` and apply the sheet flip `tau_v=sigma_0v` for
`v != 0`, with `tau_0=0`.  This gauges every root edge voltage to zero.  The
remaining voltage on edge `{i,j}` with `i,j != 0` becomes

```text
sigma_ij + sigma_0i + sigma_0j = p_0ij.
```

Thus the root triangle parities are a complete normal form.  There are
`binom(n-1,2)` such bits, giving `2^binom(n-1,2)` gauge classes.

## Significance

Doubling vertices should be read as a hidden schedule lift, not merely as a
larger tournament.  The quotient team records are all tied, so score forgets
the old hierarchy.  The old tournament is stored in the matching parity of
each doubled block.

This explains why the SC blowup erases score variation while preserving enough
structure to affect `H`: the score layer sees only the `2-2` double
round-robin aggregate, but Hamiltonian paths still see the voltage pattern.
Modulo sheet flips, the hidden voltage cube has dimension `binom(n-1,2)`, the
same dimension as the fixed-base tiling cube for tournaments on `n` vertices.

## Related

- `04-computation/double_round_robin_blowup_s501.py`
- `05-knowledge/results/double_round_robin_blowup_s501.out`
- `07-reflections/sc-blowup-and-twin-gaining.md`
- `07-reflections/lrc-tournament-first-doubling-seam-s453.md`
- THM-371
