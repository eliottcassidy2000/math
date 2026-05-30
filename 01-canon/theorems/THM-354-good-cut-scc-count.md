---
id: THM-354
name: good-cut-scc-count
status: PROVED
date: 2026-05-30
session: codex-2026-05-30-support-residue
depends_on:
  - THM-330
scripts:
  - 04-computation/goodcut_scc_defect_s354.py
results:
  - 05-knowledge/results/goodcut_scc_defect_s354.out
---

# THM-354: Good-Cut Count Equals n Minus the SCC Count

## Statement

Let `T` be a tournament on `n` vertices and let

```text
P = (v_0, v_1, ..., v_{n-1})
```

be any directed Hamiltonian path of `T`. A cut `k in {1,...,n-1}` of this path
is good if some arc crosses it backward:

```text
v_j -> v_i for some i < k <= j.
```

Let `c(T)` be the number of strongly connected components of `T`. Then the
number of good cuts is

```text
g_P(T) = n - c(T).
```

In particular, `g_P(T)` is independent of the Hamiltonian path `P`. The bad
cuts are exactly the boundaries between consecutive strong components in the
condensation order.

## Proof

The condensation of a tournament is an acyclic tournament, hence a transitive
tournament. Therefore the strong components of `T` have a unique linear order

```text
C_1 -> C_2 -> ... -> C_c
```

and every edge between `C_a` and `C_b` points from `C_a` to `C_b` when `a < b`.

Any directed path can only move forward in this component order. Hence any
Hamiltonian path visits all vertices of `C_1`, then all vertices of `C_2`, and
so on, with no interleaving between components. The cuts between consecutive
component blocks have no backward crossing edge from the suffix to the prefix,
so these `c(T)-1` cuts are bad.

It remains to show that every cut strictly inside one strong component is good.
Let such a cut split the current component `C` into a nonempty prefix part `A`
and a nonempty suffix part `B`. If no edge points from `B` to `A`, then every
cross edge points from `A` to `B`. A directed path from a vertex of `B` to a
vertex of `A` would have to cross from `B` to `A` at some first step, impossible.
This contradicts strong connectivity of `C`. Thus some edge points from `B` to
`A`, giving a backward crossing of the path cut.

So the only bad cuts are the `c(T)-1` component boundaries. Since there are
`n-1` cuts total,

```text
g_P(T) = (n-1) - (c(T)-1) = n - c(T).
```

## Consequences

- THM-330 is the special case `c(T)=1`: a tournament is strongly connected iff
  every cut of a Hamiltonian base path is good.
- The bottom bucket `g=0` is the transitive-condensation extreme `c(T)=n`,
  which forces every strong component to be a singleton.
- HYP-1764 is explained structurally: `g` is an isomorphism invariant because
  it is just `n-c(T)`. Taking the opposite tournament reverses the condensation
  order but preserves the strong components, so `g` also descends to complement
  merged classes.
- HYP-1770 and HYP-1777 follow for any move family: if `Delta g != 0`, then the
  source and target tournaments cannot be isomorphic or complement-isomorphic.

## Verification

`04-computation/goodcut_scc_defect_s354.py` performs two independent checks.

1. For all labeled tournaments on `n <= 6`, it enumerates every Hamiltonian
   path and verifies that every path has `g = n - c(T)`.
2. For the fixed base-path tiling model, it verifies all tilings through `n=7`
   and samples `n=8`.

It also samples random labeled tournaments at `n=7,8`. No violations were
found.

## Related

- THM-330: SC cut theorem in the tiling model.
- THM-336: good-cut bucket structure.
- THM-349: interval-union recurrence for the raw good-cut bucket counts.
- `05-knowledge/variables/good-cut-count.md`.
