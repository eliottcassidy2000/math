---
id: HYP-2185
status: OPEN synthesis from S616; exact boundary-count identity verified, transfer theorem open
source: user-2026-06-03; codex-2026-06-03-S616; MathOverflow 299504
tags: [infinite-go, open-games, ordinal-values, n-plus-2, partition-functions, tournaments, LRC, coimage]
---

# HYP-2185: infinite Go is the open-game face of the n->n+2 boundary partition recursion

Namespace note: this S616 packet was restored after `origin/main` had already
claimed HYP-2182 for the broader partition-function/free-energy synthesis and
HYP-2183 for the LRC `n+2` automaton. HYP-2185 is the narrower boundary-state
version: it keeps the finite-cutoff Go partition function and the exact
`2n-1` boundary tile identity as its reserved claim.

Hamkins's infinite-Go construction gives a position where a black group is
doomed, but can delay death by choosing an arbitrarily long ladder; the
life/death value is `omega`. His update says that if the main group has `n+2`
liberties, the value can be made `omega*n`: two liberties are the terminal
source/sink boundary, and the remaining `n` liberties are serial announcement
fuel.

S616 reads this as the open-game analogue of the repo's `n->n+2` recursion:

```text
two terminal boundary resources + n active fuel cells.
```

With a finite ladder cutoff `K`, one fuel cell has partition function

```text
G_K(q) = 1 + q + ... + q^K.
```

For `r = liberties - 2` serial fuel cells, the truncated partition function is

```text
Z_{r,K}(q) = G_K(q)^r.
```

The ordinal limit is not the scalar value `Z(1)`; it is the retained boundary
state / pole-order data: `r` serial geometric factors become `omega*r`.

## Exact shared boundary count

The tournament fixed-path `n->n+2` recursion adds two endpoint vertices and a
boundary shell. The number of newly exposed fixed-path tiling variables is

```text
binom(n+1,2) - binom(n-1,2) = 2n - 1.
```

This is exactly the same odd shell

```text
C = 2n - 1
```

used by THM-401/HYP-2161 as the LRC floor resonance observer. Thus the Go
`n+2` liberty recursion, the tournament `n->n+2` boundary correction, and the
LRC `C=2n-1` resonance are three faces of one boundary-partition function
phenomenon.

## What improves

This sharpens HYP-2157/HYP-2161/HYP-2181 and sits downstream of the published
HYP-2182/HYP-2183 partition-function and `n+2` recursion threads:

1. The master object is a partition function together with a boundary state,
   not a scalar count.
2. The `n->n+2` step should always be checked for the odd shell `2n-1`; this is
   not an LRC-only modulus but the boundary size of the tournament recursion.
3. Infinite games add a new diagnostic: ordinal value is a tropical/open-game
   shadow of the same partition function, where pole order or serial fuel
   depth matters more than `Z(1)`.
4. For LRC `n=14`, the familiar `C=27` is also the boundary shell of the
   `14->16` tournament recursion, not merely a residue modulus. This makes
   lift/CRT conservativity a boundary-state theorem.

## Assumption challenge

The natural vertices are not Go stones, tournament arcs, or runners. The useful
vertices are:

```text
finite-announcement fuel cells,
boundary tiles,
odd-cycle packets,
depth cells p_k,
LRC shell/lift states,
proof obligations.
```

The quotient must preserve the target predicate: life/death of a marked stone,
Hamiltonian-path count, `p_0=0`, or a proof obstruction. It may forget the
physical board geometry, the particular tournament edge, or the raw speed list
only when it retains the boundary partition function deciding the predicate.

## Artifacts

- `04-computation/infinite_go_partition_recursion_s616.py`
- `05-knowledge/results/infinite_go_partition_recursion_s616.out`
- `07-reflections/infinite-go-nplus2-partition-functions-s616.md`
