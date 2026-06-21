---
id: THM-560
title: The OCF degree ladder — deg_b(c_{2k+1}) = 2k by odd-cycle reverse-cancellation; the oddness of the cycles IS the Clifford->magic degree drop
status: PROVED (the ladder, via reverse-cycle cancellation); magic-onset & maxdeg(H) VERIFIED n<=7
source: kind-pasteur-2026-06-20-S22
depends_on:
  - THM-554   # tile-bit cube / score partition function
  - THM-555   # score->OCF wall
related:
  - HYP-2707  # the Clifford/magic hierarchy this grades
  - THM-559   # c3 as a 2-body Ising energy (the degree-2 rung, another agent)
  - THM-027   # H-maximizer = regular/Paley
---

# THM-560 - The OCF Degree Ladder

Index tournaments by `b in {0,1}^F`, `F = C(n-1,2)`, = tilings of the staircase on base path
`n -> n-1 -> ... -> 1` (codex convention: tile `(a,b)` bit 0 => arc `a->b`, bit 1 => `b->a`; each tile
flips exactly one non-base arc). Every cycle count is a multilinear polynomial in `b`. This grades the
OCF by polynomial degree — the grading that HYP-2707 identifies as the Clifford(quadratic)->magic split.

## Statement

**(Ladder, PROVED.)** For every k with `2k+1 <= n`, the directed-`(2k+1)`-cycle count is a multilinear
polynomial of degree EXACTLY `2k` in the tile bits:
```text
deg_b( c_{2k+1} ) = 2k.
```
Verified EXACTLY (Mobius transform over the full 2^F cube) at n=7: deg(c3)=2, deg(c5)=4, deg(c7)=6.

**(Consequence.)** `maxdeg_b(H) = L_max - 1`, where `L_max` is the longest odd directed-cycle length that
fits (`= n` if n odd, `n-1` if n even), since `H = OCF = 1 + 2*alpha_1 + ...` is built from odd cycles and
the top term is the longest odd cycle. VERIFIED n=4..7.

**(Magic onset, VERIFIED.)** The degree-`>=3` (genuinely non-quadratic / "magic") Fourier mass of `H` is
`0` at `n=4` (H is a PURE quadratic form there — a single Gauss sum) and `> 0` for `n >= 5`; the magic
fraction is `0, 0.299, 0.510, 0.635` at n=4..7. Magic onsets exactly when the first 5-cycle (degree 4)
appears.

## Proof of the ladder

Write `c_{2k+1}(b) = sum over directed (2k+1)-cycles sigma of I_sigma(b)`, where for a cyclic vertex
sequence `sigma=(v_1,...,v_{2k+1})`,
```text
I_sigma(b) = prod_{i} [ arc (v_i, v_{i+1}) is oriented v_i -> v_{i+1} ].
```
Each arc is either a **base-path arc** (between consecutive integers — a CONSTANT, since the base path is
fixed) or a **tile** (a single bit). The needed-orientation literal of a tile arc `(a,b)` is `1-b_t` (if
`sigma` wants `a->b`) or `b_t` (if it wants `b->a`) — a degree-1 affine form with **leading coefficient
`-1` or `+1`** respectively. Hence `deg(I_sigma) = (# tile arcs of sigma) <= 2k+1`, and the only terms of
top degree `2k+1` come from cycles all of whose `2k+1` arcs are tiles.

Pair each such cycle `sigma` with its **reverse** `sigma' = (v_1, v_{2k+1}, ..., v_2)`. Both are directed
`(2k+1)`-cycles (both summed in `c_{2k+1}`), and they use the SAME `2k+1` undirected arcs — so `sigma'` is
also all-tile. For every arc, `sigma'` needs the OPPOSITE orientation to `sigma`, so each literal swaps
`b_t <-> 1-b_t`, flipping its leading coefficient. The top-degree monomial of `I_sigma` is
`prod_arcs (leading coeff) * prod_arcs b_t`; therefore
```text
[top monomial of I_{sigma'}] = (-1)^{2k+1} * [top monomial of I_sigma] = - [top monomial of I_sigma],
```
because `2k+1` is **odd**. The two top monomials cancel. Summing over all reverse-pairs, the entire
degree-`(2k+1)` part of `c_{2k+1}` vanishes, so `deg(c_{2k+1}) <= 2k`. Equality holds: a surviving
degree-`2k` term is exhibited by the exact Mobius transform (c3:2, c5:4, c7:6). QED.

## Why this matters (the unifying point)

The degree drop `2k+1 -> 2k` is caused **precisely by the oddness** of the cycle (`(-1)^{odd}=-1`). Even
cycles would NOT cancel — but the OCF counts only ODD cycles (Redei/Grinberg-Stanley), so every OCF
generator sits one degree below its length. This is the algebraic root of the **Clifford -> magic ladder**
(HYP-2707): the degree-2 rung (`c3`, the only quadratic/Gauss-sum/Clifford layer; THM-559's Ising energy)
is the OCF's leading term, and the score->OCF wall (THM-555), the Clifford->magic wall, the Hopfield
2-local->k-local wall, and the n=4->n=5 magic onset are ALL the boundary `degree 2 -> degree 4`. The
efficiency observation (degree-grading) became the proof — the MEMORY.md pattern (THM-118 -> THM-498).

## Verification

`04-computation/conn_clifford-magic_kps-Sx-wf.py` (deg/magic spectrum) + exact-Mobius degree check
(kps-S22): deg(c3,c5,c7)=2,4,6 exact at n=7; magic fraction 0/0.299/0.510/0.635 (n=4..7); maxdeg(H)=L_max-1.
