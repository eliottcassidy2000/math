# Cake, bagel, Moser and Fibonacci are one Pascal triangle — rows vs shallow diagonals

*boxeph-2026-07-21-S207. Owner: relate the repo's polygonal/polyhedral (figurate) work to Fibonacci and
to the cake & bagel cutting sequences. Builds on opus-S317 (Vandermonde truncation law, polygonal vs
polyhedral), klein-S313 ((r,g) shadow lattice, g-bonacci kernels, deficit-1), mac-mini-S137 (the Hurwitz
golden-corner principle). All verified exactly.*

## The one object, three readings

Everything is the **Pascal (polyhedral/simplex) triangle** `C(i+j−1,j)`, read three ways:

- **ROW sums (full)** = `2ⁿ`.
- **SHALLOW-DIAGONAL (skip) sums** = **Fibonacci** `1,1,2,3,5,8,13,…` (verified: `Σ_k C(n−k,k)`).
- **TRUNCATED-ROW (fixed-depth) sums** = the **figurate cutting sequences** — this is where cake and
  bagel live.

## Cake, bagel, caterer, Moser: truncated-Pascal row sums (verified)

| cutting object | sequence | closed form |
|---|---|---|
| lazy caterer (2D disk, `n` lines) | `1,2,4,7,11,16,22,…` (A000124) | `C(n,0)+C(n,1)+C(n,2)` (Pascal to depth 2) |
| **cake** (3D ball, `n` planes) | `1,2,4,8,15,26,42,64,93,…` (A000125) | `C(n,0)+C(n,1)+C(n,2)+C(n,3)` (Pascal to depth 3) |
| Moser circle (2D disk, `n` chords) | `1,2,4,8,16,31,57,99,163,…` (A000127) | `C(n,0)+C(n,2)+C(n,4)` (the **polygonal** row-sum, opus-S317) |
| **bagel** (solid torus, `n` planes) | `1,2,6,13,24,40,62,…` (3 cuts → 13) | `C(n,3)+n(n+1)` |

The **cake** is a clean depth-3 truncation of `2ⁿ` (drop `C(n,≥4)` — a ball has no room for the higher
intersections). The Moser circle is the repo's *polygonal* triangle's row sum, which opus-S317 proved is
**the first two Vandermonde layers of the polyhedral (Pascal) triangle** — the polygonal/polyhedral gap
is a self-similar copy of Pascal two layers deep.

## The bagel's hole is the deficit-1 (the surprise)

`bagel(n) = C(n,3) + n(n+1)` — the ball's cubic term plus **twice the triangular number** (the extra
regions the torus's handle buys). Against the cake:

> **`bagel(n) − cake(n) = T_n − 1`** (verified n≥1): the torus exceeds the ball by exactly the
> **triangular numbers minus one**.

That `−1` is not incidental. It is the **same deficit-1** that klein-S313's `(r,g)` **shadow lattice**
missing-region law carries ("the deficit-1 is the g-bonacci kernel's boundary effect"). The topological
hole of the bagel and the boundary effect of the g-bonacci kernel are the *same* off-by-one — a genuine
bridge between the cutting geometry and the Fibonacci-kernel side.

## Fibonacci enters as the diagonal — and the g-bonacci kernels tie it back

The **same** triangle whose truncated rows give cake/bagel/Moser has **Fibonacci as its shallow-diagonal
sum**. And the repo's native machinery is the **g-bonacci kernel** `1/(1 − x − x^{g+1})` (klein-S313):
- `g=1`: `1,1,2,3,5,8,13,…` — **Fibonacci** exactly (verified);
- `g=2,3,…`: `1,1,1,2,3,4,6,9,…`, `1,1,1,1,2,3,4,5,7,…` — the gap-`g` shadow-lattice family.

So the cutting sequences (rows) and the Fibonacci/g-bonacci sequences (diagonals) are **two projections
of one triangle**, and the g-bonacci kernel is the generating-function bridge. This is exactly
mac-mini-S137's *Hurwitz golden-corner principle* seen figurately: the row direction counts geometry
(pieces of cake/bagel), the diagonal direction counts the golden/Fibonacci obstruction data — and both
JC₂ (the golden-degree corner) and LRC(14) (the anti-golden Eisenstein extremal, the penultimate
convergent it forbids) live on that same Pascal/Fibonacci scaffold.

## The map

```
                          PASCAL (polyhedral) triangle  C(i+j-1, j)
        full rows → 2ⁿ        truncated rows → figurate cutting        shallow diagonals → Fibonacci
                                 depth-2 = lazy caterer                       g-bonacci kernel
                                 depth-3 = CAKE (ball)                        1/(1 - x - x^{g+1})
                                 C(n,3)+n(n+1) = BAGEL (torus)                 g=1 = Fibonacci
        polygonal triangle → Moser circle A000127 (opus-S317, 2 Vandermonde layers)
        bagel - cake = T_n - 1  =  the DEFICIT-1  =  klein-S313 shadow-lattice missing-region boundary
```

## Takeaway

The cake and bagel are not exotic — they are the depth-3 and torus-corrected truncations of the same
Pascal triangle whose shallow diagonals are Fibonacci. The repo already holds the unifying theory: the
**polygonal-vs-polyhedral Vandermonde-truncation law** (opus-S317) places the cutting sequences on the
figurate side, the **g-bonacci kernels** (klein-S313) place Fibonacci on the diagonal side, and the
**deficit-1** — the bagel's hole — is literally the shadow-lattice boundary effect. So "cake, bagel,
Moser, Fibonacci" is one Pascal triangle read row-wise vs diagonal-wise, and it is the same golden/
figurate scaffold on which JC₂ and LRC(14) sit (mac-mini-S137).

Links: HYP-8820,
[[the-vandermonde-truncation-law-polygonal-vs-polyhedral-opus-S317]],
[[the-hurwitz-principle-jc2-golden-corner-lrc-ghost-convergents-gbonacci-macmini-S137]],
[[what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206]].
