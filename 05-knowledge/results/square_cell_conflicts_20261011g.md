# Exact conflicts in the two-scale square pattern

**Status: PROVED for a declared whole-cell, fixed-period model; FINITE-EXACT scale-cover certificates.** The optimum palette over all scales in this model is eight. This does not classify arbitrary colour periods, diagonal subdivisions, position-dependent rules, the Stechkin construction, or arbitrary plane colourings.

## 1. The model and its inheritance

Partition the plane into half-open macroblocks

$$[4ia,4(i+1)a)\times[4ja,4(j+1)a),\qquad (i,j)\in\mathbf Z^2,$$

where $a>0$. If $i+j$ is even, subdivide the block into four squares of side $2a$; if odd, subdivide it into sixteen squares of side $a$. Every resulting square has one constant colour. **The colour function, as well as the tile partition, is required to have periods $(8a,0)$ and $(0,8a)$.** This is the natural forty-cell motif, not a claim that the user's three-cycle diagonal prescription necessarily has this colour period.

The closest inherited mechanism is [glued chessboard rings and scaffolds](glued_chessboard_rings_scaffolds_20261006.md), especially its distinction between a lattice colour character and the actual edge graph. Its graph uses specified chess moves; the graph below uses existence of a Euclidean unit segment between two regions. These edge predicates must not be identified. The [complete boundary tables](boundary_palette_compression_20261011.md) likewise retain the original source and complete constraints rather than only local labels. The active concept board is tile ownership, actual region distance, scale, colour period, and seam refinement.

The canonical hostile here is a cell-centre test: centres can avoid distance one even when two interior points of the same colour do not. A second hostile forgets translated copies of a motif cell. The least-used useful sidecar is the **actual translation witnessing each conflict**, which can differ from edge to edge.

In units of $a$, enumerate the forty cells by macroblocks $(i,j)=(0,0),(0,1),(1,0),(1,1)$, and within each block by lexicographic lower-left coordinates $(u,v)$, where $u,v\in\{0,2\}$ for a large-cell block and $u,v\in\{0,1,2,3\}$ for a small-cell block. The source interval is $[0,8)^2$. All coordinates and side lengths are integral before rescaling.

## 2. Exact half-open rectangle conflicts

For two nondegenerate half-open rectangles $C,D$, their displacement set is an open rectangle

$$D-C=(\ell_x,u_x)\times(\ell_y,u_y).$$

For example, $[b,b+w)-[c,c+v)=(b-c-v,b-c+w)$. Thus define

$$L=\max(\ell_x,-u_x,0)^2+\max(\ell_y,-u_y,0)^2,$$
$$U=\max(\ell_x^2,u_x^2)+\max(\ell_y^2,u_y^2). \tag{1}$$

The upper squared-distance endpoint $U$ is never attained. A positive lower endpoint $L$ is never attained either. The norm's image is an interval because the displacement rectangle is connected; its infimum and supremum are exactly (1). Hence, for any positive squared target $S$,

$$\text{some }p\in C,q\in D\text{ satisfy }|p-q|^2=S
\quad\Longleftrightarrow\quad L<S<U. \tag{2}$$

The possible attainment of zero is irrelevant because $S>0$. No boundary-null-set argument replaces (2): a plane colouring must include every point, and the chosen half-open ownership is used literally.

For physical unit distance set

$$S=1/a^2.$$

For motif vertices $i,j$, place an edge when (2) holds between $C_i$ and $C_j+(8p,8q)$ for **some** integers $p,q$. A same-vertex edge is a self-conflict and forbids any colour for that cell. Proper colouring of this quotient graph is equivalent to a plane colouring with the declared period and constant-cell rule. The quotient is a constraint graph: its eight-cliques need not be realized by eight simultaneous points, or even by one consistent choice of eight cell copies.

For any fixed positive rational $S$, the API checks finitely many period shifts. If $R=\lceil\sqrt S\rceil$, it is enough to use $|p|,|q|\le\lfloor R/8\rfloor+2$: farther shifts have an axis separation exceeding $R$. The mathematical criterion (2) applies to every real $S$.

## 3. Complete scale exclusion outside a compact interval

If $S<8$, a single large square has squared diameter supremum eight and contains a forbidden pair. At $S=8$ the diameter is not attained, so this particular obstruction disappears. This is exactly the condition $2a\sqrt2\le1$, with equality permitted for half-open cells.

If $S>36$, a large square conflicts with one of its own period translates. Here is an all-scale proof rather than a finite search. Take the large cell $[0,2)^2$ and translate it by $(8p,8q)$ with $p\ge1$ and $0\le q\le p$. Its conflict interval is

$$I_{p,q}=\left((8p-2)^2+\max(8q-2,0)^2,\ (8p+2)^2+(8q+2)^2\right).$$

For successive $q<p$, the next lower endpoint minus the current upper endpoint is

$$64(q-p)+32<0.$$

The intervals for one $p$ therefore cover

$$\left((8p-2)^2,\ 2(8p+2)^2\right).$$

The next strip begins at $(8p+6)^2$, and

$$2(8p+2)^2-(8p+6)^2=64p^2-32p-28>0\qquad(p\ge1).$$

These strips overlap strictly, start at $36$, and extend without bound. Their union is $(36,\infty)$. The open endpoint $S=36$ must be handled separately: an axis translate has infimum distance exactly one there, but does not attain it.

Consequently no palette at all works in this fixed-period model unless

$$8\le S\le36,\qquad\text{equivalently }1/6\le a\le1/\sqrt8. \tag{3}$$

This small-scale obstruction comes from a fixed physical colour period shrinking with $a$. It says nothing about increasing the motif size as the cells shrink.

## 4. Five eight-cliques cover the complete remaining scale range

For $8\le S\le36$, only translations $p,q\in\{-1,0,1\}$ can contribute. Any other translation has one axis gap at least eight, exceeding $\sqrt S\le6$. All distance-window endpoints are integers. Therefore the graph is constant on each open interval $(j,j+1)$, with integer endpoints checked separately: twenty-nine endpoints and twenty-eight open strata, fifty-seven exact cases in total.

The following motif-cell sets are eight-cliques:

| Name | Cell indices |
|---|---|
| $C_0$ | $0,1,2,3,4,6,12,14$ |
| $C_1$ | $0,1,2,3,4,7,13,19$ |
| $C_2$ | $0,1,2,3,5,17,24,27$ |
| $C_3$ | $0,1,2,3,5,18,27,28$ |
| $C_4$ | $0,1,2,3,4,6,13,15$ |

Their complete real scale cover is

| Scale interval for $S$ | Clique |
|---|---|
| $[8,10)$ | $C_0$ |
| $[10,13)$ | $C_1$ |
| $[13,17)$ | $C_2$ |
| $[17,18)$ | $C_2$ |
| $[18,25]$ | $C_3$ |
| $(25,26]$ | $C_4$ |
| $(26,36]$ | $C_0$ |

For every one of the twenty-eight pairs in each row, **interval_clique_packets** retains one actual period shift and its exact integer bounds $L,U$, and verifies that the whole stated real interval lies in $(L,U)$. The two $C_2$ rows deliberately retain different shifts: for pairs $(5,17)$ and $(24,27)$ the useful windows are $(4,17)$ and $(16,37)$. Their overlap bridges the change. This is an exact finite certificate covering all real scales, not extrapolation from floating samples. The independent fifty-seven-stratum check reproduces the cover.

In particular every allowed scale needs at least eight colours. Four of the eight cells in the critical $C_0$ are the four large squares in $[0,4)^2$; the other four have lower-left corners $(0,4),(0,6),(2,4),(2,6)$ and side one. The distant pairs are constrained through translated copies. Ignoring those copies loses edges.

## 5. An eight-colour positive control, and what a seam cannot repair

At $a^2=1/8$, regard the small squares in each odd macroblock as groups of four, temporarily making all coarse squares have side $2a$. If such a coarse square has indices $(p,q)\in\mathbf Z^2$, give all its points colour

$$c(p,q)=4\bigl((\lfloor p/2\rfloor+\lfloor q/2\rfloor)\bmod2\bigr)
       +2(p\bmod2)+(q\bmod2). \tag{4}$$

This uses eight colours and has the required $8a$ coordinate periods. It is constant on every original small square. Each coarse square has diameter supremum one, unattained. If two distinct coarse squares have the same colour, their index differences are $(2r,2s)$ with $r+s$ even. If $r,s$ are odd, the distance infimum is at least $2a\sqrt2=1$; if both are even, some coordinate index difference is at least four and the gap is larger. Any equality in the infimum is unattained. Thus no same-colour unit pair exists.

Together with the obstruction, (4) proves that the minimum palette over all positive scales in this specific model is **exactly eight**. This construction is an exact boundary control, not a seven- or six-colour improvement.

**Bounded diagonal seam corollary.** Suppose a proposed six- or seven-colour rule agrees, outside a strip $|x-y|\le W$ in even one of its two complementary half-planes, with a whole-cell rule of the declared $8a$ period. Recolouring arbitrarily inside the strip cannot make it a valid unit-distance colouring.

For a fixed scale, the obstruction uses finitely many actual cells and period copies: one interior conflict, one repeated-cell conflict, or the finitely many pair lifts of one quotient clique. Translate that entire finite collection by a sufficiently large vector $(8aT,0)$ or $(-8aT,0)$. It lies completely beyond the strip, retains every metric conflict and every periodic cell address, and still gives the contradiction. The two half-planes may use independent periodic patterns; a contradiction inside one suffices.

This does **not** prove that splitting any particular eight actual cells suffices or is necessary: the quotient clique uses potentially different copies for different edges. Nor does it forbid changes of period, splitting or recolouring that reach arbitrarily far into the bulk. Those changes genuinely alter the model's constraints.

## 6. Colour-dependent distances and the next useful model

For a vector of forbidden distances $(d_0,\ldots,d_{k-1})$, a cell assigned colour $j$ uses threshold

$$S_j=d_j^2/a^2.$$

The same exact source-owned checker applies, now with a different conflict graph for each colour. **verify_period_colouring** keeps all these thresholds rather than silently replacing them by one. For example the vector $(1,1,1,1,1/2,1/2)$ gives $(8,8,8,8,2,2)$ at $a^2=1/8$. The ordinary unit-distance eight-clique proof is not, by itself, a classification of that different list-constraint problem.

Uniformly dilating the plane multiplies every $d_j$ and $a$ by the same factor, so it preserves all $S_j$ and all ratios $d_i/d_j$. A global scale change cannot turn distances one and one-half into a common forbidden distance. The historical Stechkin input is handled in the companion source review; its tiling is not identified with these checkerboard squares.

The useful next choices are explicit: enlarge the colour period (a genuine three-cycle rule can require this), split cells with a declared boundary owner, or retain position-dependent phases. Every choice must recompute actual region conflicts. Centre-distance checks, diagonal colour names, and an almost-everywhere colouring do not replace the exact pointwise test.

| Source | Map / preserved predicate | Lost information and retained sidecar | Hostile |
|---|---|---|---|
| Half-open cells | Exact displacement windows; unit-pair existence | Centres lose the region extent; retain all sides and endpoint ownership | Adjacent large-cell centres miss a real conflict |
| Infinite periodic rule | Forty-vertex quotient | A motif-only graph loses period copies; retain each witnessing translation | Cells $0,6$ conflict only after a suitable wrap at $S=8$ |
| Variable scale | Seven exact interval certificates | A sample alone loses critical endpoints | Safe diameter and separation equalities at $S=8,36$ |
| Diagonal seam repair | Translate the finite obstruction into the bulk | Changes restricted to a bounded strip cannot alter its constraints | Pairwise lifts, not one fictitious eight-point unit clique |
| Colour-dependent distances | Threshold vector $(d_j^2/a^2)_j$ | One scalar scale loses their unequal ratios | $(1,1,1,1,1/2,1/2)$ is not ordinary unit-distance colouring |

## 7. Reproduction

[Exact script](../../04-computation/experiments/square_cell_conflicts_20261011g.py) and [saved output](square_cell_conflicts_20261011g.out):

```text
python -B -X utf8 04-computation/experiments/square_cell_conflicts_20261011g.py
python -B -O -X utf8 04-computation/experiments/square_cell_conflicts_20261011g.py
```

The finite controls include complete forty-cell ownership; all fifty-seven compact scale strata; seven interval certificates with all pair shifts; the explicit eight-colouring and its infinite-period finite reduction; 1,891 repeated-cell witnesses at thresholds $109/3$ through $1999/3$; strip-overlap identities through $p=100$; finite obstruction translations beyond a diagonal strip; thirteen type/domain/false-colouring hostiles; and exact colour-distance dilation controls. There are **10,625 always-active checks**. The unbounded-scale and bounded-seam statements are proved above, separately from those finite controls.
