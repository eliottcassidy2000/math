# A pair of pants inside the four-vertex tournament

2026-11-27. The pants analogue is literal after deleting the closing
edge, and its three reversals are exact on a marked tiling. The result
also explains the first failure of the Hamiltonian-count quotient to
remember an edit response.

## 1. Why three cuffs appear

Fix the directed Hamiltonian path \(1\to2\to3\to4\). Use bits
\((a,b,c)\) for the chords \(13,14,24\), with 0 forward and 1 backward
relative to the path. The fixed path closes if and only if \(b=1\).

Delete the closing edge \(14\) temporarily. The underlying graph
\(K_4-14\) is a theta graph: between 2 and 3 there are the direct edge
and the two routes through 1 and 4. Thicken its planar embedding to
an orientable surface with boundary. It retracts onto a graph with
four vertices and five edges, so has Euler characteristic \(-1\).
It is planar, hence has three boundary circles: a pair of pants.
The two leg cuffs run along \(123\) and \(234\); the waist runs along
\(1243\). The two bounded facial cycles can be consistently signed
so that their sum is the exterior cycle, with the common edge \(23\)
canceling. Equivalently this is the boundary relation in homology.

The leg cuffs run as directed cycles exactly when \(a=1\) and \(c=1\)
respectively. The waist cannot be directed in either direction:
its fixed edges \(12\) and \(34\) run in opposite senses around it.
Reinserting \(14\) creates an additional link across this surface
model. Its orientation is not recoverable from the two leg-cycle
flags, and is exactly the distinguished path's closing bit.

## 2. The three nonidentity reversals

A hyperbolic pair of pants has three seams, cutting it into two
congruent right-angled hexagons; see
[Petri, §7.5](https://webusers.imj-prg.fr/~bram.petri/t_tt/notes_250127.pdf).
Mark a waist, order the legs, and distinguish the two sheets.
The two waist-to-leg seams correspond to \(13,24\), and the
leg-to-leg seam corresponds to the closing coordinate \(14\).
The marked three-bit seam data record a sign relative to a choice
of sheet. In this model:

| Garment operation | Tiling operation | Action on \((a,b,c)\) |
| --- | --- | --- |
| Inside out | Complement all free tiles, \(\phi\) | \((1-a,1-b,1-c)\) |
| Backwards | Reflect the path axis, \(\tau\) | \((c,b,a)\) |
| Both | \(\tau\phi=\phi\tau\) | \((1-c,1-b,1-a)\) |

The two involutions commute and generate the Klein four group.
This is a model with markings: the unmarked topological surface
does not pick a waist, a sheet, or a distinguished involution
that swaps specified cuffs. On the hyperbolic geometry the
sheet exchange is the seam reflection; the leg exchange is a
length-preserving symmetry when the two leg lengths are equal,
or a change of marking otherwise. The correspondence of seams
to tiles is separate from the cuff cycles in the theta spine.
In particular, \(\phi\) complements the *three free* arcs, not
all six. On the normalized path, \(\tau\) is obtained by taking
the converse of the tournament and then renaming \(i\mapsto5-i\).
Thus \(\tau\) preserves both path and cycle counts, and fixes
the designated closing bit. Neither of the sheet-changing
operations fixes any tiling, because both toggle \(b\).

### A second, intrinsic four-point model

Over \(\mathbb F_2\) the three cuffs give the three nonzero classes
of \(H_1(P;\mathbb F_2)\cong\mathbb F_2^2\), with their sum zero.
Label four vertices by the four points of this affine plane.
For each nonzero cuff class \(v\), translation by \(v\) pairs up
the four vertices into a perfect matching. Thus the cuffs line up
with the three pairs of opposite edges of the \(K_4\) tetrahedron.
Permuting four vertices permutes these three matchings; its kernel
is exactly the three double transpositions plus the identity,
so \(S_4\cong\mathbb F_2^2\rtimes S_3\) in this model.
This translation Klein group is a separate realization, not an
identification of the tile complement with a vertex permutation.

## 3. What the three reversals actually do

Let \(H\) count directed Hamilton paths, \(C\) directed Hamilton cycles
up to choice of start, and \(N_3\) directed triangles. Enumerating the
eight tilings gives:

| Bits \(abc\) | \(H\) | \(C\) | Sorted outdegrees |
| --- | ---: | ---: | --- |
| 000 | 1 | 0 | 0,1,2,3 |
| 001 | 3 | 0 | 1,1,1,3 |
| 010 | 5 | 1 | 1,1,2,2 |
| 011 | 5 | 1 | 1,1,2,2 |
| 100 | 3 | 0 | 0,2,2,2 |
| 101 | 5 | 1 | 1,1,2,2 |
| 110 | 5 | 1 | 1,1,2,2 |
| 111 | 5 | 1 | 1,1,2,2 |

The two \(H=3\) classes are converses, not isomorphic. Each consists
of a directed triangle and an added source or sink. The five
\(H=5\) tilings form one tournament class: the only strong
four-vertex tournament up to isomorphism. For example it consists
of a cyclic triangle with a fourth vertex winning once or twice
against the triangle; cyclic rotations identify placements and
there is an isomorphism between these two descriptions by
choosing a different cyclic triangle. The companion replay also
certifies the isomorphism classification directly.

The three Klein-four orbits are
\[
 \{000,111\},\qquad \{010,101\},\qquad \{001,100,110,011\}.
\]
Inside-out pairs are \(000/111,\ 010/101\) (blue) and
\(001/110,\ 100/011\) (black). Hence the converse-merged
metagraph is \(1\)--\(5\)--\(6\) in tiling masses, with one
blue edge from 1 to 5, a blue loop at 5, and two black edges
from 5 to 6. The mass 6 is the merged \(H=3\) converse pair,
not a single tournament with six paths.

Most notably, \(\phi\) does **not** descend to a function even
on tournament isomorphism classes. For one and the same
five-path class, its outputs can be transitive, a three-path
class, or the five-path class again. For example the three
views \(111,110,101\) give \(000,001,010\) respectively.
This is a four-vertex model of what equal path counts hide:
the way the path sits in the tournament.

## 4. A closing-edge response with explicit guards

The two new triangles involving the restored edge \(14\) close
exactly when \(b(1-a)=1\) and \(b(1-c)=1\). Hence
\[
 N_3=a+c+b(2-a-c),\qquad H=1+2N_3
\]
(the second identity follows directly from the eight-row computation).
The two pants legs are independent in the mixed-difference
sense: there is no \(ac\) term. Both interact with the closing
tile, and each interaction has mixed difference \(-2\).
Flipping \(b:0\to1\) changes \(H\) by \(4-2(a+c)\), or
\(4,2,0\) as zero, one, or two cuffs already cycle. In the neutral
case \(101\to111\), one physical arc changes whether the **marked**
path closes, while both the total path and cycle counts stay
fixed. This is the smallest concrete instance of the distinction
between marked closure and total counts explored in the
[polygon and closing-tile session](../closing_tile_polygon_dissections_and_metagraph_20261127/README.md).

Run `node replay.js` for exact checks.
