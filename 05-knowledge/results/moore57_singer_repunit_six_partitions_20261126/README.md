# Repunit 57, an orthogonal array on 54 labels, and the six-partition Moore obstruction

2026-11-26. Continues the [two sharp prism shells and boundary placements](../moore57_sharp_prism_boundary_20261126/README.md), using the archive's [moving-base repunit constructions](../session218_smith_repunit_totient_20261107/smith/SMITH_REPUNIT_BRIDGES.md), [Smith-base completion versus congruence testing](../session218_smith_repunit_totient_20261107/smith/BASE_COMPLETENESS.md), and [THM-4940's distinction between integer targets and base-dependent odd-cycle reads](../../../01-canon/theorems/THM-4940-odd-cycle-packing-counts-in-every-base.md). The strongest actual six-partition witness below fits **2,477 of the 2,917 vertices**, including the three distinguished blocks. It is neither a full incidence design nor a Moore graph.

## 1. An exact repunit construction for all six unconstrained columns

Write $R_b(\ell)=1+b+\cdots+b^{\ell-1}$. In the moving-base construction every odd prime $p$ is $R_{p-1}(2)$, but an equality of repunit values is not itself a common coordinate system. Here $R_7(3)=57=3\cdot19$ is composite for the familiar cyclotomic reason: $7\equiv1\bmod3$, leaving the other prime $19\equiv1\bmod3$. What matters is the geometry of the *base*: $R_7(3)=(7^3-1)/(7-1)$ counts the projective points over $\mathbb F_7$. This is the same length-three polynomial that gives 7 in base 2 and 21 in base 4, but the Chan–Pak holes concern binary tournament path counts, not existence of this geometry.

**Proved.** Delete any three collinear points from $\mathrm{PG}(2,7)$. On the 54 survivors the remaining lines are one block of size 5, 21 blocks of size 7 and 35 of size 8. The chosen line loses all three points; for each removed point the other seven lines lose exactly one; the rest lose none. Every pair of surviving points still lies in exactly one block.

This gives a literal strength-two orthogonal array on 54 labels in **six** columns:

1. At each size-7 or size-8 block use the local field array $a+s_i t$ with six distinct $s_i$ and $t\ne0$. The excluded $t=0$ rows are exactly the diagonal rows. Thus there is one row for every *distinct* ordered pair of symbols in every pair of columns, and none for equal ones.
2. On the unique size-5 block use all 25 rows from the five slopes and the parallel direction $t$. This covers both diagonal and off-diagonal pairs in that block. A six-column field array of order 5 cannot reserve a shared diagonal: any two distinct rows in a six-column strength-two array of order 5 agree somewhere (for each fixed row, the six columns already account for all 24 other rows).
3. Add one constant row for each of the 49 labels outside the size-5 block. Since every unequal pair has a unique block, there are exactly 2,916 rows and all 43,740 oriented column-pair cells once.

This is the positive bridge from repunits: the three collinear deletions give exactly the three block sizes that support this patch. The replay uses explicit homogeneous coordinates and polynomial arithmetic over $\mathbb F_8$. It does not assume an order-54 field.

The connection is **not** an inherited cyclic symmetry. In [Singer's projective indexing](https://doi.org/10.1090/S0002-9947-1938-1501951-4) every line is a difference set with each nonzero ordered difference used once. Therefore no nonidentity Singer translation fixes a subset of a line of size at least two: for each element of a stable subset, applying that shift would repeat the same ordered difference. In particular the three deleted points, and their complement of size 54, have trivial stabilizer in the Singer cycle. There is also no order-53 rotation fixing the six parent labels on a 54-sheet prism shell: a nontrivial order-53 permutation fixes one sheet in a fiber, and the matchings propagate these unique fixed sheets to every fiber, giving a literal triangle over a triangular face.

## 2. Two exact restrictions on using an ambient orthogonal array

### For the four fixed boundary inputs: deletion alone cannot work

For each pair of columns put an edge on a cell **not** allowed for a new ordinary tuple. This includes both holes from the shell and ordinary pairs already used by the distinguished-block witnesses. The resultant graph has 2,430 forbidden pairs for each concurrent case and 2,415 for each separated case. In any ambient 2,916-row orthogonal array, a row containing a forbidden pair has to be deleted. Each forbidden pair must be accounted for by exactly one deleted row, but the rows cannot concentrate the forbidden pairs enough.

The replay exhausts every six-tuple containing at least 13 forbidden pairs: such a tuple has a four-clique in the forbidden graph (the six-vertex Turán bound for a $K_4$-free graph is 12). Enumerate all four-cliques from four different columns, fill the last two columns, and deduplicate. If there are $h_k$ high tuples of weight $k$, then $m$ deleted rows can account for at most $12m+\sum_{k>12}(k-12)h_k$ forbidden pairs. This upper bound even allows simultaneous use of high tuples sharing a pair, so it is conservative.

| Fixed input | 13, 14, 15 forbidden pairs in a tuple: counts | Minimum OA rows deleted | Upper bound surviving | Needed ordinary rows |
|---|---|---:|---:|---:|
| Equal triangles, concurrent | 56, 5, 1 | 197 | 2,719 | 2,754 |
| Equal triangles, separated | 59, 6, 0 | 196 | 2,720 | 2,755 |
| Paired squares, concurrent | 33, 1, 0 | 200 | 2,716 | 2,754 |
| Paired squares, separated | 41, 7, 0 | 197 | 2,719 | 2,755 |

This excludes deletion-only constructions from **any** order-54 six-column array on the prescribed boards, not just the repunit array or particular labelings. Mixed trades replacing surviving tuples are not excluded.

### For flexible boundary inputs: three coherent removed layers cannot work

Suppose the concurrent boundary is built from three *distinct* transversal layers of an ambient array, one per distinguished block. Each layer consists of 54 six-tuples covering every symbol of each column once. Suppose deleting those three layers gives the exact ordinary design and each exclusive distinguished group is the projection of its assigned deleted layer to its four ordinary columns.

Label each layer by its distinguished opposite-parent group $g$. For opposite columns of group $g$, shell cells are precisely the direct matching, since the two other boundary layers consume the two other forbidden cells in each row. Thus each direct opposite edge occupies layer $g$. For a pair of columns from groups $g,h$, their shell cells occupy layers $g,h$, because the remaining layer supplies their boundary pairs. Pick a horizontal prism edge $x-y$ between two such groups. Both $x'-y$ and $x-y'$ are length-two shell cells, where $x',y'$ are the vertical mates. If $x-y$ were in layer $h$, either $x'-y$ is in that layer (forcing $x,x'$ in layer $h$, a contradiction), or $x'-y$ is in layer $g$ (forcing $x-y$ in layer $g$, a contradiction). The other endpoint symmetrically forbids layer $g$. Neither is available. This rules out globally coherent three-layer deletion even if both shell and boundary can change. It does not exclude a non-layered hole decomposition, or a design not coming from an ambient array.

## 3. A moment law for a regular abelian development

Assume a regular abelian sheet action on all six fibers, fixing parent indices. Its 54-element sheet group is written additively. A separated boundary is impossible under the action: the unique member of each intersection of two distinguished blocks would be fixed, but has ordinary coordinates in the third opposite group, on which the action is free. In the concurrent case each of the three exclusive boundary groups is one 54-tuple orbit.

Gauge the three vertical prism matchings and the first two top triangular edges to shift 0 (these five edges form a tree). In cyclic orders $(0,2,4)$ and $(3,5,1)$, write top shifts $(0,0,t)$ and bottom shifts $(r_0,r_1,t+r_2)$, so the bottom holonomy is $u=t+r_0+r_1+r_2$.

In any four-column circuit the sum of ordinary differences over *all* ordinary tuple orbits is zero by telescoping. Use the circuit $0,1,3,4,0$. Its boards are all between the same two opposite groups, so the boundary-orbit differences telescope on the circuit too. In each board the sum of every possible difference is the sum $\sigma$ of the abelian group; $4\sigma=0$ since negation pairs its elements, leaving only the elements killed by doubling. The sums of the two shell-hole shifts around this circuit are
\[
 -2t-r_2,\quad t+r_2-r_1-r_0,\quad -2t-r_2,\quad t.
\]
These too must sum to zero, so **$t+u=0$**. In particular, a regular abelian full incidence design cannot have all squares paired: an abelian group of order 54 has just one nonzero element of order 2, so three swap shifts sum to that element, while $r_0+r_1+r_2=-2t$ is in the odd-order doubling image.

This moment condition concerns hypothetical translation shells. The two pinned sharp splices are more rigid: an explicit centralizer test finds that their only sheet action fixing all six parent labels is identity. A fiber-preserving action in tree gauge must be the same permutation everywhere and commute with the four fundamental transports, so trying the image of a single sheet suffices by transitivity. The argument does not assume that a graph isomorphism must fix the hexagon parents.

## 4. Verified six-way placement, beyond the distinguished blocks

A best-of-250 relabeling run on the repunit array retains 1,619, 1,606, 1,636 and 1,631 ordinary rows respectively after deleting the rows with forbidden pairs. Because of the obstruction in §2 it is useful only as an initial reservoir. A separate search using direct six-column common-neighbor intersections and local exchanges yields the stored witnesses below. No search result is an upper bound on the true packing size.

| Same input order as §2 | Ordinary tuples fitted | Boundary tuples retained | Total fitted of 2,917 | Left per pair board |
|---|---:|---:|---:|---:|
| Equal triangles, concurrent | 2,312 | 163 | 2,475 | 442 |
| Equal triangles, separated | 2,308 | 162 | 2,470 | 447 |
| Paired squares, concurrent | 2,314 | 163 | **2,477** | **440** |
| Paired squares, separated | 2,304 | 162 | 2,466 | 451 |

The exact replay finds **zero** complete six-tuples using only the leftover pairs of each stored packing. This is local maximality, not optimality. Even the next single gain requires exchanging existing rows. An appropriate next target is an exchange that simultaneously frees a clique and leaves a second one, or a redesigned boundary whose surplus holes align beyond pair marginals.

## Replay

- Node.js **node boards.js --check** reconstructs all four exact inputs from the preceding repository packet; **node symmetry.js** checks the two fiber-preserving automorphism statements.
- Python 3 **python3 singer_seed.py** constructs the full repunit seed; pass **b0.txt b1.txt b2.txt b3.txt** to rerun the finite relabeling control.
- Python 3 **python3 highholes.py** exhausts the high-weight forbidden tuples and verifies the ambient-OA upper bounds; **python3 verify.py** checks admissibility, pair uniqueness, and local maximality of the stored six-way placements.
- **pack.cpp** is the heuristic continuation, not a completeness test. It takes a board, output path, RNG seed and iteration cap, optionally followed by a text seed (one six-tuple per line). Each character of a **witnessN.b64** file records one numerical coordinate 0–53 under the URL-safe base64 alphabet; these are plain text, not a compressed executable.
