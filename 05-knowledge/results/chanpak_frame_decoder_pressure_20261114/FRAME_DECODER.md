# Decode a locked frame and the whole count jointly

**Scope.** This is a target-first decoder for arbitrary **fixed** two-hub
templates, not just the two special rails of
[Session 336](https://github.com/eliottcassidy2000/math/blob/8d2d70178f9d17b191b93dcaffa6c1a8e25dbcb6/05-knowledge/results/session336_overlap_inverse_and_intercyclic_20261114/TARGET_RAILS.md).
It specifies which unbounded variables can be removed, inverts a chosen frame
and ordinary count simultaneously, and returns a single isomorphism class
within any given finite bank. This is not a statement that a finite bank, or
even all two-hub graphs, attains every permitted integer.

## 1. The nonnegative coordinates of a block

In a two-hub word (hubs $x\to y$), choose $d$ disjoint marked homogeneous
contiguous blocks of widths $k_1,\ldots,k_d\ge1$, with other sites and their
digits fixed. After deleting the hubs, each Hamiltonian path has three ordered
slots. Inside each slot the spine is ascending. Fix:

* the hub order;
* the slot of every unmarked site;
* the nonempty set of slots $S_j$ occupied by block $j$.

This specification is either legal for every realization of those slot sets,
or illegal for every realization: hub incidences depend on the letters at
extreme slot sites, a whole marked block has the same letter, and the positions
of all blocks and unmarked sites are fixed. Empty slots and the order of
consecutive hubs are also determined before assigning individual sites.

Let $c_{s_1,\ldots,s_d}$ count the legal specifications having
$|S_j|=s_j$. Then
$$
 H(k_1,\ldots,k_d)=\sum_{\boldsymbol s\in\{1,2,3\}^d}
   c_{\boldsymbol s}\prod_j U_{s_j}(k_j),\qquad
 U_1(k)=1, U_2(k)=2^k-2, U_3(k)=3^k-3\,2^k+3 .
$$
The $U_s(k)$ are surjections onto $s$ labelled slots, or $s!S(k,s)$.
They make visible nonnegative coordinates that can be obscured by signed
coefficients in the usual exponential basis. In one direction the exact
inverse of the evaluation map at widths 1,2,3 is
$$
  c_1=H(1),\qquad c_2=(H(2)-H(1))/2,\qquad
  c_3=(H(3)-3H(2)+2H(1))/6 .
$$
Tensor this inverse in all $d$ directions. The resulting integer tensor can
be compiled from $3^d$ calls to the inherited exact 36-state counter. Width
zero is a different skeleton (a block can no longer occupy any slot); do not
evaluate a positive-width polynomial at zero. The combinatorial expansion
also holds for any fixed number $h$ of hubs, with $h+1$ slots,
surjections $s!S(k,s)$, and the evaluation grid $1,\ldots,h+1$.

## 2. Detecting every true width variable

Test a marked segment at width 1, with all other variable widths 1. If its
vertex lies in an oriented triangle, duplication makes another triangle with
the same other two vertices. The
[Grinberg–Stanley odd-cycle packing expansion](https://arxiv.org/abs/2307.05569)
includes every old packing unchanged, plus the new singleton-cycle packing
which contributes 2. Hence the coefficient with $s_j=2$ and all other
indices 1 is at least 1. This argument excludes a live segment with no
two-slot coefficient at the all-one base.

Conversely a triangle touches each duplicated block in at most one vertex,
since any two twins have identical arcs to any outside vertex and have a
transitive mutual orientation. Any triangle at any widths through the segment
would collapse to a triangle in the all-one skeleton. Thus if the base test
finds none, every vertex of that segment remains in no directed triangle and
is a singleton strongly connected component. (A nonsingleton strongly
connected vertex lies on a directed triangle: on a shortest directed cycle
through the vertex, a chord from it to the third vertex or from the third
vertex to it either shortens the cycle or completes a triangle.) Because
paths multiply along the ordered components, the segment is truly count-neutral
for every choice of the other widths. The decoder normalizes it to width one,
meaning it returns every active solution, not infinitely many inert variants.

In a legitimate proper frame on four marked fixed-width-one sites $i<j<k<l$,
the two backwards hub bridges overlap, and every vertex between $i$ and
$l$ is in the same strong component as the hubs. In fact between the low
and high anchors $a,b$ of a bridge, a site $v$ lies on
$a\to v\to b\to\text{hub}\to a$, and the two bridge intervals intersect.
Thus no width counted by $\alpha$ or $\beta$ is inert.

Every live $k_j$ contributes the positive term $U_2(k_j)$, even when all
others are minimal; the other terms are nondecreasing. Consequently at a
target $N$, $k_j\le L=\lfloor\log_2(N+2)\rfloor$, and $H$ increases
strictly with that width under any fixed choice of the others.

## 3. Invert marked and whole targets together

The two pairs of extreme anchors supply independent two-by-two permanents
$A,B$, each counting distinct assignments of its two anchors to two hubs.
Zero permanent rejects the frame. With $\alpha$ the number of sites in the
two outer wings and $\beta$ the number in the central overlap, the entire
contribution of **that quadruple**, summing both hub orders and all its seam
matchings, is
$$
 M=2AB\,2^\alpha3^\beta .
$$
This is never substituted for the whole count. Given both $N$ and $M$,
divide by $2AB$; the remaining factor must be 2,3-smooth, whose valuations
give the required $\alpha,\beta$. Each is a disjoint linear equation in
the marked widths, plus the fixed sites in its zones. Eliminate one width
per nonempty graded group, rejecting impossible grades. If $r$ live
widths remain after those eliminations, enumerate $r-1$ on $[1,L]$
and binary-search the remaining width using ordinary $H$; absent free
widths need just one ordinary check. The binary-searched width is in an
unconstrained zone when a marked target is active, so no eliminated value
moves against it. With fixed template the search uses
$O(L^{r-1}\log L)$ exact tensor evaluations for $r\ge1$.
Invalid seam labels, $M>N$, bad factorizations, grades, and final
absolute-count mismatches each reject independently.

For example:

| exact whole $N$ | marked $M$ | template | anchored solution |
|---|---:|---|---|
| 385 | 36 | $1\,0^a\,30\,2^{b+1}\,1$ | $a=b=1$, 1030221 |
| 573 | 108 | $2\,0^a\,30\,2^{b+1}\,12$ | $a=b=1$, 20302212 |

A deterministic class selector can choose minimum surviving order followed
by the least adjacency code among all lexicographic refinements of the
degree classes, which is an isomorphism invariant. Search marks are not
part of the output class. This selection certifies **one** class for every
target the decoder hits; it does not turn a failed finite bank search into
a proof of global nonattainment. In passing from a fixed number of blocks to
unbounded shapes, the tensor is local and reusable, but one still needs a
coverage or descent theorem selecting the shape.

Run `python replay_frames.py` to verify the decoder, examples, substitutions
in different block shapes, and the 63 obstruction below. See
[RESIDUAL_SELECTION.md](RESIDUAL_SELECTION.md) for the weaker scalar
selector that fails and an exact residual-based notion that does compose.
