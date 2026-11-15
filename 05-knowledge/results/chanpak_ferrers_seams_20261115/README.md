# Quantitative Ferrers rigidity of small-excess Busch seams

**Stronger follow-up:** [exact Ferrers rigidity and the Douglas sorting theorem](../chanpak_exact_seams_douglas_20261115/README.md) show the cross matrix is already Ferrers at surplus below n−2, and rule out the first positive even excess for strong tournaments at orders n≥8.

**All-order result.** Split a tournament into induced strong factors A and B meeting at v, and put L = H(T) − H(A)H(B). Choose *any* directed Hamiltonian cycle of each factor and start both at v. In the resulting cycle order the cross-arc matrix is within
`L + floor(L²/4)` entry flips of a Ferrers matrix (ones increasing along rows and decreasing down columns). Thus for excess 2 the distance is at most 3, and for excess 4 it is at most 8. No uniqueness of the factor cycles is required. This applies to every split in the two stable cases of the [Busch stability theorem](../chanpak_busch_stability_20261115/README.md).

There is also a linear-in-matrix-size collision test and an all-order exclusion for single flips away from a canonical one-way floor seam. Neither result excludes every small-excess mixed seam or completes Chan–Pak.

## 1. A two-sided lattice bound

Write the factor cycles as `v,a_1,...,a_p,v` and `v,b_1,...,b_q,v`, and D[i,j]=1 if a_i beats b_j. Let `X[i,j],Y[i,j]` count valid shuffles of the two initial lists with counts i,j, ending on the A and B sides respectively. Start with X[i,0]=1 for i>0, Y[0,j]=1 for j>0 and the opposite counts zero. For i,j>0,
```
X[i,j] = X[i-1,j] + (1-D[i,j])Y[i-1,j]
Y[i,j] = Y[i,j-1] + D[i,j]X[i,j-1].
```
Put F=X+Y with F[0,0]=1. Every pair of lists has a shuffle: choose the winning head of the two lists repeatedly, and append the remainder if one list runs out. Each chosen head beats both heads available at the next step. Thus F≥1, and the recurrence uses O(pq) integer operations. X is nondecreasing down columns and Y along rows.

For the reverse calculation in the converse tournament swap the sides. Its matrix is `R[i',j']=D[p+1-j',q+1-i']`. Write G for its total shuffle table. A pair of rotations of the two cycles determines precisely the i A-vertices and j B-vertices following v. The shuffles following v and those before v can be chosen independently; v beats the first vertex of either nonempty initial list and is beaten by the last vertex of either nonempty terminal list. Hence the fiber over these two rotations has exactly `F[i,j]G[q-j,p-i]` paths. All other pairs of factor Hamiltonian paths have at least one lift by the same two-list merge on both sides of v (the [Busch merge](../chanpak_near_floor_crossing_20261115/GLOBAL.md)). Paths with an invalid projection can add surplus but cannot cancel it. In particular,
```
L >= sum_{0<=i<=p,0<=j<=q} (F[i,j]G[q-j,p-i] - 1).
```
This is a necessary test on the cycle orders, not a computation of the entire Hamiltonian count.

## 2. A defect budget forces a nearby interleaving order

Call (i,j) a prefix bad cell when F[i,j]>1; put these in P. Put the converse bad cells in S in original *edge coordinates*, so (i,j) belongs to S if `G[q+1-j,p+1-i]>1`. All bad cells have both indices positive in their own matrix. In the sum in §1, a bad prefix and a bad suffix contribute at least one each even if charged to the same rotation fiber, since xy−1≥(x−1)+(y−1) for x,y≥1. Thus `|P|+|S|<=L`.

At any good interior prefix cell F=1, the orientation D is exactly its Y state. Indeed D=0 forces X=F[i−1,j]≥1, while D=1 forces Y=F[i,j−1]≥1. The corresponding statement holds in the converse.

The monotonicity required of a Ferrers matrix is violated by a 1 at (i,j) and a 0 at (u,w) with u≤i and w≥j. Consider pairs whose two cells lie outside P union S. At their southeast corner (i,w), the 1 propagates a positive Y along its row and the 0 propagates a positive X down its column. Hence (i,w) is in P. Applying the same reasoning in the converse puts their northwest corner (u,j) in S. These two bad corners uniquely identify the violating pair. We have at most `|P||S|` such pairs.

Erase the cells in P union S and one endpoint of every surviving violating pair. The number erased is at most
```
|P| + |S| + |P||S| <= L + floor(L²/4).
```
The retained colors extend to a Ferrers matrix: force a 1 wherever required by some retained 1 to the southwest (including the cell itself), and set other cells to 0. No retained 0 is forced to 1. This proves the edit bound. In particular H(T)=H(A)H(B) requires D itself to be Ferrers; conversely a Ferrers matrix gives exactly one lift of every pair of **cycle rotations**, not necessarily a lift of multiplicity one for other factor paths or the absence of broken paths.

A Ferrers matrix is the cross orientation of one interleaving of two formal chains of lengths p and q. Consequently there are at most `binom(p+q,p)*sum_{j=0}^r binom(pq,j)` possible cross matrices with surplus at most L, where r=`L+floor(L²/4)`. This is exponential in p+q times a polynomial for fixed L rather than exponential in pq. It is a necessary description even for mixed-direction seams; it does not assert the existence of a physical seam realizing a specified L.

## 3. A single reversal has only one possible low-cost location

Suppose D is all ones except a 0 at (i,j). If j>1, every prefix-length pair `(k,l)` with k≥i and l≥j admits two paths: one puts A before B; the second puts the first j B-vertices immediately before a_i, using the flipped arc b_j→a_i. Follow the rest of the A-prefix, then any remaining B-prefix. If j=1 use the same proof with k>i, so the first shuffle avoids the flipped arc. In the converse the analogous rectangle is disjoint in rotation-fiber coordinates: it has A-prefix length at most i−1 and B-prefix length at most j−1 in the original tournament. Thus
```
L >= (p-i+[j>1])(q-j+1) + (j-1+[i<p])i.
```
If i<p and j>1 this lower bound is at least p+q, by xy≥x+y−1 on each product. If precisely one of i=p or j=1 holds it is at least p+q−1, by minimization of the resulting affine function. The sole exception is (p,1), the pair immediately preceding v in A's cycle and immediately succeeding v in B's.

For the ordinary one-way join of exact floor chains with v at A's overall bottom and B's overall top, the exceptional edge goes from the last interior vertex of A's final link to the first interior vertex of B's initial link. This follows by merging cycles along the chain: at a common boundary the incoming arc from the preceding link's last interior vertex is replaced by the cycle through the later link and back to that boundary. Thus the exception belongs entirely to the two adjoining links. For the two-link core, with row and column indexed by link orders 3,4,5, the extra path counts are
```
2   4   6
4   6   8
6   8  10
```
These are exact via the subset recurrence in the reproducible script: initialize d({u},u)=1 and sum d(U−{w},u) over incoming neighbors u in U−{w} to obtain d(U,w). The arc flipped in the natural chain order, zero-indexed for link orders a,b, is (a−2,a). If there are any exterior links, use Busch multiplication with these unchanged strong links. The factors' original exact chain product is H(A)H(B), so the total increase is **at least** the table entry times the product of the exterior link counts, thus at least 6. Equality in that multiplication is not needed and can fail. Therefore no single cross-arc flip at a *canonical* floor split achieves excess 2 or 4 in the stability band. This last assertion does not classify arbitrary rooted special-chain factors with the lone reversal in the exceptional cell.

[Replay](verify.py): 8,848 matrices compared to literal shuffle enumeration, including the two-sided edit certificate and the rectangle bounds; all nine local values and 36 whole-tournament subset-DP controls. Next: use the actual off-cycle rotations and nonrotation factor paths to rule out (or construct) a few exceptional edits against a general Ferrers interleaving. Looking at just cycle rotations cannot do so: every Ferrers interleaving has all of those fibers simple.
