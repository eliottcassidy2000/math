# The second Busch gap and a growing forbidden band at fixed order

Write H(T) for the number of directed Hamiltonian paths, C(T) for the number of directed Hamiltonian cycles counted up to rotation, and f(n) for the minimum of H(T) over strong n-vertex tournaments.

**Results.**
- No strong tournament on n≥9 vertices has H=f(n)+4. The cutoff is sharp: the eight-vertex backward spine with gap word (2,0,2) has H=49=f(8)+4.
- At every order n≥3, if a strong tournament is not at the floor, its excess K=H−f(n) satisfies K≥n/5. For n≥9 it also satisfies
  `3K + 3 log_3 K >= n − 2`.
  In particular the first positive excess is at least n/3−O(log n).

These are exclusions at a **fixed order of a strong tournament**, not holes in the full Chan–Pak spectrum. The proof combines exact interleaving structure, a broken projection, and many possible replacements for the long block supporting that projection. It does not use extrapolation of a coverage census.

## 1. Inherited normalization

Use the [Busch stability theorem](../chanpak_busch_stability_20261115/README.md), the [sharp cycle inequality](../chanpak_near_floor_crossing_20261115/GLOBAL.md), and [exact seam rigidity and Douglas sorting](../chanpak_exact_seams_douglas_20261115/README.md). In particular:
- In the stable band K<2f(n)/25, all split-floor products are exact, all decomposition leaves are links of order 3, 4, or 5, and any nonfloor node with local surplus L costs at least wL, where w is the product of the floors of omitted leaves.
- An exact-floor strong factor is a [Moon–Yang special chain](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf) with the same link choices and at most one order-5 link. Each link count is 3, 5, or 9. Its unique cycle lists the junctions, including the two extreme marks, in reverse natural order, then all link interiors in natural order.
- In a strong split A,B meeting in v, every pair of valid factor path projections has at least one merged path. The broken paths, with an invalid projection, all add to surplus L=H(T)−H(A)H(B).
- If L<n−2, two factor cycle orders read after v give a Ferrers cross matrix. Equivalently there is one chronological interleaving of those orders with every cross arc directed forward. When both factors are exact-floor, exactly one Hamiltonian cycle has two valid projections.
- If C(T)=1, either H=f(n) or H≥47f(n)/45 (n≥5). Also H≥f(n)+3(C−1) in general.

## 2. Replacement lemma for a bad cycle in a Ferrers seam

Suppose both factors are exact-floor and L<n−2. If C(T)>1, take a cycle with an invalid projection. Assume first L<(n−1)/2. A projected invalid arc of A has a supporting segment of this cycle with t vertices strictly between its endpoints, all in B−v. All cuts outside the segment give broken Hamiltonian paths, so
```
t >= n−1−L.
```
Two invalid gaps in the same projection cannot coexist: their supporting edge segments are disjoint, and any cut would leave at least one invalid arc, giving n broken paths. Nor can both projections be invalid: each would require at least n−1−L intervening vertices from the other factor, a total greater than n−1. Thus B has its valid unique cycle order, and A has exactly one invalid projected arc. Write a=|A|, b=|B|, D for the t-vertex B-block inside the invalid gap, and O for other vertices of B−v. Set r=|O|. Then
```
a >= 4,       r = b−1−t <= L+1−a.
```
Indeed a directed triangle has no cyclic ordering with exactly one invalid arc. In B's own cycle, the complement O∪{v} of D is a contiguous run of r+1 vertices.

**Nonroot gap.** If neither endpoint is v, Ferrers chronology makes the entering and exiting arcs uniform to and from all D, because its cycle order is inherited from B. Every Hamiltonian path of D can replace the block to give a distinct bad Hamiltonian cycle. Each yields at least one broken rotation, and two different cycles cannot give the same broken path since the path uniquely determines its closing arc. Hence L≥H(D).

**Root gap.** Up to converse the invalid arc leaves v and enters an A vertex u. Then, reading the bad cycle from v, D comes first, before u and the remaining A vertices (including the held-back vertices of O in their original positions). Ferrers chronology guarantees all vertices of D beat u. Put R=B−O. For every Hamiltonian path of R with v not last, split it into P,v,Q with Q nonempty. Use P,v,Q,u, then the original suffix of the bad cycle from u to just before v. All arcs hold, and the resulting linear path's A projection starts with the invalid arc v→u. Distinct R paths give distinct broken paths in this one A projection. Thus
```
L >= H(R) − E_R(v),
```
where E counts paths ending at v. This count includes v with an arbitrary nonempty prefix, not just paths that start there. If R is strong, its Hamiltonian cycle gives at least |R|−1=t paths with v not last, which is already >L. If R is not strong, the original v,D path starts at v and covers R, so v is in the first strong component. All paths of R then qualify, giving L≥H(R).

## 3. The replacement bound grows with n

In the unique cycle of a special chain, a contiguous run of s≥1 vertices touches at most s+1 links. The first vertex touches at most two. Each next interior vertex can introduce at most one new link; each next junction shares a link with its predecessor (another junction, or the last interior vertex of the last link). Thus deleting a contiguous run of r+1 vertices leaves at least k_B−r−2 intact links if B has k_B links.

The intact links form one or more vertex-disjoint special-chain pieces after grouping successive links. Each piece contributes its chain product. The product over disjoint pieces is a lower bound in the whole remainder: combine their disjoint odd-cycle packings independently in the [Grinberg–Stanley formula](https://arxiv.org/abs/2307.05569), leaving all unused vertices as singletons. This argument does not require uniform cross-arcs between pieces; a deleted interior vertex can leave a surviving reversal between adjacent pieces. The lower bound is at least 3^k for k intact links. The same bound is valid for R since it contains D. Also b≤3k_B+2, because every link costs at most three vertices beyond its first, except for at most one order-5 link. Therefore
```
k >= (b−2)/3 − r − 2
  >= (n+2a−10−3L)/3
  >= (n−2−3L)/3.
```
In the nonroot case, or the nonstrong-R root case, §2 gives L≥3^k. This establishes 3L+3 log_3 L≥n−2. Using also 3^k≥2k+1 for integer k≥0 gives 9L≥2n−1 and thus L≥n/5 for n≥5. If L≥(n−1)/2 both conclusions already hold. If C(T)=1, the 47/45 ratio and parity imply L≥n/5 for all n≥5, and L≥n/3 for n≥9, since f(n)≥(15/2)n there. Thus the conclusions hold for every nonfloor union of two exact-floor factors whose product is f(n).

To pass to arbitrary strong T, outside the stable band K≥2f(n)/25≥n/3 for n≥9; at orders 3–8 the uniform bound K≥n/5 follows from evenness of each positive K. Inside the band, choose a deepest nonfloor node U in a Busch decomposition. Its children are exact floors; write m for its order and w for its omitted-leaf multiplier. The previous paragraph yields local K_U≥m/5. Each omitted leaf contributes a path count at least 1 plus its order cost, so w≥1+(n−m). Thus K≥w K_U≥m(n−m+1)/5≥n/5. Finally if n≥9 and K<n/3, no such nonfloor node can be proper: for 2≤m≤n−1, m(n−m+1)≥2(n−1) and 2(n−1)/5>n/3. Consequently the top factors are exact, and the two-factor result gives the logarithmic inequality for K itself. If K≥n/3 the inequality is immediate.

## 4. Exclude exactly four at every order n≥9

For K=4 and n≥9 we have 4<2f(n)/25, and the root factors are exact. The unique-cycle case is excluded by the [previous exact-seam and Douglas theorem](../chanpak_exact_seams_douglas_20261115/README.md), whose stronger integral consequence is H≥f(n)+6 if C=1 and H≠f(n) at n≥9. The sharp cycle inequality forces C=2. By Ferrers rigidity there is exactly one bad cycle. Its invalid projection is on A of order 4 or 5. In the boundary n=9 where 4=(n−1)/2, the two-invalid-projections argument still holds: both would need at least four vertices from the other side, and each would have size at least 5. This would require two factors of order 5; their floor product 81 exceeds f(9)=75. Thus B has one valid projection as before.

For K=4, r≤1 for a=4 and r=0 for a=5. In a nonroot gap the replacement bound is stronger than just L≥H(D): there is at most one bad cycle, so H(D)=1, which makes D transitive (a directed triangle would give at least three paths by insertion). Every link of B must then meet the deleted run O∪{v}. Its at most two vertices touch at most three links; for r=0 at most two. Hence n≤14 in this case (if a=5, n≤12).

In a root gap with r=0, R=B is strong and |R|−1≥5 because b≥6. This violates the root replacement bound. If r=1, then a=4 and O={z}, the predecessor of v on B's cycle (up to converse). If R is strong and b≥7, it gives at least five qualifying rotations. If R is not strong and B has at least four links, deleting z destroys at most two links; therefore H(R)≥9. Since b≥12 forces at least four links, n≥15 is impossible. This is the missing root-incidence case; the v-start gate alone would not exclude it.

The remaining orders 9–14 admit a finite structural replay, checking only Ferrers splits between exact floors. Two separate programs use the necessary condition that an invalid projection lies in a factor of order 4 or 5. The [subset-DP route](verify.js) filters by literal candidate cycles, using every possible position of the zero or one outside B vertex; it makes **no transitivity assumption on the block**, including at v. The [independent C++ route](independent.cpp) uses no candidate-cycle filter at all: it enumerates every rooted interleaving for these two factor sizes and counts actual directed paths by DFS, stopping at f(n)+5. Of 250,316 rooted matrices, 54,476 pass the candidate-cycle filter. Neither route finds f(n)+4. This finishes the exclusion. The eight-vertex boundary follows from the gap word (2,0,2) in the [exact spine transfer](../chanpak_exact_seams_douglas_20261115/README.md).

**Reproduce:** `node verify.js`; `g++ -O3 -std=c++17 independent.cpp -o /tmp/busch-gap-replay && /tmp/busch-gap-replay`. Construction counts include labeled rooted repetitions, not isomorphism classes.
