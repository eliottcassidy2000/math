# All orders: the pairwise-intersecting triangle class

A hub means a vertex on every directed triangle, equivalently a vertex whose deletion leaves a transitive tournament. Call a tournament *intercyclic* if it has no two vertex-disjoint directed triangles (equivalently no two vertex-disjoint directed cycles).

**Theorem.** For every $n\ge6$, the strong intercyclic tournaments fall into exactly $2^{n-3}+n-3$ isomorphism classes. Of these $2^{n-3}-1$ have a hub. The $n-2$ hub-free classes split as follows:

| two-hub word | parameters | ordinary Hamiltonian count | number of classes |
|---|---|---|---|
| $1^a2\,1^b$ | $a,b\ge1,\ a+b=n-3$ | $1+3\cdot2^{n-3}$ | $n-4$ |
| $1^{n-4}20$ and its converse $321^{n-4}$ | both | $9\cdot2^{n-4}-3$ | 2 |

In particular at $n=10$ the total 135 is 127 hub classes plus 8 hub-free ones. **Six** of the latter have count 385 and **two** have count 573. At $n=5$ the hub-free part consists of one class of count 13 and two of count 15 (one of the latter is the regular five-tournament). At orders at most four there is no hub-free strong intercyclic tournament of order at least three.

This proves the all-order prediction P16 of the [earlier two-hub waves](https://github.com/eliottcassidy2000/math/blob/f66f56b51fdfeff009159575736e4f7932815eea/05-knowledge/results/chanpak_waves_20261113/README.md) in its relevant range. It does not use the finite census as an induction hypothesis.

## 1. An exact regular language for all admissible ordered hub choices

[Bang-Jensen and Bessy, *Small cycle transversals in tournaments*, Theorem 1.8](https://www.lirmm.fr/~bessy/publis/SmallTransversTournament.pdf) supply the all-order transversal bound two. Thus if no single hub works, choose two covering hubs $x\to y$. With an ascending spine, put $\delta(v)=2[x\to v]+[y\to v]$. The following elementary predicates are necessary and sufficient.

1. **Strongness.** The first letter is not 0 (otherwise a source), the last is not 3 (otherwise a sink), and $y$ can reach $x$. The latter holds iff there is a letter 1 or a subsequence 30: a shortest route between the hubs through the increasing spine needs at most two spine vertices. When the hubs reach each other, the first and last spine conditions attach every site to their strong component.
2. **No disjoint triangles.** An $x$-only triangle is a bit-$x$ inversion at an ordered pair of sites; analogously for $y$. A disjoint pair must use distinct hubs, so consists of an $x$-inversion and a $y$-inversion on disjoint sites. The union of their four sites gives a forbidden length-four subsequence. There are exactly 72 such digit strings: generate them by the six choices of the two $x$-sites and test both inversions.
3. **Neither named hub suffices.** Each of those two types of inversion must occur.
4. **No site is a hub.** After deleting a site the remaining tournament is transitive iff the new word lies in $0^*2^*3^*$. In fact transitivity eliminates letter 1 (the triangle using both hubs) and both types of inversion. Therefore no deletion of one site can give that language.

All predicates are finite-state. The forbidden-subsequence test remembers which proper prefixes of the 72 strings have already occurred as subsequences. A one-deletion transitivity test needs only a no-deletion state and a set of possible states after spending the deletion. The conjunction recognizes exactly the union
$$
1^+20,\quad321^+,\quad1^+21^+,\quad1^+3^+0,\quad2^+10,\quad30^+1^+,\quad312^+,\quad310,\quad320.
$$
The certificate constructs the product of the predicate automaton and a nondeterministic pattern automaton represented by subsets. It visits **220 reachable paired states**, proves acceptance agrees in each, and checks all four successor letters. Since every word reaches one of these states, this is an all-length language equivalence, not evidence from a large cutoff.

## 2. Remove the dependence on a named pair of hubs

For the graph $A_{a,b}$ with word $1^a21^b$, let $E,L$ be its early and late blocks of 1s and $z$ its single 2. Rerooting at $(z,y)$ gives spine $(E,L,x)$ and word $1^a3^b0$. At $(x,z)$ the spine is $(y,E,L)$ and the word is $30^a1^b$. If $b=1$, let $v$ be the sole late vertex; rerooting at $(y,v)$ gives $2^a10$ on $(E,x,z)$. Taking converses gives the $312^+$ family. Thus all seven unbounded word patterns reduce to the two rows of the table.

For $a+b\ge3$, the three vertices lying on more than two directed triangles in $A_{a,b}$ are exactly $x,y,z$. Their triangle counts are $a+2b,\ 2a+b,\ a+b$, respectively; the last one is strictly smaller than the other two. Each early leaf is in the two triangles on core pairs $xy,yz$, and each late leaf in the two on $xy,xz$. The arc $x\to y$ orders the two remaining core vertices, so the pair $(a,b)$ is an isomorphism invariant. This explains why a single scalar can hide $n-4$ distinct overlaps. In the $A$ family the minimum and maximum outdegrees are $1,n-2$; in $1^{n-4}20$ they are $2,n-2$ for $n\ge6$, and in its converse they are $1,n-3$. The two $B$ graphs are distinct from one another and from every $A$ graph. The short exceptions at order five are directly replayed.

## 3. Exact path counts and the 127 hub classes

For fixed homogeneous blocks of positive length $r$, a Hamiltonian path can use a block in one, two or three run slots. Fix which slots are occupied, including the hub order and the rest of the assignments. Seam eligibility is independent of the internal distribution among those slots. The number of surjections onto them is $1,\ 2^r-2,$ or $3^r-3\cdot2^r+3$. Therefore counts lie in the span of $1,2^r,3^r$ in each length independently. The exact subset DP in the replay gives the $3\times3$ coefficient grid for $A$ and three coefficients for $B$. Invertibility of the sampling basis at lengths $1,2,3$ proves the formulas in the table at all lengths.

A strong one-hub tournament has a binary spine word starting 1 and ending 0, giving $2^{n-3}$ marked words. Every other candidate hub on the spine must be one of those extreme sites, because they form a triangle with the marked hub. Deleting the first site then forces the remaining word to be all zeros; deleting the last forces all ones. For $n\ge4$ these two distinct words are the same one-backward-arc tournament with the other hub as root. All other words have a unique hub and are recovered by its deletion. This proves the subtraction of exactly one.

## 4. Prescribed target and why this band cannot select a ternary overlap

In the hub-free intercyclic subclass a prescribed $N$ is now solved completely. Test if $(N-1)/3$ is a power of two at least $2^3$, or if $(N+3)/9$ is a power of two at least $2^2$, and handle the order-five values separately. The exponent gives the order. In the first case **any** positive split $a+b=\log_2((N-1)/3)$ constructs the desired tournament, and all its hub-free intercyclic witnesses occur in the described table. In the second case output either of the two converse words. Thus 385 gives the six-way split of 7, and 573 gives the two directions with exponent 6.

No frame with a triply eligible free site can occur here. Its four anchors would give disjoint triangles, one for each hub: with a path ordered $R_0,x,R_1,y,R_2$ put $a=\max R_0$, $b=\min R_1$, $c=\max R_1$, $d=\min R_2$. If site $i$ belongs to all three eligibility intervals then $b,d<i<a,c$. So $x\to b\to a\to x$ and $y\to d\to c\to y$ are disjoint. Hence $\beta=0$ in every frame of every such chart. The 385 and 573 examples are exact *binary* controls, not hidden three-choice devices: a universal exact-target strategy that needs genuine ternary overlap must leave this band.

This also distinguishes the small numbers: 3 path slots; 6 permutations of three nonempty pieces; 12 potential full orders after choosing a hub order; 36 states in the transfer counter (2 orders, the two 3-status intervals and one 2-status interval). 127 is the marked-word quotient $2^7-1$, not a ternary-branch count. No domination conjecture is assumed in the proof.
