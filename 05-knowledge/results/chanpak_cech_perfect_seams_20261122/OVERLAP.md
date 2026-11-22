# Three overlapping faces and the path interaction they miss

**PROVED**, with explicit exact examples. In the fixed path $1\to\cdots\to n$, write a free tile as $(a,b)$ where $a-b\ge2$. The three faces from [THM-553](../../../01-canon/theorems/THM-553-tile-address-recurrence.md) and [THM-801](../../../01-canon/theorems/THM-801-mobius-cech-descent-and-pure-cubic-colour-law.md) are $A:a<n$, $B:a-b\ge3$, $C:b\ge2$. The upward embeddings for the faces and overlaps are
$$
 (a',b')\longmapsto (a'+y+z,b'+z),\quad
  (x,y,z)\in\{0,1\}^3\setminus\{(0,0,0)\},
$$
where $x,y,z$ select $A,B,C$. At degree one the lower order is $n-1$, at degree two it is $n-2$, and $G$ is their order-$n-3$ common overlap. The tiled orientations glue uniquely if the three pair overlaps agree. The sole order-three tile is the seed. This is the original seven-term recursion on *addresses*. The B operation is the gap-window tournament of [THM-830](../../../01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md), not ordinary vertex deletion.

## A local cubic term of the actual path count

Fix all tiles except the three unique-to-one-face corner tiles $\{1,3\},\{1,n\},\{n-2,n\}$. For $n\ge$\{1,6\}$ these form one undirected path. Code bit 1 for each edge in the direction $3\to1\to n\to n-2$. The full path count as a multilinear polynomial in these three bits has coefficient of their product
$$
      \kappa = \Gamma_{3,n-2}-\Gamma_{n-2,3}.
$$
Here $W=[n]\setminus\{1,3,n-2,n\}$, and
$$
 \Gamma_{p,q}=\sum_{S\subseteq W} L_p(S) R_q(W\setminus S),
$$
where $L_p(S)$ counts paths ending at a vertex beating $p$, and $R_q(S)$ counts paths starting at a vertex beaten by $q$, both 1 on the empty set. Each counted path contributes the product of its required corner literals. The third mixed derivative retains precisely paths containing all three corner edges, which have to be traversed consecutively in one direction or the other. The two orientations have cubic signs $+1$ and $-1$. The two slots count the exterior portions with no assumption of adjacency between the slots when the block is removed. This is a direct bounded-boundary contraction, not a count of paths in a tournament of order $n-3$. If the two outer vertices have identical contacts on $W$, the cubic vanishes.

The coefficient can grow exponentially even with transitive $W$. For $n\ge9$, let $W$ have its increasing orientation, let every member except 4 beat $p=3$, and let $q=n-2$ beat every member except $n-3$; these contacts respect the forced spine. Write $m=n-4$. Each subset has a unique path. The forward count forbids left slot ending at 4 (two subsets) or right slot starting at $n-3$ (two, disjoint), so it is $2^m-4$. In the reverse count both slots must be nonempty, the left must end at $n-3$ and the right start at 4. Four specified vertices then have forced slots (2 and $n-3$ left; 4 and $n-1$ right), leaving $2^{m-4}$ choices. Thus $\kappa=15\cdot2^{m-4}-4$. On the cube of three corners this prohibits any additive formula built from functions each ignoring at least one of the corners, independently of how much common overlap data each function retains.

## Seven scalar reads fail even to tell a unit move from a non-unit move

Sort tile bits by increasing $a$ then $b$, least significant bit first, with 1 for $a\to b$, and fix $1\to2\to\cdots\to n$. At $n=$\{1,6\}$, masks 287 and 406 have the same seven lower counts. Reversing corner ${1,6}$ yields identical lower count changes too. The actual upper moves are respectively $23\to25$ and $19\to29$. Using face order $C,B,BC,A,AC,AB,ABC$, both pairs have seven lower reads $(11,11,5,9,5,5,3)\to(11,13,5,9,5,5,3)$. At $n=5$ masks 18 and 19 give $13\to15$ and $11\to11$ with all seven reads identical and unchanged. These examples disprove a function of the seven *before and after* counts deciding the +2 property. They do not dispute exact gluing of marked tiles. See [check_faces.py](check_faces.py).

## A multiplicative repair on path-independent cycle objects

In the commutative squarefree algebra $\mathbb{Z}[x_v:v\in V]/(x_v^2:v\in V)$, define $Z(D)=\prod_C(1+2x_{V(C)})$ over odd directed simple cycles of a directed graph $D$, counted up to rotation. The sum of all coefficients of $Z(T)$ equals $H(T)$ **for tournaments** by the [Grinberg–Stanley expansion](https://arxiv.org/abs/2307.05569). For incomplete digraphs $Z$ is a cycle ledger, not that graph's path count. In particular, retain all $n$ vertices and the forced spine in each face graph $D_A,D_B,D_C$, but retain the optional arcs only on that face. Their union is the full marked tournament. For each intersection use the actual common arcs. Then
$$
 Z(T)=Z(D_A)Z(D_B)Z(D_C)\,Z(D_{AB})^{-1}Z(D_{AC})^{-1}Z(D_{BC})^{-1}
       Z(D_{ABC})\prod_{C\text{ mixed}}(1+2 x_{V(C)}),
$$
where a mixed cycle is contained in none of the three face graphs. For any nonmixed cycle its exponent is $k-\binom{k}2+\binom{k}3=1$; each inverse of an individual cycle factor is $1-2x_{V(C)}$. This works equally with induced covers. Disjoint mixed cycles can be selected first, deleting their vertices from all the face calculations. Summing coefficients only at the end is essential: ordinary multiplication after summing would no longer enforce disjointness. Here the exact seven-term *cycle-interaction* repair is division by overlap ledgers and inclusion of mixed cycles, not an additive recurrence for the seven lower Hamilton counts.

For an induced two-bag cover $A,B$ of a tournament, with separator $W$, suppose a changed arc is in $A$ but not $B$, and neither orientation admits any mixed odd cycle. The sign correction on $W$ is $Z(W)^{-1}$. If $b_Y$ is the coefficient of $x_Y$ in $Z(B)$, and $m_Z$ that of $x_Z$ in $Z(W)^{-1}$, then the normalized seam read satisfies
$$
 \delta_{A\cup B}(\varnothing)
 =\sum_{Y,Z\subseteq W\cup(B\setminus A),\; Z\subseteq W,\;Y\cap Z=\varnothing}
 b_Y m_Z\,\delta_A((Y\cup Z)\cap A),
 \qquad \delta_A(S)=(H(A'-S)-H(A-S))/2.
$$
All coefficient indices are vertex sets; the formula follows by fixing the supports from the unchanged factors. Consequently the gcd of the normalized deletion defects on the interface divides the glued unit read, and mod 2 the glued read equals the original read (nonconstant factors are even). This stronger *labelled* bookkeeping is the relevant interface data for a recursive arc move.
