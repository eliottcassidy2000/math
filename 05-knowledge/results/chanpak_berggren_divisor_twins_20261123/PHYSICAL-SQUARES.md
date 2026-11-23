# Two literal squares: the forbidden prime centers and the 3–4–5 mixed differences

The two prime centers $2\cdot3+1=7$ and $2\cdot11-1=21$ themselves remain absent. The result here is an exact joint placement of both *bridges* on one physical tournament square, followed by an independent strong scaffold realizing the 3–4–5 relation.

## The two prime bridges share a marked arc

On six vertices use lexicographic pair mask $4$: bit one for pair $i<j$ means $i\to j$. Mark edges $e=(3,4)$ and $f=(0,5)$. The literal square is
$$
\begin{array}{c|rr}
 &e=0&e=1\\\hline
f=0&5&9\\
f=1&19&23
\end{array}
$$
So the horizontal reversal jumps each forbidden midpoint by four at the same physical edge; the vertical reversal translates both endpoints by fourteen. The mixed height difference is zero. The two upper corners are strong.

**All-order classification.** Up to simultaneous relabeling preserving the two individually marked edges, there are exactly two reduced squares, with masks $3655,3997$ in the normalization $e=\{0,1\},f=\{2,3\}$. Converse exchanges the two. Every square with these four heights is common ordinal singleton padding of one of the two. Its marked counting series is $2z^6/(1-z)^2$.

To prove reduction, the [gap-four classification](../chanpak_gap_cubical_geometry_20261122/LITERAL-BRIDGES.md) puts the upper row $19\leftrightarrow23$ on the same strong component of order six, with singleton components outside. The edge $f$ from the 19 vertex cannot change an outside interval disjoint from this core, because the factor 19 would survive and the required count 5 is not divisible by 19. If $f$ crosses component boundaries spanning any of the core, its reversal merges the whole interval into a strong component of order at least seven, whose minimum is 25 by [Busch](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3). Thus $f$ is internal, and all four corners have the same outside padding. Exhaustive enumeration of labeled tournaments and marked edges at order six gives 1440 admissible ordered squares and the two canonical masks. The code also checks that in either reduced square the 9 corner has one outside singleton and its strong five-core has canonical mask 8, the transitive deletion word OOOI (O indicates an arc from the added vertex to the ordinary vertex). That core is not OIOI, whose canonical mask is 10.

## A strong physical 3–4–5 identity

Start with the strong order-six tournament of mask 81 (row masks $34,8,3,5,15,30$) and height 23. Independently reversing one of $(2,5),(2,4),(1,5)$ produces another strong tournament of height 29, 31 or 33, respectively. The three positive gaps are therefore $2(3,4,5)$. Two ordinal copies with two independent reversals of the same type form a square on twelve vertices. Its mixed difference is the square of the corresponding gap, because $H(U\to V)=H(U)H(V)$. This gives the literal relation $36+64=100$ on a common labeled scaffold, with the edge type selecting the term.

There is also an all-strong version on a fixed fifteen-vertex scaffold. Put the interior $B$ between two singleton anchors $a,b$ in a directed cycle of blocks $a\to B\to b\to a$. Add singleton $z$ with $B\to z\to b$. The switchable third edge is $az$. Write $Q=Q_2(B)$ for the count of coverings of all vertices of $B$ into two ordered possibly empty directed paths. In the $a\to z$ orientation the height is $H(B)+2Q$; in the $z\to a$ orientation it is $-H(B)+3Q$. Both shells are strong for every nonempty $B$: the cycle through $a,B,b$ connects every interior vertex, and $B\to z\to b$ attaches $z$ in both directions.

The two shell formulas are exact for arbitrary B. Removing B from any full path gives an ordered outside skeleton with up to four slots for nonempty paths of B. Hence the shell read is a fixed linear combination of the allowed-empty counts $Q_j(B)$, $1\le j\le4$. At transitive B of sizes one to four, $Q_j(B)=j^{|B|}$. The literal four checks give $(5,9,17,33)$ and $(5,11,23,47)$ in the two orientations; the $4\times4$ Vandermonde at bases $1,2,3,4$ uniquely determines the coefficients. The replay independently checks all eight physical corners of all three strong cubes.

For the two copied interior switches let $D_0,D_1$ be their mixed differences read in the two wrapper orientations. Then
$$3D_0-2D_1=5(\Delta H)^2,$$
since the combination cancels the auxiliary $Q$ before taking the mixed difference. The three strong cubes give raw pairs $(4644,6876),(6336,9344),(15588,23132)$, and therefore the exact identity $180+320=500$. This is a relation between signed differences of actual Hamilton heights, not an assertion that even numbers are attainable as Hamilton heights. The third reversal exposes the arithmetic register without assuming it is an absolute count.

Run `python3 pythagorean-physical-check.py` for the classifications and strong-cube controls.
