# Linek's tree problem at activity two

[Law's rooted theorem](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v17i1n18/pdf/)
realizes every pair (root-omitting count, total independent-set count)
modulo every modulus by trees at activity one. Linek's exact question,
whether trees miss only finitely many positive integers, remains open
in [Chan–Heilman–Panova's analysis](https://arxiv.org/abs/2604.19094),
which connects the inverse problem to continued fractions. The activity
for Hamiltonian paths is instead **two**, by the
[Grinberg–Stanley odd-cycle expansion](https://arxiv.org/abs/2307.05569).
The change of activity matters.

## Scalar and rooted local behavior part ways at activity two

For a rooted tree let a be the weighted count omitting its root and b
the total, with weight 2 per chosen vertex. Attaching a new root gives
F(a,b)=(b,b+2a). Joining two roots to a new root gives
(bd,bd+2ac)=F(ac,bd) for pairs (a,b),(c,d).

**Every scalar residue is realized modulo every odd modulus.**
Fix an odd modulus M. The map F is invertible modulo M, and its inverse
is physically a nonnegative power of F. In particular the available
pairs are closed under coordinatewise multiplication by the two-branch
construction followed by F inverse. Starting from (1,3), successive
inverse steps yield (1,1),(0,1),(1/2,0),(-1/4,1/2).
An attainable second coordinate u gives diagonal scaling by u:
multiply its pair by (0,1) and apply F to get (u,u). Thus we get
diagonal scaling by 1/2 and, by its finite order, by 2. We obtain (1,0)
and then (-1/2,1) by F inverse. Multiply these last two coordinatewise
and scale by 2 to obtain (-1,0). Applying F gives (0,-2), so diagonal
scaling by -1 is allowed as well.

Write c=-1/2 and D=diag(c,1). Multiplying by (c,1) realizes D,
and diagonal scaling by c inverse gives diag(1,c inverse) after D.
The conjugation
$$
 \operatorname{diag}(1,c^{-1})FDF^{-1}(a,b)
    =(a,b+((1-c)/c)a)=(a,b-3a)
$$
is consequently physical modulo M. The second coordinates of
(1,0),(1,1),(-1,-1) represent the three classes modulo 3,
and in each case the first coordinate is a unit. Repeatedly applying
the shear gives every residue if 3 does not divide M, and every residue
in the chosen class modulo 3 otherwise.

**But the two-coordinate statement is already false modulo 3.**
The exact pair set there is
$$
 \{(0,0),(1,0),(0,1),(1,1),(-1,0),(0,-1),(-1,-1)\}.
$$
The six nonzero pairs form the F-orbit of (0,1); their coordinatewise
products, together with (0,0), stay in the displayed set. These two
operations generate rooted trees: for a root with any number of
children first multiply the children's pairs, then apply F. The
singleton is in the set, and (0,0) is the product of (1,0) and (0,1).
Thus opposite nonzero coordinates are impossible, even though the
scalar projection is onto.

## Exact scalar test at powers of two

For any nonempty tree on n vertices, the empty set, singletons,
and independent pairs give
$$
 b\equiv 1+2n+4\left(\binom n2-(n-1)\right)
   =2(n-1)^2+3\pmod 8;
$$
larger independent sets have weights divisible by 8. So the possible
classes are 3 for odd n and 5 for even n. Conversely, at every
modulus 2^e with e≥3, stars have b=3^d+2 for d leaves and attain
every class 3 or 5 mod 8. Indeed the order of 3 is 2^(e-2),
since successive squaring gives v₂(3^(2^k)-1)=k+2 for k≥1.
Its powers lie in the two classes 1 and 3 mod 8, which together
have the same size. The assertion for e<3 is the projection of
the e=3 result. No simultaneous claim for mixtures of odd and
even moduli is needed.

In particular 7 is forbidden for connected activity-two trees
by a simple congruence, but 21 is attained by the four-vertex path
(1+8+12). A weighted variant of Linek's question is eventual
coverage of every integer in the admissible 3 or 5 mod 8 classes.
The separate odd and even tests do not prove simultaneous saturation
or exact coverage.

## A quadratic clique bound prevents the literal tournament transfer

Let Ω(T) be the graph whose vertices are directed odd cycles,
adjacent when they share a tournament vertex. Its components
are precisely the nonsingleton SCCs of T, by the
[universal-vertex argument](CONFLICT_GRAPHS.md).

In a strong SCC of size m, any two odd cycles of length >m/2
intersect. If r is the number of permissible odd lengths j≥3
with m/2<j≤m, [Moon's minimum cycle count, Corollary 2.1](
https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf)
gives a clique of size at least r² for m odd or r(r+1) for m
even, by summing m-j+1 over these lengths.

A sharper asymptotic bound comes from a different clique: all odd
cycles through one fixed vertex. Average vertex–cycle incidence gives
$$
 \omega(\Omega(T))\ \ge\ {1\over m}\sum_{3\le j\le m,\ j\ {\rm odd}}
 j(m-j+1) = m^2/12+O(m)
$$
for a strong SCC of size m. The quadratic scale is tight: orient the
consecutive pairs on 1,...,m as i→i+1 and all other pairs backward.
This tournament is strong, with the directed spanning cycle that
closes from m to 1. In any cycle, start at its least vertex; the
first downward step must close the cycle, since otherwise it
revisits an interior vertex. Thus the odd cycles are exactly the
consecutive intervals of odd length at least three. Pairwise
intersecting intervals have a common vertex. The maximum clique
is therefore the maximum number of odd intervals through one
vertex. At position t this count is t(m-t+1)/2+O(1):
choose a start at or before t and an end at or after t of the
same parity, then remove the singleton. Maximizing at the center
gives m²/8+O(m).

It follows that uniformly bounded conflict clique number bounds
strong SCC orders, though not total tournament order. The
Hamiltonian path count is the product of the SCC path counts, so
the possible values use a fixed finite set of s primes and number
at most O((log X)^s) up to X. Also, if p divides H(T), there is
some SCC of size m with p≤H(SCC)≤m!. Thus the full conflict
graph for any such witness needs clique number at least
(1/12-o(1))(log p/log log p)^2 as p grows.

## Sharp classification when the clique number is at most three

In fact Ω(T) with clique number at most three is a disjoint union
of isolated vertices and isolated edges. For m≥7 Moon's long-cycle
clique above already has at least four vertices. At m=5 there are
at least three triangles, all intersecting one another, and a
Hamiltonian cycle. At m=6 there are at least two 5-cycles and
four triangles; some two triangles intersect, and both intersect
both 5-cycles, again giving a four-clique. For m=3 a strong component
has one triangle. At m=4 its scores are 1,1,2,2, so there are
exactly two triangles (subtract the two transitive triples from
the four triples); they intersect.

Consequently these restrictions are equivalent: the full conflict
graph is a forest; it is triangle-free; it has clique number at
most three; the nontrivial SCCs all have size 3 or 4. The exact
Hamiltonian spectrum in this class is
$$
 \{3^a5^b:a,b\ge 0\},
$$
including 1 from a transitive tournament. Indeed, the independent-set
polynomials of an isolated vertex and edge evaluate to 3 and 5
at activity two, and ordinal sums realize all such products.
In particular no full conflict graph has clique number exactly 3.
The four-vertex path's value 21 was an independent-set
calculation on a tree; no tournament can use that tree as its full
conflict graph. The restriction here is on the full graph, not on
chosen induced conflicts.

Thus the analogue with Linek is genuine at the level of generating
functions and inverse enumeration, but the unrestricted transfer
would need to control extra odd cycles, not merely encode the input
tree. The bridge and initial clique count came from the background
research pass; the weighted modular analysis, incidence bound and
sharp clique-at-most-three classification extend it. See
[the finite control script](replay_linek.py).
