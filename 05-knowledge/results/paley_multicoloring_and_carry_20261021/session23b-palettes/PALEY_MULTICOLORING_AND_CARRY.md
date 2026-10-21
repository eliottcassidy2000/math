# Paley acyclic palettes, substitution, and an aperiodic bounded carry

Status: PROVED in the stated scopes, with an independent audit in the sibling session23b-palettes-review directory. This concerns acyclic vertex coloring, not proper unit-distance coloring and not distinguishing colors. No literature-priority claim is made. Ordinary dichromatic coloring of Paley7 and general lexicographic-product coloring are existing topics; the explicit constructions below are proved directly.

Let P be the Paley tournament on Z/7, with i→j when j−i belongs to {1,2,4}. A color class is legal when its induced tournament is transitive. A k-fold acyclic coloring assigns each vertex exactly k distinct colors, with every color's support transitive. Repeated supports are allowed because colors remain distinct. Write χ_k(P) for the least number of colors, with χ_0=0.

## 1. Exact integer multicoloring law

**Theorem.** For every integer k>=0,

χ_k(P)=ceil(7k/3).

P has no transitive four-vertex set: any first vertex in such an order would dominate the other three, but each out-neighborhood {i+1,i+2,i+4} is a directed triangle. Thus each color covers at most3 vertices and ≥ceil(7k/3) colors are necessary.

The following base palettes attain the bound:
- k1, three colors: {0,1,2}, {3,4}, {5,6}.
- k2, five distinct colors with supports {0,1,2}, {0,1,2}, {3,4,5}, {4,5,6}, {3,6}. Every vertex occurs exactly twice. The duplicate triple carries two different colors; deleting the duplicate would be an error.
- k3, seven colors: the translates of {0,1,2}. Every vertex occurs exactly three times.

All displayed triples are transitive, since their cyclic differences1 and2 are forward. Combine disjoint copies of the k3 palette with the appropriate k0/1/2 remainder. For k=3m+r the number is7m+0,3,5 respectively, exactly the ceiling lower bound.

The exact budget generating function is

sum_(k>=0) χ_k(P) z^k = (3z+2z²+2z³)/[(1-z)(1-z³)].

Its reduced denominator has degree4: the numerator is nonzero at1 and at both primitive cube roots. Hence the minimal eventual recurrence has order4, with χ_k=χ_(k-1)+χ_(k-3)-χ_(k-4).

## 2. A genuine substitution law, not multiplication of chromatic numbers

For tournaments T,U let T[U] be lexicographic substitution. An induced vertex subset of T[U] is transitive exactly when (i) its occupied base vertices form a transitive subset of T, and (ii) its part in every U-block is transitive. Necessity follows by taking representatives; sufficiency concatenates the two transitive orders.

Therefore, for any nonempty tournament U with dichromatic number k,

χ_d(P[U]) = χ_k(P) = ceil(7k/3).

For the lower bound, a coloring restricts to at least k colors per U-block. Choose k of those per block, retaining legal color supports on P; this produces a k-fold coloring. Conversely take a k-fold coloring of P and a fixed k-coloring of U, and biject each block's local k colors to its assigned global colors. The support criterion proves correctness.

In particular χ_d(P[P])=7, not χ_d(P)^2=9. Palette incidence, not just the scalar number of colors, is the preserved sidecar. Here the fractional dichromatic number is the minimum total nonnegative weight on transitive vertex subsets covering each vertex to weight at least1. In that sense χ_(d,f)(P)=7/3, with the uniform seven translated triples at weight1/3 attaining the counting lower bound. These are elementary cover certificates, not a numerical optimization assertion.

## 3. Iteration changes the recurrence class

Let P^[r] denote the r-fold iterated lexicographic power with P^[0]=K1. Its dichromatic numbers a_r obey

a_0=1, a_(r+1)=ceil(7 a_r/3).

They start1,3,7,17,40,94,... . Despite the rational generating function in the palette budget k above, the depth generating function sum a_r z^r is NOT rational over Q.

More generally, if A>B>1 are coprime integers and x_0 is a positive integer, x_(n+1)=ceil(A x_n/B) has a nonrational ordinary generating function.

Proof: set e_n=B x_(n+1)−A x_n, so e_n lies in {0,...,B−1}. If x were eventually C-finite over Q, then e would also be. A bounded integer C-finite sequence is eventually periodic: its finite windows of recurrence length range over a finite set, and the solved leading recurrence determines the next value from each window. Thus e has some eventual period L.

On each residue class after the onset, iteration gives x_(n+L)=r x_n+c with r=(A/B)^L and a fixed rational c for that residue. Subtract the rational fixed point α=c/(1−r): x_(n+tL)−α=r^t(x_n−α). Choose a fixed positive integer M clearing α's denominator. The left side times M is always an integer, while the right side is A^(Lt) N/B^(Lt) for the fixed integer N=M(x_n−α). Coprimality forces B^(Lt)|N for all t, hence N=0. This would make x bounded and eventually periodic, contradicting x_(n+1)≥(A/B)x_n and x_0>0. Therefore the depth sequence cannot be C-finite.

This proof also establishes that its bounded carry sequence e_n is not eventually periodic. A finite alphabet of carries is not an autonomous finite-state description of the integer dynamics.

For P powers, a_r/(7/3)^r increases to a finite positive limit: the increment is e_r/[3(7/3)^(r+1)] and 0<=e_r<=2. Thus1<=limit<=3/2. No exact value or transcendence claim for that limit is made, and no claim of non-D-finiteness is made for the depth series.

## 4. Fano and open-problem boundaries

If the seven Fano lines are all directed triangles, every two-coloring has a monochromatic Fano line, so the tournament has dichromatic number at least3. This uses the classical failure of Property B for the Fano triple system; an exact128-coloring control accompanies the proof. It does not say all triples are edges of a unit-distance graph, nor that Fano lines are Euclidean straight lines.

Distinguishing colors kill automorphisms. Acyclic colors forbid monochromatic directed cycles. Proper graph colors forbid monochromatic edges. They are three different predicates. Paley7 has distinguishing number2 but dichromatic number3; its underlying undirected graph K7 has proper chromatic number7 and cannot be a seven-vertex unit-distance clique in the Euclidean plane. Hence neither the two-color rigidity theorem nor the three-color acyclic theorem is a Hadwiger–Nelson lower bound.

The exact carry construction is relevant as a hostile to a proposed proof shortcut: bounded integer debts alone do not imply a periodic clock word or rational depth series. It does not model Collatz's native sources, guards, arrival clocks, or source budget, and proves no Collatz/LRC conclusion. A transfer would need an explicit map preserving those predicates.

## Literature context and provenance

A source search locates the established study of the minimal3-dichromatic tournaments (Neumann-Lara, 1994) and general acyclic coloring of digraph products (2024). The self-contained argument here makes no claim to originate those subjects or the abstract multicoloring substitution principle. Primary records: https://www.sciencedirect.com/science/article/pii/0012365X93E0113I and https://www.sciencedirect.com/science/article/abs/pii/S0166218X24000283 . The author preprint https://arxiv.org/abs/2204.14212 explicitly states the corresponding substitution principle in its abstract. These records were located, not treated as a newly audited proof of our exact formula. The explicit certificates and bounded-carry argument above stand independently.
