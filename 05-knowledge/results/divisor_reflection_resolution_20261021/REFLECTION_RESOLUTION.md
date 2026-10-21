# Two divisor reflections, graph complement, and the minimal path-chart resolution

2026-10-20/21. **PROVED elementary**, with finite exact controls. No literature-priority claim. The inherited divisor chart and ordered-graph correspondence are not claimed as new. The new contribution here is the joint-action obstruction, its sharp twelve-state repair, and the exact capacity criterion. No divisibility-preserving graph equivalence or unrestricted tournament-class bijection follows.

## 1. Retain the two actual reflections

For distinct primes p,q,r and integer k>=2, let D_k be the proper nontrivial divisors of p^k q r. Use exponent triples (a,b,c), with 0<=a<=k and b,c in{0,1}, omitting (0,0,0),(k,1,1). The inherited chart in trijection_reset_20261020/arithmetic.md is

 f_k=2a-k(b+c),  g=b-c.

It is injective on D_k; the only collision in the full box is its two removed endpoints. Divisor complement C and prime exchange Q act by

 C(a,b,c)=(k-a,1-b,1-c),  Q(a,b,c)=(a,c,b).

They commute. The two coordinate reflections are Q:(f,g)->(f,-g) and CQ:(f,g)->(-f,g); their product is C. The inherited fixed-point counts are

 |D_k|=4k+2, Fix(C)=0, Fix(Q)=2k,
 Fix(CQ)=2 if k is even, and0 if k is odd.

At k=2, sign(g) has sizes3,4,3 and sign(f) has sizes4,2,4. Their middle points are respectively Fix(Q) and Fix(CQ). These are exact reflections of one decoded set, not interchangeable numeric partitions. C is divisor complement, NOT the different prime-to-extra involution from session17_prime_shapes_20261020.

## 2. The ordered-graph dictionary transports the requested operations

Let G be a labelled simple graph on vertices0,...,m-1, m>=2. Make a tournament T_G on0,...,m by fixing the Hamiltonian path0->1->...->m. For every nonconsecutive pair i<j, orient j->i exactly when {i,j-1} is an edge of G. This is the inherited staircase bijection between all graphs on m vertices and tournaments with this fixed ordered path.

Write C_G for ordinary graph complementation in K_m and rho for vertex reversal i->m-1-i. Write K for reversal of all tournament arcs NOT on the fixed path, and R(T)=r(T^op), where r(i)=m-i. The full converse in R is essential: reversing the vertex order alone does not retain the forward fixed path. Directly on each pair,

 T_(C_G G)=K(T_G),   T_(rho G)=R(T_G).

Indeed (i,j) corresponds to graph pair(i,j-1), and tournament reversal sends it to(m-j,m-i), corresponding to graph pair(m-j,m-1-i), exactly rho(i,j-1). Complement toggles precisely the free bits.

Thus graph complement and actual fixed-path anti-reversal give commuting involutions on the same Boolean chart. R-fixed charts are the symmetric (blue) tilings for this fixed-path reflection. This does not identify all other paths of their underlying tournament classes.

## 3. A fixed-pair obstruction rules out the ten-state map

**Lemma.** If sigma is any involution of the vertices of a simple graph on at least two vertices, no graph satisfies sigma(G)=C_G(G).

Proof. If sigma has a transposition, its two vertices form a fixed unordered pair. If it has none, sigma is the identity and every pair is fixed. On a fixed pair, the equation would require its edge bit to equal its opposite. Contradiction.

In particular C_G rho has no fixed graphs. In the tournament chart, KR has no fixed states either.

**Theorem 1.** For even k, there is no map D_k into any such fixed-path chart that intertwines both C with K and Q with R. This is stronger than failure of a bijection: even a many-to-one equivariant map is impossible.

Proof. Choose x fixed by CQ. Its image would have to be fixed by KR, contrary to the lemma. If the two source axis reflections are interchanged, each source reflection still has fixed points, so the same obstruction applies.

Consequently the k=2 ten-divisor picture cannot be carried literally to graph complementation plus an involutive vertex/path reflection. The obstruction is fixed incidence of an unordered vertex pair, not insufficient graph cardinality. Self-complementary graphs do not contradict it: a complementing permutation need not be an involution.

## 4. A minimal twelve-state repair, with every source retained

Resolve only the two CQ-fixed divisors pq and pr. Replace each by two marked states (d,epsilon), epsilon in{0,1}; leave the other eight divisors single. Call the twelve-state set D'_2 and let pi:D'_2->D_2 forget epsilon. On the split states put

 C(d,epsilon)=(Cd,epsilon),
 Q(d,epsilon)=(Qd,1-epsilon).

Elsewhere retain the inherited actions. These are commuting involutions; pi respects both; CQ now exchanges the two states in each formerly fixed fibre. No original divisor is deleted or replaced by a different source.

Encode graphs on four vertices with edge order(01,02,03,12,13,23), first edge the low bit. Ordinary complement is x->63 xor x. Reversal rho=(03)(12) swaps bit positions0,5 and1,4, fixing2,3. The following explicit map is a bijection onto a twelve-graph subset:

| source | graph code |
|---|---:|
| p |0|
| pqr |63|
| p^2 |4|
| qr |59|
| q |1|
| r |32|
| p^2 r |62|
| p^2 q |31|
| (pq,0) |2|
| (pq,1) |47|
| (pr,0) |61|
| (pr,1) |16|

Substitution gives Phi C=C_G Phi and Phi Q=rho Phi for all twelve entries. Through section2 this is a concrete twelve-tournament marked-path carrier on five vertices, with the same operation laws. Its four R-fixed states give two blue K-pairs; the other eight states give four black K-pairs. These are selected labelled marked-path states, not a complete isomorphism-class metagraph.

**Minimality.** In any equivariant surjection Y->D_2 with CQ fixed-point-free on Y, each of the two CQ-fixed source fibres has even positive size, hence at least two. Every other fibre has at least one point. Therefore |Y|>=12, attained above. Four or fewer tournament vertices have at most2^3=8 fixed-path states, so five tournament vertices are also minimal for such a faithful realization.

The two additional states are an orientation sidecar, not two newly discovered divisors. The resolved source retains its actual divisor through pi. It is not an order embedding: p maps to the empty graph, contained in the image of q although p does not divide q. Multiplication, magnitude, and graph-edge inclusion are not transported. The precisely preserved structure is the joint involution action, source projection, and blue/black fixed-path predicate.

## 5. All-exponent orbit capacities and the parity boundary

This gives an exact action-level extension, not a canonical arithmetic embedding. For even k, split the two CQ-fixed points just as above; for odd k no split is needed. Write D'_k for the result. Its C2xC2-orbits are:

- k two-point orbits with stabilizer generated by Q;
- floor(k/2)+1 free four-point orbits;
- no orbit fixed by C or CQ.

For odd k, D'_k=D_k has4k+2 points. For even k it has4k+4 points and the two extra states are necessary by the same fibre argument.

On graphs with m>=2 vertices let N=binom(m,2), f=floor(m/2), t=(N+f)/2. Vertex reversal fixes exactly f edge coordinates; its edge-orbit count is t. Thus rho fixes2^t graphs, complement and C_G rho fix none. The graph action has

 2^(t-1) two-point orbits with stabilizer <rho>,
 (2^N-2^t)/4 free four-point orbits.

**Theorem 2 (capacity iff).** D'_k admits an injective equivariant map to this graph chart iff

 k <= 2^(t-1),
 floor(k/2)+1 <= (2^N-2^t)/4.

Necessity counts the two orbit types. For sufficiency choose the required distinct target orbits and map one representative of each source orbit to a target representative; matching stabilizers make its equivariant extension well-defined and injective. This choice is generally noncanonical. Before resolution an even-k source admits no equivariant map at all, by Theorem1. For odd k the iff applies to the unmodified divisor set.

## 6. Connection contract and reproduction

Source: full exponent-labelled divisors and the actual two chart reflections. Target: labelled simple graphs, then fixed-path tournaments. Original map: impossible for even k. Minimal repaired map: source-preserving two-fibre resolution followed by the explicit table. Preserved: both involutions, their product, original source through pi, and fixed-point colours. Lost: divisibility order, multiplication, numeric magnitude, unmarked class identity, and all-path incidence outside the selected carrier. Necessary sidecars: prime labels, vertex order, and exactly one bit on each of two central divisors.

Run python3 verify_reflections.py and python3 -O verify_reflections.py; compare the JSON receipts. The checker exhausts the displayed graph/tournament table, both operations, actual Hamiltonian-path orientation, blue status and pair counts; verifies divisor orbit types for k=2..12; and checks graph orbit counts on m=2..6 vertices. The uniform proofs above do not depend on a finite census. No relationship to Collatz chronology or idoneal completeness is asserted here.
