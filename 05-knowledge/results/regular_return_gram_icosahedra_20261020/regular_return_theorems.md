# Endpoint response and the regular two-complement return graph

**Status: PROVED.** The all-order statements below are elementary deductions with complete proofs. The order-five graph is classified by a structural argument and a finite explicit adjacency certificate, not by a large census. Independent executable checks are supplementary. None of these results proves the ordinary Hamiltonian-path spectrum conjecture or connectivity of the black merged metagraph.

## 1. Scope, inheritance, and objects

Write H(T) for the number of directed Hamiltonian vertex words, before quotienting by automorphisms. For a Hamiltonian path P, K_P(T) reverses every arc outside P and preserves the path arcs. Our graph G_n has the labelled regular tournaments on a fixed n-element set as vertices. Distinct T,U are adjacent when

    U = K_Q(K_P(T))

for directed Hamiltonian paths P in T and Q in K_P(T), with U regular. The intermediate tournament is not required to be regular. This is a new derived two-step graph, not the original metagraph and not a claim that either constituent move is black.

Already-established inputs, explicitly not claimed new here:

* `05-knowledge/results/spectrum_basins_20261020e/regular.md` gives Ryser's score-fibre connectivity, its elementary cycle-to-triangle compiler, the regular order-seven census, the one-step endpoint defect, and the fixed-label R-parity obstruction.
* `THM-781-hamiltonian-path-inverse-metagraph-fibre.md` proves the free path-orbit quotient H/|Aut(T)|. Raw H and this quotient are not interchangeable.
* `PROP-001-arc-flip-identity.md` contains the arc-flip/odd-cycle identity. Our punctured-triangle endpoint formula below is a separate direct counting decomposition.
* `THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md` distinguishes a response lattice from an actual interval.

The concept board is endpoint charge; score fibre; labelled orientation parity; successor permutations; endpoint-refined counts; raw H versus path-orbit mass. The hostile is the regular order-seven spectrum {171,175,189}: connected triangle motion has increments with gcd two but leaves holes. The useful additional coordinate is the ordered endpoint pair and the full response kernel, not a scalar degree or count.

## 2. Exact two-step regular return criterion

For any tournament T and Hamiltonian path P from s to t, put delta(P)=e_s-e_t. Reversing all arcs and then restoring the path arcs gives

    d(K_P(T)) = (n-1)1 - d(T) + delta(P).                 (1)

For a second path Q in K_P(T), therefore,

    d(K_Q(K_P(T))) = d(T) - delta(P) + delta(Q).          (2)

For n>1, equality of endpoint dipoles means equality of the ordered endpoint pair. Thus a two-step complement returns to the same labelled score vector if and only if P and Q have the same start and the same end. In particular this is necessary and sufficient for return to regularity when T is regular.

The net flip set is the symmetric difference of the two undirected path-edge sets:

    F = E(P) symmetric-difference E(Q).                  (3)

A Q-edge shared as an unordered pair with P must have the same direction, because K_P(T) preserves P. Every other Q-edge points opposite its original direction in T. Consequently F, oriented as in T, is the union of the P-only arcs and the reversals of the Q-only arcs. Equation (2) makes this directed difference balanced. It decomposes into directed cycles, but this does not make the two complement moves identical to arbitrary cycle reversals.

### Parity and loss at isomorphism quotient

Both path-edge sets have n-1 edges, so |F| is even. Every G_n edge preserves

    epsilon(T) = sum_{i<j} [i -> j in T] mod 2.           (4)

For odd n>=3 both signs occur among regular tournaments: any transposition of the labels changes epsilon. Hence G_n has at least two labelled components. A directed triangle reversal changes epsilon and cannot be one such edge, or a sequence of such edges.

Under relabelling by a permutation pi,

    epsilon(pi T) = epsilon(T) + parity(pi) mod 2.       (5)

Indeed the only sign changes come from pairs whose order is reversed by pi; their number is its inversion parity. Thus epsilon is not an invariant of an unlabelled isomorphism class. Every ordinary class contains both signs. This obstruction supplies no disconnection theorem after quotienting by isomorphism, or after converse merging.

## 3. Minimum support and a precise six-arc move

**Theorem.** Every nontrivial G_n edge reverses at least six tournament arcs. If exactly six arcs change, its defining paths exchange two consecutive nonempty interior blocks:

    P = A B C D,       Q = A C B D,                     (6)

where A and D are nonempty and retain the common endpoints. Both middle blocks cannot be singletons. Conversely, any such exchange for which P is directed in T and Q is directed in K_P(T) gives a score-preserving two-step move with six changed arcs.

**Proof.** Append the same artificial closing edge t->s to both paths. Their successor maps sigma and tau are n-cycles. Because no common underlying edge is traversed oppositely, if k path successors change then exactly 2k tournament edges change. One successor cannot change in a permutation. Two changed successors exchange their two values, splitting a single n-cycle into two cycles, so tau would not be an n-cycle. Therefore k>=3.

For k=3, cutting P at the three changed arcs makes four nonempty directed blocks A,B,C,D. Every internal edge of each block is unchanged, the first block still begins at s, and the last still ends at t. The only different possible block order is A,C,B,D. If B and C were singletons, their common edge would be traversed in reverse, contrary to Q being directed after K_P. Conversely the three boundary replacements of any legal exchange change three old edges and three new edges, proving the assertion. This also proves the support statement without relying on the sign of a successor permutation. (That sign independently implies sigma^{-1}tau is even.)

If B and C both have at least two vertices, the six changed arcs form one directed six-cycle. To see this, write x=end(A), b=first(B), b'=end(B), c=first(C), c'=end(C), y=first(D). In T the six-cycle is

    x -> b -> c' -> y -> b' -> c -> x.                   (7)

If exactly one middle block is a singleton, the repeated boundary vertex in (7) splits the support into two directed triangles sharing that vertex. These are the only support types for a six-arc return. In particular two vertex-disjoint directed triangles cannot be a six-arc return.

## 4. Exact endpoint-refined return mass

Let H_{s,t}(S) count directed Hamiltonian words in S starting at s and ending at t. Count path-pair presentations of two-step regular returns, with multiplicity, by

    M(T) = sum_{P in Ham(T)} H_{start(P),end(P)}(K_P(T)). (8)

This is exact by (2). Exactly H(T) presentations return to the identical labelled tournament: by (3), this requires E(P)=E(Q), hence Q=P. Thus

    nonidentity presentation mass = M(T)-H(T).           (9)

Distinct path pairs can yield the same neighbour, so (9) is not a graph degree. If w(T,U) is the number of presentations leading to U, then w(T,U)=w(U,T): reverse a presentation by exchanging P and Q around its common intermediate tournament. Aut(T) acts freely on its outgoing presentations, since it already acts freely on P. Hence M(T)/|Aut(T)| is an integer presentation-orbit mass, not raw H and not a count of distinct neighbouring classes.

These statements apply to any fixed labelled score fibre, not only regular tournaments, when the same-score two-step graph is used.

## 5. The regular order-five graph is exactly two icosahedra

**Theorem.** G_3 consists of two isolated vertices. G_5 consists of two disjoint icosahedral graphs, one on each orientation-parity sheet. At order five every vertex has H=15, |Aut|=5, M=25; each of its five distinct neighbours has presentation weight two, and its identity-return weight is fifteen. The ordinary-isomorphism quotient, and hence the converse-merged quotient, has only one vertex.

### 5.1 Structural classification of regular five-tournaments

For a regular five-tournament T, define f(v) to be the outneighbour that beats v's other outneighbour. This map is injective. If f(u)=f(v)=w, orient u->v without loss; then N^+(u)={v,w}, so f(u)=w forces w->v, contradicting v->w. There are no fixed points or two-cycles because f(v) is an outneighbour. An injective permutation of five elements without one- or two-cycles must be a five-cycle.

The edges v->f(v) form that cycle. The remaining pairs form its chord five-cycle, with one outgoing and one incoming chord per vertex. There are only two possible directed chord orientations, f^2 or f^3. The f^3 choice would make f(v) lose to the other outneighbour f^3(v), contradicting the definition of f. Therefore

    v -> f(v), f^2(v)                                  (10)

describes T. Conversely every five-cycle f gives such a tournament and is recovered intrinsically. Thus there are 4!=24 labelled regular tournaments, all isomorphic. Its automorphisms commute with the intrinsic f, so Aut(T)=<f>, of order five. The even relabellings form two orbits of size |A_5|/5=12, distinguished by epsilon.

### 5.2 Why adjacency is exactly six-edge difference

Every nonidentity return has at least six changed arcs by Section 3 and at most 2(n-1)=8. It has even support. Eight changed arcs between regular five-tournaments are impossible: their complement pair would disagree in two arcs, and a nonempty balanced directed graph has at least three arcs. Therefore the support is exactly six.

A balanced six-edge oriented graph on five vertices decomposes into two directed triangles; a directed six-cycle needs six vertices, and no other partition of six into cycle lengths at least three is possible. The triangles are edge-disjoint and meet in exactly one vertex.

For the carousel R with i->i+1,i+2 modulo five, the five directed triangles are

    013, 023, 024, 124, 134.

Exactly five pairs meet in one vertex:

    (013,024), (013,124), (023,124), (023,134), (024,134).

These form one orbit under the rotation i->i+1. The legal pair

    P=(0,2,3,4,1),       Q=(0,3,4,2,1)                  (11)

realizes one of these reversals; rotating it realizes all five. Thus every potential six-edge neighbour is an actual return, and each vertex has exactly five neighbours. No large enumeration is required.

### 5.3 An explicit icosahedral certificate

Encode a tournament by the ten bits [i->j] for pairs i<j in lexicographic order, with the first pair as the low bit. One parity component is the following twelve states:

    N=947, S=76,
    U=(613,217,710,341,668),
    L=(313,682,355,410,806).

Their adjacency is exactly

    N ~ every U_i,           S ~ every L_i,
    U_i ~ U_{i-1},U_{i+1},   L_i ~ L_{i-1},L_{i+1},
    U_i ~ L_i,L_{i+1},                              i mod 5. (12)

For example the neighbour criterion is simply that the XOR of the two codes has six set bits. The codes are obtained from R and (11) by rotations and converse; all are regular and on the same sign sheet. They are distinct, so they exhaust its twelve states. Equation (12) follows by checking the representative U_0 row and the poles, then rotating indices; it is the standard two-pentagon construction of the icosahedral graph. Its connectivity proves that this whole sign sheet is one component. An odd relabelling transports it to the other sheet.

The three directed Hamiltonian paths of R starting at zero, and their matching-endpoint paths after complement, are

| P | Q in K_P(R), with the same endpoints |
|---|---|
| (0,1,2,3,4) | P only |
| (0,2,4,1,3) | P only |
| (0,2,3,4,1) | P; (0,3,4,2,1); (0,4,2,3,1) |

Rotation gives H=5*3=15 and M=5*(1+1+3)=25. The five neighbours are a single stabilizer orbit, so the ten nonidentity presentations distribute equally, two per neighbour. The path-orbit mass is H/|Aut|=3; the pair-presentation orbit mass is M/|Aut|=5. They are different objects.

Converse swaps N with S and each state with its antipode in (12). It does not switch the sign sheet, since binomial(5,2)=10 is even. This describes the derived graph only; the original tournament metagraph is not claimed to be an icosahedron.

## 6. Exact local H variation through endpoint kernels

For a digraph D, a vertex x, and a set A not containing x, let

    L_D(A,x) = number of Hamiltonian paths of D[A union {x}] ending at x,
    R_D(x,A) = number of Hamiltonian paths of D[A union {x}] starting at x.

Both equal one when A is empty. For x,y not in W, define

    B_D(x,y;W) = sum_{A subset W} L_D(A,x) R_D(y,W\A).    (13)

It counts choices of a prefix ending at x and a suffix beginning at y, covering W outside the two endpoints. In particular it counts completions of a prescribed consecutive directed block whose first/last vertices are x/y, once its internal vertices are removed from W.

### Single-arc formula

If a->b in T and T' reverses just this arc, W=V\{a,b}, then

    H(T')-H(T) = B_T(b,a;W)-B_T(a,b;W).                  (14)

Each path that changes validity uses the flipped adjacency. The subtournaments in (13) never contain both a and b, so the same kernel is valid before and after reversal. This is an exact endpoint-refined formula, not an inference from the scores.

### Directed-triangle formula

Suppose a->b->c->a. Delete all three internal triangle arcs to form the digraph D, and set W=V\{a,b,c}. If T^triangle reverses that triangle, then

    H(T^triangle)-H(T)
      = sum_cyclic(a,b,c) [ B_D(b,a;V\{a,b}) - B_D(a,b;V\{a,b}) ]
        + sum_cyclic(a,b,c) [ B_D(c,a;W) - B_D(a,c;W) ].  (15)

**Proof.** A Hamiltonian path uses zero, one, or two triangle arcs; it cannot use all three directed cycle arcs. The zero-arc paths are exactly those of D and cancel. Exactly-one-arc paths give the first sum: every other arc must be in D. Any two triangle arcs form one consecutive three-vertex block. The old blocks are abc,bca,cab; the new ones are cba,acb,bac. Their first/last endpoint kernels give the second sum. These classes are disjoint and exhaustive. This proves (15) in every order, with no regularity assumption and no hidden third-order term.

The degree vector and directed-triangle count are preserved by the switch, but none of the endpoint terms in (15) is fixed by that fact. This identifies the actual missing coordinate behind H changes within a score fibre. A regular two-step return can likewise be evaluated by applying (14) successively to its at least six changed arcs; that telescope is exact but its intermediate single-arc states need not be regular.

## 7. What does and does not follow for the spectrum

The inherited regular order-seven triangle graph has raw values

    175 -- 171 -- 189

and nonzero directed increments +/-4,+/-18. Their gcd is two. The regular state graph is connected under triangle reversals, but 173 is absent even though it lies inside that integer lattice and between attained values. Thus neither score-fibre connectivity nor parity/gcd compatibility supplies an intermediate-value theorem. A construction needs actual legal response increments and stateful overlap, not just their gcd.

All tournament H values are odd; Aut(T) has odd order and acts freely on Hamiltonian paths. Thus |Aut(T)| divides H(T), and H/|Aut(T)| is odd. Vertex transitivity further forces n to divide H(T) by uniform starting-vertex counts. Regularity alone does not: the inherited regular seven-tournament with H=171 is a counterexample to 7-divisibility. None of these elementary congruence facts yields the global ordinary spectrum.

The conjecture that the global ordinary spectrum is exactly the positive odd integers except 7 and 21 remains OPEN. Its two omissions are established, and the current source audit is in `spectrum_basins_20261020e/spectrum.md`. The derived G_n graph is a useful regular transport object; it loses raw path multiplicities under isomorphism and has no asserted equivalence with the black merged metagraph.

## Sources and reproducibility

* H. J. Ryser's classical triangle-reversal theorem is reproduced as Theorem 1 in Guiduli, Gyarfas, Thomasse and Weidl, *2-Partition-Transitive Tournaments*, JCTB 72 (1998), 181-196: https://users.renyi.hu/~gyarfas/Cikkek/85_2partitiontranstourn.pdf . The current repository already includes the elementary proof.
* H. W. Willie Wong and E. G. Tay, *Arc reversals of cycles in orientations of G vertex-multiplications*, Discrete Mathematics 345 (2022), 112999: https://doi.org/10.1016/j.disc.2022.112999 . Its introduction explicitly states Ryser's score-sequence theorem. No stronger labelled cycle-length equivalence is imported here.
* Chan and Pak, *Computational complexity of counting coincidences*, Theoretical Computer Science 1015 (2024), 114776, Section 6.5: https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf . This is the published OPEN spectrum conjecture, not a result proved by the present transport.

Run `python3 check_regular_transport.py` and `python3 -O check_regular_transport.py` from this directory. The checker independently compares literal path enumeration against the endpoint formula for every triangle in every labelled tournament of orders three through five, checks the order-five structural certificate and weighted returns, and checks (15) on each triangle of the three supplied regular seven representatives. The root session separately reproduced all order-three/five return edges and the two-icosahedra certificate. No order-seven return-graph census or all-order component classification is claimed.
