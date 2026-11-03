# Shared-root diamonds: exact unit-distance compilation, rigidity, and the loss of a unit spine

**Status.** Proved elementary classification and constructions, with exact independent finite controls. The seven-vertex member is the classical Moser spindle, not a new four-chromatic graph. Literature priority is not asserted for the longer family or its properties. No improvement to the Hadwiger–Nelson bounds or planar unit-distance edge maxima is claimed.

The concrete bridge has two parts. First, a familiar local module admits a complete classification of which abstract gluing graphs have a genuine injective Euclidean realization. Second, the Moser module occurs literally as an induced subgraph of the single-arc-reversal metagraph of strong six-vertex tournaments. This metagraph is **not** the repository's complement-line/converse-merged metagraph.

## 1. The rooted diamond compiler

For a finite simple graph H with m vertices, form D(H) as follows. There is one shared root O. For each v in V(H), introduce three vertices t_v,a_v,b_v and the five edges

    O a_v, O b_v, a_v b_v, a_v t_v, b_v t_v.

Thus O and t_v are the nonadjacent terminals of a diamond K4 minus an edge. Also include t_u t_v for every uv in E(H). No other edges are specified. A unit realization here is injective: all distinct graph vertices occupy distinct points. A faithful realization additionally has no unlisted unit distances. Crossings of drawn edges are allowed.

**Theorem 1.** D(H) has an injective planar unit realization if and only if H is a linear forest (a disjoint union of paths, with isolated vertices allowed). Whenever it has one, it has a faithful one.

**Proof, necessity.** The points O and t_v are the two distinct intersections of the unit circles about a_v and b_v. Since |a_v b_v|=1, their distance is

    |O t_v| = sqrt(3).

All tips therefore lie on one circle of radius sqrt(3). Two points of this circle are a unit distance apart exactly when their angular difference is plus or minus theta, where

    cos(theta)=5/6,  q=exp(i theta)=(5+i sqrt(11))/6.

The rotation q has infinite order. Indeed, if q were a root of unity, q+q^(-1)=5/3 would be a rational algebraic integer, hence an integer, a contradiction. Consequently the unit-distance graph on each orbit z q^Z of this circle is a two-sided infinite path. Different orbits have no unit edges. Every finite subgraph is a linear forest, proving necessity.

**Proof, sufficiency.** For a path of length m in vertices, use the coordinates in Section 2. Different path components may be rotated about O independently. There are only finitely many prohibited relative rotations producing either a collision or an unlisted unit distance between two component vertices. Avoiding that finite set produces a faithful realization of every linear forest. The same argument applies to isolated diamond components. ∎

This conclusion is stricter than the ordinary plane-wide K2,3 obstruction. The radius forced by the retained shared root rules out every tip cycle, including an even cycle of any length, and every tip of degree at least three.

**Relaxing injectivity changes the answer.** If each individual diamond remains injective but different diamond copies may occupy the same points, D(H) is realizable if and only if H is bipartite. Necessity follows from the infinite-path orbit description. Conversely, map the two sides of a bipartite component to two unit-separated tips, and superimpose corresponding diamond copies. Thus the weak constraint language admits all bipartite graphs, while the physical injective language admits only linear forests. This is an explicit example where forgetting global distinctness changes the realizability class, rather than merely changing a count.

## 2. An exact infinite family

Write G_m=D(P_m), with vertices O and t_k,a_k,b_k for 0<=k<m. Put

    O=0,
    t_k=sqrt(3) q^k,
    a_k=((sqrt(3)+i)/2) q^k,
    b_k=((sqrt(3)-i)/2) q^k,
    q=(5+i sqrt(11))/6.

All real Cartesian coordinates lie in Q(sqrt(3),sqrt(11)). The only unit distances are the five diamond edges at every k and t_k t_(k+1). Hence

    |V(G_m)|=3m+1,       |E(G_m)|=6m-1.

Here is a complete check for the absence of accidental edges. Tips have radius sqrt(3), interiors a_k,b_k have radius 1, and O is the sole center. Tip-tip unit distances require angular difference plus or minus theta, giving consecutive indices only. Two interior points are unit-separated exactly at angular difference plus or minus pi/3; irrationality of theta/pi leaves just a_k b_k. A tip and an interior point are unit-separated exactly at angular difference plus or minus pi/6; again only the two members of its own diamond qualify. The same irrationality excludes collisions among distinct interiors; different radii exclude all other collisions.

G_1 is the diamond and G_2 is the classical Moser spindle. This is a *shared-root chain*, not a claim that every construction called a chain of diamonds in the literature uses the same gluing.

**Chromatic ceiling.** chi(G_1)=3 and chi(G_m)=4 for every m>=2. In a three-coloring, the adjacent vertices a_k,b_k force both O and t_k to use the remaining color. The edge t_0 t_1 forbids that. A four-coloring is explicit: use color 4 at O, colors 1 and 2 at every a_k,b_k, and alternate colors 3 and 4 on the tips.

More generally, for every realizable D(H) with H nonempty, its chromatic number is 3 when H has no edges and 4 when H has an edge. Therefore this compiler cannot yield a five-chromatic unit-distance graph. The obstruction is a proved limit of this specific gluing language, not a limit on arbitrary Moser-spindle compositions or on the Hadwiger–Nelson problem.

## 3. Rigidity is compatible with exponentially many labeled realizations

**Theorem 2.** If H is a nonempty linear forest with m vertices and c components, every injective unit realization of D(H) has independent bar constraints and rigidity-matrix rank

    rank R = 5m+|E(H)| = 6m-c.

Its infinitesimal motion space modulo Euclidean motions has dimension c-1. In particular G_m is infinitesimally rigid and isostatic: its 6m-1 edges give rank 2(3m+1)-3.

**Proof.** Each diamond is two nondegenerate equilateral triangles sharing an edge, so it is infinitesimally rigid and contributes five independent constraints. Pin O. Each diamond then has one residual angular velocity omega_v about O, with no other internal infinitesimal motion. For a unit tip edge uv, differentiate its squared length. The result is a nonzero scalar multiple of omega_u-omega_v, since the tip angle is plus or minus theta and sin(theta) is nonzero. The tip-edge constraints are therefore a rescaled incidence matrix of H. A forest has rank m-c. The pinned motion dimension is c, and restoring translations then quotienting all three Euclidean motions gives c-1. This also proves independence of every bar row. ∎

The geometric count needs a qualification that plain rigidity does not provide. With O and t_0 fixed, successive tips admit two local angle choices, plus or minus theta. A length-(m-1) sequence therefore has 2^(m-1) locally admissible words. Any change of sign immediately repeats a tip from two steps earlier. Global injectivity leaves exactly the two constant-sign words. Reflection identifies these two point-set configurations.

Each diamond still admits an independent exchange of its labeled twins a_k,b_k. For m>=2:

- The underlying point configuration is unique up to Euclidean congruence and independent twin exchanges.
- There are exactly 2^m congruence classes of **labeled** realizations. For example, normalize O=0, t_0=sqrt(3), a_0=(sqrt(3)+i)/2. There are two choices of tip-chain direction and 2^(m-1) choices for the remaining twin labels.
- The abstract automorphism group is (C2)^m semidirect C2, of order 2^(m+1). The final C2 reverses the tip path.
- The Euclidean symmetry group of the point set has order 2. Its nontrivial symmetry reflects k to m-1-k and exchanges a with b.

For the automorphism statement, O is the unique vertex of degree 2m (also unique of degree 4 when m=2). Its nonneighbors are exactly the tips; their induced graph is P_m. After deciding path reversal, only the independent twin exchanges remain. The point-set symmetry statement follows because an infinite-order circle orbit has no nontrivial finite rotational symmetry, and the finite consecutive tip interval has exactly one reflection axis.

Thus **infinitesimal rigidity, labeled global rigidity, abstract automorphism, and actual Euclidean symmetry are four different notions** in one completely explicit family. In particular, the family is not globally rigid as a labeled framework when m>=2.

## 4. A genuine unit-spine failure at 16 vertices

Prior repository work distinguishes intrinsic traceability of the unit graph from failures caused only by the external ordering used to build a point tournament. The family G_m gives an exact, intrinsic failure while retaining rigidity and four-chromaticity.

**Theorem 3.** For m>=1:

    G_m has a Hamiltonian path  if and only if m<=4;
    G_m has a Hamiltonian cycle if and only if m<=2.

For m>=2 it has minimum degree 3 and vertex connectivity 2.

**Negative path direction.** For m>=5, delete

    S={O,t_1,t_3}.

There are exactly five remaining connected components: the triangle at index 0, the private edge a_1 b_1, the triangle at index 2, the private edge a_3 b_3, and the nonempty chain beginning at index 4. A Hamiltonian path would leave at most |S|+1=4 components after deleting S, a contradiction.

**Negative cycle direction.** For m>=3, deleting {O,t_1} leaves three components, whereas a Hamiltonian cycle could leave at most two.

**Positive directions.** The diamond has the cycle O,a_0,t_0,b_0,O. G_2 has the cycle

    O,a_0,b_0,t_0,t_1,a_1,b_1,O.

For G_3 and G_4, respectively, Hamiltonian paths are

    a_0,b_0,t_0,t_1,a_1,b_1,O,a_2,b_2,t_2;
    a_0,b_0,t_0,t_1,a_1,b_1,O,a_2,b_2,t_2,t_3,a_3,b_3.

Their listed consecutive pairs are graph edges. Deleting O leaves a connected chain of triangles. Deleting any other single vertex leaves every remaining diamond fragment connected to O. Thus there is no cut vertex. Deleting {O,t_0} disconnects the graph for m>=2, proving vertex connectivity 2. The degree formula gives minimum degree 3. ∎

The first failure has exactly **16 vertices and 29 unit edges**. It is isostatic and four-chromatic, and its failure cannot be repaired by relabeling points or changing a tournament flip convention. However it is not proposed as an edge-extremal configuration. In particular it does **not** refute the repository's conjectural traceability for planar unit-distance extremizers. It isolates why rigidity, minimum degree 3, and two-connectivity are insufficient hypotheses for such an extremal theorem.

Independent subset dynamic programming gives the following exact unoriented counts, including both positive and negative cases:

| m | Hamiltonian paths | Hamiltonian cycles |
|---|---:|---:|
| 1 | 6 | 1 |
| 2 | 48 | 4 |
| 3 | 80 | 0 |
| 4 | 64 | 0 |
| 5 | 0 | 0 |

## 5. The exact path-cover threshold, not only traceability

Let pc(X) denote the minimum number of vertex-disjoint nonempty undirected paths covering every vertex of X; singleton paths are allowed. Then

    pc(G_m-O)=ceil(m/2),
    pc(G_m)=max(1,ceil(m/2)-1).

**Proof.** After deleting O, consecutive pairs of diamonds can each be covered by the path

    a_(2j),b_(2j),t_(2j),t_(2j+1),a_(2j+1),b_(2j+1).

If m is odd, the final triangle has its own path a_(m-1),b_(m-1),t_(m-1). This uses ceil(m/2) paths. For a lower bound, delete the tips t_1,t_3,... with index strictly less than m-1. If their number is s=floor((m-1)/2), the remaining graph has 2s+1 components. Deleting s vertices from a cover by p paths produces at most p+s path components, so p>=s+1=ceil(m/2).

Deleting O from any path cover of G_m increases its number of paths by at most one. Hence pc(G_m)>=pc(G_m-O)-1, and a nonempty graph always needs at least one path. Conversely, when the displayed cover of G_m-O has at least two paths, each of the first two has an interior diamond vertex at an endpoint. Connect those two endpoints through O, reducing the number of paths by one. When there is only one path, prepend O at an interior-vertex endpoint. This attains the lower bound. ∎

This gives a literal connection to the S165 cover-interface definition, which applies to general directed graphs, not only tournaments. Replace each unit edge of G_m by its two directed arcs and call the resulting digraph X_m. For its unordered directed-path-cover coefficients p_k(X_m),

    p_k(X_m)=0 for k<max(1,ceil(m/2)-1),
    p_k(X_m)>0 at k=max(1,ceil(m/2)-1).

Consequently the r-slot set C_(X_m)(r) is empty exactly below that threshold. The asymptotic slot cost is m/2, even though the graph has one shared geometric root and every module is individually rigid. This is a proved resource bottleneck, not a statement that the bidirected graph is a tournament or that an arbitrary tournament completion retains the unit-edge path-cover counts.

## 6. A two-state chromatic interface with an exact transfer

The same rooted module has a small, genuine boundary interface for colorings. Fix the root color and let q be the number of colors. Define

    A=(q-1)(q-2),       B=(q-2)(q-3).

If a tip equals the root color, there are A ordered choices of colors for a,b. If the tip differs from the root color, there are B. Track two states E and N according as the final tip equals or differs from the root. The initial weighted row vector and transfer matrix are

    v_1=(A,(q-1)B),
    M = [[0,(q-1)B], [A,(q-2)B]].

Then the chromatic polynomial of G_m is

    P_(G_m)(q) = q v_1 M^(m-1) (1,1)^T.

This identity holds for all positive integer q by the coloring decomposition and hence as a polynomial. If H has path components of vertex sizes l_1,...,l_c, its shared-root polynomial is

    P_(D(H))(q) = q product_j [v_1 M^(l_j-1) (1,1)^T].

There is one factor q, not q^c, because the root color is shared. For Z_m=v_1 M^(m-1)(1,1)^T, Cayley–Hamilton yields

    Z_(m+2) = (q-2) B Z_(m+1) + (q-1) A B Z_m.

At q=4 this is Z_(m+2)=4 Z_(m+1)+36 Z_m, with Z_1=12 and Z_2=96. For m=2 one obtains the explicit Moser polynomial

    P_(G_2)(q)=q(q-1)(q-2)^2(q-3)(q^2-3q+4).

The transfer is an actual two-state count, not an analogy. By contrast, deciding whether arbitrary local rotation words give an injective realization also requires remembering occupied tips. On this particular path language, the freshness guard collapses the 2^(m-1) local words to two; in the entire compiler, dropping that guard enlarges linear forests to all bipartite tip graphs. This is a precise instance of the distinction between a formal interface and a physically realizable object emphasized in S165.

## 7. A literal Moser spindle inside a tournament flip metagraph

Let M_6^strong have one vertex for each isomorphism class of strong tournaments on six vertices. Two distinct classes are adjacent when reversing one arc in a representative of one yields a representative of the other. Forget loops and edge multiplicities.

Encode a labeled tournament by the 15 bits indexed by lexicographically ordered pairs (i,j), 0<=i<j<6; bit 1 means i->j. The seven canonical representatives

    (O,A,B,C,D,E,F) = (18,22,26,181,81,146,145)

are all strong and pairwise nonisomorphic. They induce precisely

    OA, OB, AB, AC, BC, OD, OE, DE, DF, EF, CF.

This is G_2 under C=t_0, F=t_1, (A,B)=(a_0,b_0), and (D,E)=(a_1,b_1). Thus its induced graph has an exact faithful planar unit-distance realization and chromatic number 4. In particular chi(M_6^strong)>=4.

The independent audit enumerates all 720 vertex permutations for each required canonicalization, checks all 15 flips of each representative, verifies all seven strong-connectivity conditions, and compares the resulting 21 adjacency decisions to exact squared distances in Q(sqrt(3),sqrt(11)). It does not assume that a canonical labeling preserves the original vertex names. The JSON contains a representative arc-flip witness for each edge.

This is an actual common finite object in the two settings. It does not imply that the whole metagraph is a planar unit-distance graph; that G_m embeds in it for m>2; that the metagraph's vertices are Euclidean points before the embedding is chosen; or that the complement-line merged metagraph has the same edges. Those would be separate assertions.

## Verification and source boundaries

Run:

    python3 check.py > results.json
    python3 -O check.py > results-optimized.json
    cmp results.json results-optimized.json
    python3 audit_metagraph.py > metagraph-audit.json
    python3 -O audit_metagraph.py > metagraph-audit-optimized.json
    cmp metagraph-audit.json metagraph-audit-optimized.json

Both pairs of runs agree byte-for-byte. Controls include:

- 6,846 exact squared-distance checks for m=1,...,12,16,24, proving the asserted full unit graph for those coordinates, with no floating-point tolerance.
- Exact modular rigidity-rank certificates for the same 14 frameworks. A rational algebraic coordinate reduction to a finite field is used only for rank, with explicitly recorded prime and square roots. A nonzero modular minor certifies a nonzero characteristic-zero minor; the general rank statement has the analytic proof above.
- Independent direct color enumeration against the transfer for m=1,2,3 and q=1,...,5.
- Independent Hamiltonian path/cycle dynamic programming through m=5, plus explicit optimal cover/separator certificates for m=1,...,24.
- All 1,099 labeled simple graphs through five vertices: the linear-forest characterization agrees with an independent injective integer-address solver. Linear-forest counts are 1,2,7,34,206; bipartite counts are 1,2,7,41,376.
- Exhaustive local sign words through m=12, leaving exactly two injective tip words in every nontrivial case.
- A separate canonical-permutation audit of the seven-class tournament spindle, with exact geometric verification.

Relevant earlier repository work:

- [S165 FI interfaces](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/fi_interfaces/FI_INTERFACES.md), for the explicit separation of shared resources, finite interfaces, and literal tournament realizability. No category-equivalence claim is needed here.
- [Unit-spine/gauge separation](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/hypotheses/HYP-2203-unit-distance-moser-unit-spine-gauge-separation.md) and [endpoint-compatible ear recursion](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/hypotheses/HYP-2620-unit-distance-spine-ear-recursion.md), for the intrinsic-vs-order-dependent traceability distinction. Their dense-extremizer conjectures are not resolved here.
- The earlier `hadwiger_dehn_niven_spindle_s685.out` identified the spindle angle cos(theta)=5/6. Only that elementary angle observation is used. This note does not use or endorse that file's broader Dehn-invariant or lattice-escape identifications.
- Aubrey de Grey, [The chromatic number of the plane is at least 5](https://arxiv.org/abs/1804.02385), is the established five-chromatic finite-unit-distance result. The present family remains four-chromatic.

The new metagraph witness was supplied by the concurrent planar-metagraph investigation; `audit_metagraph.py` is an independent verification, not a claim to have discovered the witness independently. The family classification, transfer, and rigidity/spine calculations are given with full proofs so their validity does not depend on finite extrapolation.
