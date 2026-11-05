# A two-sided physical gateway: H149/H151 escapes, exact wrapper calculus, and infinite continuation rank

Session 197, 2026-11-05. **All-order physical constructions and a sharp linear-interface boundary.** The unrestricted Hamiltonian-path spectrum remains open. The finite calculations below certify specified graphs and specified rightward languages; no finite support interval is extrapolated.

## 1. Change the direction, not the number of exceptional vertices

The Session197 two-port construction retains two ordered actual vertices. A rightward birth adds v beating a chosen subset of the two ports, with every other old vertex beating v; one of the three possible ports is then forgotten, without deleting its vertex. Its exact nonnegative cover interface has 27 coordinates and minimal scalar linear continuation dimension 17.

Now also permit its **leftward counterpart**: the new vertex is beaten by a chosen subset of the two ports and beats every other old vertex. As before, forget any one of the three candidate ports. Thus a leftward birth has at most two old in-neighbors, whereas a rightward birth has at most two old out-neighbors. These are actual tournament operations; they preserve every old arc.

This very small-looking addition has two different effects:

- A leftward birth escapes a concrete obstruction that survives arbitrary-depth rightward growth from every possible initial exposed pair.
- Repeated source-side and sink-side composition admits an exact all-order path-cover differential operator, but its continuation language has **infinite** rational linear dimension. The 17-dimensional state cannot merely be augmented by a fixed number of counters to implement every two-sided word.

The nearest inherited mechanisms are the [S195 endpoint-typed path-cover algebra and its arbitrary-context rank obstruction](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md), and the [S196 clone-run matching theorem](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md). The present deductions concern this explicitly restricted singleton-birth language and a particular reusable two-vertex wrapper. No literature-priority claim is made.

## 2. Strict directional escapes at both H149 and H151

Use the named S196 primitive input

    T=[122,112,67,102,108,4,32], H(T)=49.

For **every** choice of two distinct initial vertices as the rightward ports, the complete arbitrary-depth rightward descendant spectrum through 151 is contained in

    {49,85,91,115,137,147}.                            (1)

In particular neither 149 nor 151 occurs. This statement allows arbitrary later port changes and does not insist that the original root5 remain exposed.

The 21 unordered initial pairs, with their exact capped state counts, are recorded in `rightward_certificates.json`. The file contains the complete state sets and physical initial cover vectors, not only numerical summaries. Port ordering does not affect scalar reachability: maintain the bijection between the two orderings, translate the chosen subset and forgotten old position at each birth, and update that bijection on the surviving ports. The same actual vertex word and graph then result. Thus the21 unordered pairs cover all 42 ordered choices.

For clarity, the rightward state consists of unordered path covers with at most one non-port initial vertex, encoding each component by its initial and terminal port label, or a non-port label. There are at most three components and 27 endpoint signatures. To insert a rightward-born vertex, select at most one old component ending at a vertex that beats it and at most one old component starting at a selected port, excluding using the same component on both sides (which would make a cycle). Append, prepend, merge through the new vertex, or leave it singleton, then forget and relabel the chosen port. Discard a result with two non-port starts. Conversely deleting the new vertex from a counted output cover creates only an old-port start, so this rule recovers exactly the counted old covers. These inverse operations prove the universal nonnegative transfer represented by the included matrices.

Completeness uses these rightward nonnegative transfer matrices. Cap every coordinate at 152. Capping commutes with these matrices, and a capped H-readout at most151 is exact. Old Hamilton paths inject into enlarged paths, so a state with H>151 can be discarded forever. For each initial pair the verifier checks that its listed finite set contains the actual capped seed, is closed under all permitted outputs with H<=151, and consists entirely of states reachable from that seed. This is an arbitrary-depth proof, not a graph-order cutoff. A separate C++ program independently reproduces every sorted state set in full.

In contrast, choose the two old vertices 1 and 4 as leftward ports. Add one vertex that loses to those two and beats every other old vertex. The literal output is

    [122,240,67,102,236,4,32,109], H=151.              (2)

The new vertex has out-neighbor mask 109, namely `{0,2,3,5,6}`, and exactly the two in-neighbors `{1,4}`. Independent subset DP checks 151. An exhaustive check of all 128 possible one-vertex extensions of this named input finds no H149 extension and finds exactly this out-neighbor mask for H151.

H149 also has a short physical leftward construction. Start from ports `(0,1)`, add a universal source7, and forget old port1. With ports now `(0,7)`, add vertex8 that loses only to old port0 and beats every other old vertex, then retain both ports. The intermediate source addition preserves H49. The final graph is

    [378,112,67,102,108,4,32,127,254], ports=[0,7],
    H=149.                                           (2a)

The operation word, using a leftward mask for the old vertices beating the new vertex, is `[(L,0,1),(L,1,2)]`. Both subset DP and literal permutations verify (2) and (2a), and both final graphs are strong. `escape149.json` records the entire physical construction.

Thus both literal values `3*49+2` and `3*49+4` are actual induced extensions of the same old graph, although their displayed constructions choose different initial port pairs. This is not a claim that both children return from the same marked interface, and not a universal per-input digit theorem. Directional asymmetry, rather than the old induced core itself or insufficient freedom to choose an initial pair, is responsible for this particular rightward obstruction. General leftward/rightward spectrum coverage remains open.

## 3. The cyclic wrapper uses only two singleton births

For any nonempty tournament T, define W(T) by adding u,v with

    u -> T -> v -> u.                                (3)

This is physically available in the two-sided two-port language:

1. add a universal source u by a leftward birth with no exceptional in-neighbors, forgetting the first old port;
2. add v by a rightward birth whose only old out-neighbor is u, forgetting the remaining original port.

The final ports are `(u,v)`. Repeat the same two moves to iterate W indefinitely. The initial ports of T can be any distinct vertices; no special scalar or endpoint hypothesis is needed. The actual wrapper graph is also defined when T has one vertex, even though the two-port machine itself starts at order2.

Every W(T) is strong: u reaches every old vertex, every old vertex reaches v, and v reaches u.

If tau(T) is the minimum number of vertices whose deletion makes T transitive, then

    tau(W(T))=tau(T)+1,       T nonempty.              (4)

For an upper bound, delete either u or v and a minimum feedback vertex set of T. For a lower bound, if any T vertex remains, keeping both u and v would leave the directed triangle `u,t,v`, so at least one outer vertex must be deleted as well as enough old vertices to make the surviving old tournament transitive. If all of T is deleted, the cost is `|T|>=tau(T)+1`. These alternatives prove (4). Consequently repeated wrappers give strong tournaments with unbounded feedback vertex number, not a fixed-apex family.

## 4. Exact all-cover differential operator

Let

    P_T(z)=sum_(k>=0) f_k(T) z^k

be the ordinary unordered directed path-cover polynomial. A cover has nonempty path components, no cycles, and uses every vertex. The empty tournament has `P=1`. For nonempty T, `f_1(T)=H(T)`.

**Theorem1 (wrapper calculus).** For every tournament T,

    P_(W(T))(z)
      = z(z+1) [P_T''(z)+P_T(z)]
        + (2z^2+3z) P_T'(z).                         (5)

This is an actual graph identity. Its derivative notation is a compact encoding of positive fragment counts, not a formal signed construction.

**Proof.** Fix one actual old cover with k distinct path fragments `A_1,...,A_k`. Maximal old runs in a new cover can only connect through u,v. In the fragment graph, the allowed edges are

    u -> A_i, A_i -> v, v -> u.

A fragment cover can select at most one edge of each of these three types. Selecting `u->A_i`, `A_j->v`, and `v->u` simultaneously is legal exactly when i and j differ; if they coincide it forms a directed cycle, which must be excluded.

The numbers of allowed selections with 0,1,2,3 edges are respectively

    1, 2k+1, k^2+2k, k(k-1).

These have `k+2`, `k+1`, k, and `k-1` components. Therefore the contribution of the old cover is

    z^(k+2)+(2k+1)z^(k+1)+(k^2+2k)z^k+k(k-1)z^(k-1). (6)

The last term is zero for k=0 or 1. Cut a new cover at every old/new edge to recover exactly its old cover and its selected acyclic fragment edges. Conversely concatenate the chosen fragments. These maps are inverse and retain named components, so there is no missing factorial and no cycle contribution. Summing (6) over old covers gives (5). QED.

In particular,

    H(W(T))=3 f_1(T)+2 f_2(T).                       (7)

The scalar first-wrapper formula is the inherited `3c_1+c_2` rule, where `c_k=k!f_k` is the ordered k-cover count. The full polynomial operator and its repeated restricted-context consequences are what are used here.

## 5. Explicit arbitrary-order strong families and matching geometry

Put `T_(m,n)=W^m(Tr_n)`, with m>=0,n>=1. The initial transitive chain has

    P_(Tr_n)(z)=sum_k S(n,k) z^k,

because each unordered set partition determines one increasing path in each block. Repeated use of (5) constructs every graph and counts its covers exactly.

The first scalar kernels in ordered cover coordinates are

    H(W(T))   = 3c_1+c_2,
    H(W^2(T)) = 15c_1+11c_2+2c_3,
    H(W^3(T)) = 123c_1+133c_2+52c_3+6c_4.             (8)

Consequently, for every positive n,

    H(T_(1,n)) = 2^n+1,
    H(T_(2,n)) = 2*3^n+5*2^n-1,
    H(T_(3,n)) = 6*4^n+28*3^n+13*2^n-11.             (9)

These follow from the proved operator and the Stirling expansion, not fitted values. For m>=1 all these graphs are strong and have feedback vertex number exactly m by (4).

There is a geometric explanation of the growing highest exponential base. Collapse the initial transitive module to a vertex r and label the new pair at wrapper step i as `(u_i,v_i)`. The out-neighbors of r are the v_i, and its in-neighbors are the u_i. Their cross arcs satisfy

    v_i -> u_j iff i=j.                              (10)

For i<j the later source u_j beats the earlier v_i; for i>j the later sink-side vertex v_i loses to the earlier u_j. Thus the root-triangle incidence graph is exactly a matching of size m. The S196 matching theorem gives maximal clone-run count m+1. At that maximum, the exterior gaps are the m matched pairs in arbitrary order, giving coefficient m! for the ordered `(m+1)`-cover count. In particular

    H(T_(m,n)) = m! (m+1)^n + O(m^n)                 (11)

for fixed m>=1 as n tends to infinity. This asymptotic is only a growth description, not a primality or spectrum-support claim.

## 6. Infinite linear continuation rank already in this tiny two-sided language

**Theorem2.** No fixed finite-dimensional rational state with homogeneous linear generator updates and a linear H-readout can implement every two-sided word from the two-vertex seed. The obstruction already uses neutral sink prefixes and the two-birth wrapper W.

Let `L` denote the polynomial operator in (5). The coefficient of `f_(m+1)(T)` in `H(W^m(T))` is

    m!(m+1)!,                                       (12)

and no `f_k` with k>m+1 occurs. To reach the coefficient of z from degree m+1 in m applications of (6), every step must lower the degree by one. Multiplying the coefficients `k(k-1)` for k=2 through m+1 proves (12). Therefore in ordered cover coordinates the continuation kernel is triangular with leading coefficient m!.

Fix any N. Use as prefixes the transitive inputs `Tr_2,...,Tr_(N+1)`, all reachable from the two-vertex seed by universal rightward-sink births that retain the original ports. Use as suffixes `W^0,...,W^(N-1)`, each a fixed physical word. The resulting N-by-N continuation matrix C satisfies

    det C = N! product_(j=1)^N j!(j-1)! !=0.          (13)

To see the extra N! factor caused by starting the inputs at order2, let `S0_(i,j)=S(i,j)` for 1<=i,j<=N. This is lower triangular with diagonal 1. The matrix `S1_(i,j)=S(i+1,j)` equals S0 times the upper-bidiagonal matrix with diagonal 1,2,...,N and superdiagonal 1, by the Stirling recursion. Hence `det S1=N!`. Multiplication by the factorial column weights and the triangular wrapper kernels gives (13).

Every d-dimensional homogeneous linear machine would factor C through a d-dimensional space, giving rank at most d. Since N is arbitrary, no finite d works. This is a statement about the specified linear representation model; it does not prohibit unbounded-degree polynomials, nonlinear encodings, exact graph algorithms, or a more restricted two-sided sublanguage.

This boundary explains why merely adding the converse of a missing birth operation does not preserve a fixed continuation interface. For the wrapper sublanguage, the full polynomial (5) gives an exact replacement, but its degree and the number of relevant cover counts grow. No transfer on the untyped polynomial alone is claimed for arbitrary port-dependent leftward/rightward words. The existence of that replacement does not solve the inverse target problem.

## 7. Certificates and reproduction

First run `python3 build_rightward.py` in this directory to regenerate the large state-set and seed files; see `REPRODUCE.md`.

The packet is self-contained. `rightward.py` and the 27-state serialized matrices are unchanged copies of the independently audited one-sided construction, with their provenance stated above. `core.py` implements the new wrapper and exact polynomial operator.

    g++ -O3 -std=c++17 census.cpp -o /tmp/s197_twosided_cap
    /tmp/s197_twosided_cap rightward_seeds.txt > independent_cpp.json
    python3 verify.py --cpp independent_cpp.json > receipt.json
    python3 -O verify.py --cpp independent_cpp.json > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The verifier checks (5) against literal unordered-cover enumeration for every labelled input of orders 0 through 4, including the empty case: 76 actual graphs and their wrapped outputs. It constructs the 5-by-5 continuation matrix by independent vertex-subset Hamilton counting and checks determinant 1,194,393,600. It checks the explicit scalar families, the central matching graph, the exact feedback-number increment on all 75 nonempty small inputs, both child graphs by subset DP and permutations, and all 128 one-vertex extensions of the named H49 input. For all 21 initial rightward port pairs it independently reconstructs physical cover profiles using a different component-order convention, validates entire finite closure certificates, and compares the complete sorted state sets with C++.

The all-order formulas, feedback-number statement, and infinite-rank theorem follow from the stated bijections and determinant factorization, not from these controls. The inverse target problem and unrestricted spectrum coverage remain unresolved by this work. No dirty checkout or worktree was modified.
