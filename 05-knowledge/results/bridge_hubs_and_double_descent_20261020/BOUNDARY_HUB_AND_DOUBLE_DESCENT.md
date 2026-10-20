# Boundary hubs, double descent, and bridge-cut restrictions

Date: 2026-10-20. Status: PROVED after independent recovery review; see [the independent audit](AUDIT.md).
The frozen order-nine packet is not modified.

## 1. Graph reconciliation

Q here is the **paired-black** multigraph: reflected black complement lines are identified in pairs, loops and remaining parallel multiplicities retained. K is the literal black multigraph, which doubles every edge of Q. F is K plus the actual blue lines. The incoming packets sometimes use Q for the blue quotient obtained after contracting K components; that is NOT our Q.

At order nine the two independently computed censuses agree: 97,144 merged vertices, 97,138 black-active vertices, one active black component and six pure-blue isolates. Our Q has seven bridges, all Q-leaf edges. F has five bridges, the five pure-blue leaf edges. The actual blue subgraph has 24 bridges, each incident to a blue-degree-one node. Nineteen of these are mixed-to-mixed and have no endpoint of full F-degree one. Thus the incoming phrase '19 non-leaf blue-only bridges' measures leaf status in F, not in the blue subgraph. There is no numerical contradiction.

Exact incoming sources read using git show, without checkout edits:
- Commit 126830b14ba947b188f5fa6f1c7a080ffe636a82, [order-nine census](../blue_quotient_20261020b/census.md).
- Commit 7d21b9fdcfaa30e723b92a51ebe02d35ff57b250, [cycle-changing descent](../black_descent_20261020c/cycles.md), [SCC terminals](../black_descent_20261020c/scc.md), and [original boundary hubs](../defect_transitions_20261020/blue.md).
All source paths are under 05-knowledge/results/. Exact commit, Git-blob, byte-count and SHA-256 pins are in [DEPENDENCIES.json](DEPENDENCIES.json).

The following arguments use their proved all-order black descent and prime-orbit block lemma. They do not infer connectivity or bridge-freeness from descent.

## 2. A common black hub certifies Q nonbridges

Let U be a nonstrong tournament of order n>=3. Suppose two black Hamiltonian paths P,P' have the same ordered first, second, and last vertices (a,b,z), and represent distinct Q edges incident to [U]. Write d for outdegree in U. Assume

    d(a)+d(b) != n,
    delta = 2d(a)-d(b)-d(z)-3 != 0.                 (1)

Then BOTH Q edges lie on a cycle of length at most four.

Proof. Let K_P complement all nonpath arcs while retaining labels, and put S=K_P(U), S'=K_P'(U). Both are strong: a nonstrong Hamiltonian path starts in the first SCC and ends in the last, so U has a->z; complementing that nonpath arc closes the preserved path into a directed Hamiltonian cycle.

Rotate the preserved path in S to (b,...,z,a), and do the corresponding rotation in S'. Both are directed because the preserved path closes. The two resulting complement partners are the same labelled tournament

    V = U with ab and az reversed.                 (2)

The rotation endpoints in V are b,a and have score sum d(b)+1+d(a)-2=d(a)+d(b)-1. The first condition in (1) rules out a positional anti-automorphism, so both rotation edges are black.

The score changes U to V are -2 at a and +1 at b,z. The cyclic-triangle identity c3=binom(n,3)-sum binom(d,2) gives c3(V)-c3(U)=delta. Thus [V]!=[U], even modulo converse. Also [S],[S'] differ from [U] because they are strong. The walk

    [U] --e'-- [S'] --black-- [V] --black-- [S]

avoids e: its first edge is explicitly different, and neither remaining edge can equal e because none has endpoint [U]. This is an alternate route of length at most three. Vertex identifications or additional parallelisms may shorten the cycle; they cannot reintroduce e. Interchanging P and P' proves the claim for e'. QED.

This is not the arc-disjoint commuting-square criterion. The two source paths may share almost all their arcs. The changed triangle count at the common hub is the additional anti-folding coordinate.

### General shared-edge version

It suffices instead that P,P' share endpoints a,z and a directed path edge u->v, represent different incident Q edges, and satisfy

    d(u)+d(v)-[u=a]+[v=z] != n-1,
    -sum_i d(i)q_i - (1/2)sum_i q_i^2 != 0,
    q = -e_a+e_z-e_u+e_v.                            (3)

Cut both partner cycles immediately after u->v. The common hub is U with az and uv reversed. The first expression is the hub's endpoint score sum, and the second is its cyclic-triangle change. The preceding proof is otherwise unchanged. The endpoint pairs az and uv are distinct because n>=3 and uv is a consecutive path edge.

## 3. Uniform nonbridge theorem for rigidly guarded nonstrong nodes

Call an ordered SCC word rigid if every letter is a singleton or a directed triangle. Let

    U = A -> B -> D,

where A and D are nonempty rigid words and B is any tournament, possibly empty. Assume max(|A|,|D|)>=2. If deg_Q([U])>=2, EVERY incident Q edge lies on a cycle of length at most four. Consequently any Q bridge incident to such a U is a leaf edge.

Proof. Taking converse if necessary makes |A|>=2. Every singleton/C3 component has one Hamiltonian-path orbit under its automorphism group. Independent automorphisms of the outer components therefore normalize the outer restrictions of EVERY Hamiltonian path to fixed literal paths in A and D. In particular every incident Q edge has a representative with one common boundary triple (a,b,z). Choose any other incident Q edge; nonstrong U has no loops because all partners are strong.

If B also consisted only of rigid SCCs, U would have exactly one ordinary Hamiltonian-path orbit. Its merged Q degree would be either 0 (self-converse palindrome) or 1 (nonpalindrome). The degree hypothesis therefore forces B to contain a nonrigid strong SCC, whose order is at least four. Hence n>=|A|+|D|+4>=7.

The first vertex a lies in a first singleton or C3, so d(a)>=n-2. Since b lies within the rigid prefix, d(b)>=n-3; also d(b)<=n-2 because a->b. The last vertex lies in a last singleton or C3, so d(z)<=1. Therefore

    d(a)+d(b) >= 2n-5 > n,
    2d(a)-d(b)-d(z)-3 >= n-6 > 0.

Section 2 applies to the chosen pair of edges. Since the first edge was arbitrary, all incident edges are nonbridging. QED.

This includes infinitely many non-self-converse and mixed-self-converse nodes, with arbitrary middle tournament B. The case 1->B->1 is deliberately excluded unless B itself supplies further rigid outer SCCs. The theorem does not address edges with two strong endpoints.

## 4. Every strong non-self-converse node has two descending Q edges

**Theorem.** At every order n>=3, every strong non-self-converse tournament has at least TWO distinct Q edges whose other endpoint has strictly fewer cyclic triangles.

All its Hamiltonian paths are black. A directed Hamiltonian-cycle cut beginning at v_j and ending at v_(j-1) descends precisely when

    d(v_j) <= d(v_(j-1)).                              (4)

At least one cut satisfies (4), since consecutive score differences sum to zero.

Let G<=Z/nZ be the subgroup of rotations of a chosen cycle that are automorphisms of T. Two cut paths give the same descending Q edge if and only if they differ by a rotation in G. Indeed, within a non-self-converse source ordinary path-orbits are the only possible identifications: reflection lands in the distinct converse fibre, and complementation lands at a strictly smaller c3. Any automorphism mapping one complete cut word to another is precisely the corresponding cycle rotation.

If T is regular, every cut descends. A single descending Q edge would therefore force every rotation to be an automorphism, so T would be circulant. Negation of the cyclic labels is an anti-automorphism of every circulant tournament, contrary to non-self-converseness. Thus regular T has at least two descending edges.

Suppose T is nonregular and, for contradiction, has only one descending Q edge in the entire graph. On the chosen cycle all descending cuts then form a single coset of G. Put g=|G|. Tournament automorphism groups have odd order: an involution would interchange some pair of vertices and reverse their arc, and Cauchy's theorem excludes every even-order group. Thus g is odd. If g=1, the remaining n-1 cut positions are strict score increases. Cutting the sole descending boundary produces n distinct increasing scores, impossible because a strong tournament has all scores in {1,...,n-2}.

Choose an odd prime p dividing g. The rotation by n/p is an automorphism of order p. Its vertex orbits have constant scores.

The case ell=n/p=1 is impossible here: the rotation is then transitive and forces T to be regular, contrary to the standing nonregularity assumption.

### At least three prime orbits

If ell=n/p>=3, apply the incoming prime-orbit block lemma ([cycles.md](../black_descent_20261020c/cycles.md), Section 3). It gives a directed Hamiltonian cycle in which each p-vertex orbit occurs consecutively. Nonregularity gives at least two orbits with different scores. Cut an internal edge in each of these two blocks. Both resulting paths have equal endpoint scores and hence descend by exactly one triangle. Their endpoint scores differ between the two paths, so no automorphism can identify them. They give two distinct Q edges, a contradiction.

### Exactly two prime orbits

Here n=2p and the rotation by two is an automorphism. Write the alternating original cycle as U_0,W_0,U_1,W_1,... with subscripts in F_p. Each induced orbit is a circulant tournament; let A,B be its respective sets of forward increments. Both have size (p-1)/2. Write

    S={d : U_i -> W_(i+d)}.

This is independent of i. The original alternating cycle gives 0 in S and -1 not in S, so S is nonempty and proper. The two orbit scores are distinct because T is nonregular.

If there are alpha in A and beta in B with alpha+beta!=0, some x in S has x-alpha-beta not in S. Otherwise S would be invariant under a nonzero translation, impossible in F_p. Put t=x-alpha. Then

    U_0,U_alpha,...,U_(-alpha),
    W_t,W_(t+beta),...,W_(t-beta)

is a directed Hamiltonian cycle: the two block-boundary tests are t+alpha=x in S and t-beta=x-alpha-beta not in S. Cutting one internal edge in each block again gives two inequivalent equal-score descending paths.

If no such alpha,beta exist, A and B must both be singletons. Thus p=3 and B=-A. Every nonempty proper S subset F_3 is invariant under a reflection d->2c-d (take c to be its singleton element, or the singleton element of its complement). The map

    U_i -> W_(i+c),    W_j -> U_(j-c)

is an anti-automorphism: the two induced directed triangles have opposite increments, and the cross-arc test is exactly S=2c-S. This makes T self-converse, again a contradiction. The two-orbit case is exhausted. QED.

The theorem is genuinely stronger than existence of descent, but its non-self-converse hypothesis is essential. Strong self-converse examples at n5 have Q degrees 4 and 3 but only one descending edge; at n6 there is such a node of Q degree 6. Thus no universal two-exit claim is made.

## 5. A new restriction on bridge cuts

Let e=uv be a Q bridge and delete it from its active component. Suppose one of the resulting sides A contains no nonstrong tournament.

Then:
1. A's endpoint u is its UNIQUE minimum of c3.
2. c3(v)<c3(u), and e is the ONLY strictly descending Q edge at u.
3. u is self-converse.
4. The opposite side contains a nonstrong tournament.

Proof. Every vertex of A is strong and black-active, so the incoming all-order black-descent theorem gives a strictly lower neighbor. At a minimum w of c3 in A, that neighbor must be outside A. The only boundary edge is e, so w=u and that descending edge is e. Any second minimum would have the same contradiction without access to e. Every descending edge at u must leave A and is therefore e. Section 4 rules out a non-self-converse u. Starting at v and repeatedly descending cannot cross e upwards to u, so termination supplies a nonstrong vertex on the opposite side. QED.

Equivalently, any bridge of equal endpoint potential has nonstrong terminals on BOTH sides. So does any bridge whose higher-potential endpoint is non-self-converse. A hypothetical terminal-free bridge branch is forced to hang from a self-converse single-descent bottleneck.

This does not prove that every bridge meets a leaf. A bridge could separate two nonstrong terminal basins, or have the indicated self-converse bottleneck. Contracting an arbitrary descent forest preserves connected components but can erase bridges; it is not a valid bridge-classification shortcut.

## 6. Remaining obligations

- Completed in the independent recovery audit: Q-edge orbit identification, prime-orbit block lemma, two-prime-orbit case, and bridge-cut consequence.
- Completed finite controls: literal shared-edge hub checks through order seven and replay of the rigid-guard certificates through order eight; these do not replace the all-order proof.
- Resolve paired-black Q short-cycle coverage for singly padded 1->B->1 and general nonrigid endpoint SCCs.
- Exclude or construct bridges separating two terminal basins.
- Classify strong self-converse single-descent bottlenecks; existing examples prohibit assuming they are leaves.

All-order black-core connectivity, universal short-cycle coverage, and the vanishing of the nonleaf-bridge series remain open.
