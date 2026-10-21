# Port-faithful metagraph realizability and synchronized triangle actions

Session25, 2026-10-21. **PROVED:** elementary certificate and transfer statements below. **FINITE-EXACT:** independent agreement/tiling checks recorded in `check.out`. This is an exact finite wiring algorithm, not an all-order structural classification, efficient isomorphism algorithm, or non-leaf-bridge theorem. No literature priority is claimed.

## Inheritance and scope

The carrier and color conventions are those of THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md, section 0, and THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md, statement 3. The closest mechanism is Session24's block-contiguous triangle lift; the canonical hostile is C3[C3], whose complemented lifted path has H=2983, versus H=27 for the substituted complemented base. The corrected near miss is that an injection of path orbits and a black-degree inequality do not transport adjacency. The least-used sidecar is the identity of every run's entry and exit vertex.

Concept board: (1) class wiring versus path states; (2) run covers versus scalar cover counts; (3) genuine S3 actions versus six-state cardinalities; (4) loops/multiplicity versus simple graph shape; (5) exact simulation versus insufficient quotient coordinates. The test below recovers actual adjacency, not just degree statistics.

All theorems assume n>=3. There are no loops or repeated vertices in a Hamilton path. Tournament converse reverses every arc. For a Hamilton path P of A, K_P(A) preserves its n-1 arcs and reverses all other arcs.

## 1. Necessary-and-sufficient pair realization

**Theorem 1 (agreement certificate).** Let A and B be specified n-vertex tournaments. For a bijection f:V(A)->V(B), define D_f to contain precisely those directed arcs u->v of A for which f(u)->f(v) is an arc of B. Then B is isomorphic to K_P(A) for some Hamilton path P if and only if D_f is exactly a directed Hamilton path for some f.

Equivalently the agreement digraph must have n-1 arcs, every indegree and outdegree at most one, and be weakly connected. One may replace weak connectivity by acyclicity in this criterion. These are sufficient conditions, not just invariant matching conditions.

*Proof.* If f identifies B with K_P(A), exactly the path pairs agree, and all other pairs disagree. Conversely, when the agreeing pairs form P, A and the pullback of B agree on P and disagree on its complement, which is exactly the definition of K_P(A). The equivalent graph conditions hold because a digraph with both degrees at most one is a disjoint union of directed paths and directed cycles; n-1 arcs plus no cycles gives exactly one component. Weak connectivity gives the same conclusion. ∎

**Hostile.** Arc count and indegree/outdegree bounds without connectivity/acyclicity are NOT sufficient. On four vertices take A=C3->K1 and B=K1->C3, keeping the same labels on C3. Their agreement digraph is C3 plus one isolated vertex: three arcs and both degrees at most one, but not a Hamilton path. This example also makes clear why a cyclic trace or an undirected support cannot silently replace directed path data.

For merged classes U and V, apply the criterion to either of their one or two ordinary classes. Thus it is an iff decision procedure for every merged adjacency. This is a verifiable realization certificate for named endpoints, not a theorem characterizing arbitrary abstract graphs as entire metagraphs.

### Exact multiplicities, colors, and loops

Let a_A=|Aut(A)|. Let N_AB count the bijections in Theorem 1. Let N_AB^blue count those for which the unique path P=(p_1,...,p_n) has reversal r_P(p_i)=p_(n+1-i) an anti-automorphism of A. Set

    I_AB^blue = N_AB^blue/(a_A a_B),
    I_AB^black = (N_AB-N_AB^blue)/(a_A a_B).

**Theorem 2.** These are nonnegative integers, symmetric in A,B, and give the exact directed incidence matrices of the unmerged colored line multigraph. An off-diagonal entry is its line multiplicity; a diagonal entry is twice its loop multiplicity.

*Proof.* Each actual path P with K_P(A) isomorphic to B gives exactly a_B identifying bijections f. Automorphisms of A act freely on ordered Hamilton paths, so division by a_A gives the distinct pinned-path presentations of A with partner class B. The blue predicate is invariant under this action. In path coordinates, reversal anti-automorphism is exactly the mirror-fixed condition. Complement preserves that predicate. Complement pairs the presentations bijectively and is fixed-point-free for n>=3, proving symmetry and the diagonal factor two. ∎

For merged nodes U,V define J_UV^c=sum I_AB^c over ordinary A in U and B in V. The full merged blue multiplicity and the **fundamental paired-black** multiplicity are

    U != V:   blue = J_UV^blue,       fundamental black = J_UV^black/2;
    U  = V:   blue loops = J_UU^blue/2, fundamental black loops = J_UU^black/4.

The raw merged black graph has twice the fundamental multiplicity. Loops contribute twice to degree in both conventions. Proof of the additional division by two is the inherited mirror-twin theorem; it is not an extra unidentified normalization. These formulas determine the complete colored multigraph from tournament representatives, including all parallel edges and loops. They do not infer graph shape from degree counts.

## 2. Exact port-aware substitution transfer

Let T=A[B_1,...,B_m] be a tournament substitution with named disjoint module vertex sets V_i. All cross arcs follow A. Consider the following finite data:

1. A word w=(i_1,...,i_r) with adjacent letters different and each transition i_j->i_(j+1) an arc of A.
2. For each i, an ordered list C_i of nonempty vertex-disjoint directed paths in B_i, whose union is V_i, with one path for each occurrence of i in w.
3. At each occurrence use the next path in C_i; retain its actual entry and exit vertices.

Concatenating these paths gives a Hamilton path P of T. Conversely every Hamilton path of T has a unique such encoding: take its maximal constant-module runs. No state is counted twice.

Write E_i for the set of internal run arcs in C_i and F for the set of actual exit-to-entry arcs at consecutive runs.

**Theorem 3 (port-faithful complement).** K_P(T) is obtained by the following explicit construction:

- Inside V_i start with B_i^op and restore precisely the arcs of E_i to their original direction.
- Between modules start with A^op[B_1^op,...,B_m^op]'s uniform cross directions and restore precisely the individual port arcs in F.

In particular the run data give all partner tournaments and all their incidences, with multiplicities, upon enumeration and ordinary isomorphism reduction. The construction intertwines complement exactly on this marked state space.

*Proof.* The n-1 path arcs partition disjointly as (union E_i) disjoint union F. Every such arc is restored after global converse and no other arc is restored. This is K_P by definition. Uniqueness of maximal runs proves exhaustive, duplicate-free path enumeration. ∎

At class level one may count partners in this complete run enumeration and divide by |Aut(T)|, because the source action on full paths is free. This remains valid if automorphisms move a chosen module decomposition: division is performed on the full path count, not on an assumed faithful action on arbitrarily named blocks. To mark blue states, test the actual full reversal, not a guessed block reflection.

### Block-contiguous special case and exact missing coordinate

When each module occurs once, each C_i is one Hamilton path. Internal partners are K_(C_i)(B_i). Cross directions are the converse quotient directions except for the ONE exit-to-entry arc at each consecutive pair. Equivalently start with the uniformly substituted complemented quotient: on each consecutive block pair one must reverse every cross arc EXCEPT the designated boundary port. Exactly |B_i||B_j|-1 arc changes are required at that boundary. Thus the uncorrected operations commute literally precisely when every adjacent size product is one. For m>=2 this forces every module to be a singleton. This is a statement about the named block decomposition, not an assertion that accidental unlabelled isomorphisms are impossible in every instance.

The arbitrary-run theorem is stronger than this special case: a module may be visited repeatedly, so it is a path cover and its individual attachment vertices, not just one internal Hamilton path, that must be retained. Numerical cover vectors alone preserve the number of paths but do not supply the orientations of the cross-port corrections.

## 3. Genuine S2 and S3 carriers; the synchronization obstruction

On pinned-path presentations the two operations kappa (nonpath complement) and sigma (converse followed by reversed path labels) are commuting involutions. They generate a C2 x C2 action. The endpoint pair is fixed under sigma and toggled by kappa, so neither kappa nor sigma*kappa fixes a presentation. Therefore each blue edge is a two-state orbit, while each fundamental black edge is a four-state orbit. This is the exact S2-level origin of the factors in Theorem 2. It is not S3.

A genuine S3 action appears on an oriented triangle WITH a path frame. Fix C3 on labels 0,1,2, and form

    F = {(epsilon,p): epsilon in {+1,-1}, p a Hamilton path in C3^epsilon}.

There are six states. For g in S3,

    g.(epsilon,p) = (sgn(g)*epsilon, g(p)).

This is a free transitive action: a full ordered path fixes all three labels. It is regular S3, not regular C6. The three rotations preserve orientation, while each transposition reverses it. On this same carrier signed path reversal (epsilon,p)->(-epsilon,reverse(p)) commutes with left relabeling: it is a right action, not another independent identification of all six states.

**Theorem 4 (synchronized local signs).** In T=A[C3], for |A|>=2, every signed automorphism preserving the bottom triangle blocks has the form

    (g; pi_1,...,pi_m),
    g in Aut^epsilon(A),  pi_i in S3,  sgn(pi_i)=epsilon for EVERY i.

Conversely all these maps are signed automorphisms. For this substitution the bottom triangle blocks are intrinsic, by Session24's cyclic-three-module recovery, so these are all signed automorphisms. In particular Aut(T)=C3^m semidirect Aut(A); if A is self-converse, exactly another equally large coset has negative sign.

*Proof.* Cross-block arcs force g to preserve or reverse all quotient arcs with the global sign epsilon. Internal arcs in every triangle force the same sign for pi_i. These conditions also suffice pairwise, internally and externally. Intrinsic bottom-block recovery is proved in Session24 CONNECTIONS section 2: cyclic triples spanning two blocks are impossible; one spanning three blocks is distinguished by an unused vertex of one of its blocks. ∎

Thus independent arbitrary S3 actions on the modules are invalid: a reflection in just one triangle reverses internal arcs while preserving cross arcs, hence is neither an automorphism nor an anti-automorphism. For the two-vertex transitive base there are 18 signed block permutations, versus 72 total block permutations (S3 x S3) semidirect S2. The latter carrier is four times too large.

This is a genuine shared-character/fiber-product phenomenon with the earlier THM-3173-six-state-free-factor-actions-and-pointed-frame-cube.md: both require equality of actual binary characters rather than treating two six-state sets as interchangeable. The actions and carriers differ, and no identification with the Farey/C6 carrier is claimed.

**Important boundary.** Nonpath complement of a directed triangle in one of its Hamilton-path frames is transitive, not cyclic. Hence the six-state triangle carrier is not closed under kappa. An S3 triangle operation is not automatically a metagraph symmetry. It acts on signed triangle frames and compatible substituted states; transporting it through complement requires the port-aware target construction above.

## 4. Verification and remaining frontier

Run `python3 check.py` from this directory. The code has no repository imports or third-party dependencies.

Two implementations are deliberately different:
- enumerate class-pair bijections and inspect their agreement digraphs;
- enumerate every pinned-path tiling, complement its free bits, and canonicalize its endpoints.

They compare every ordered pair and both colors, not just aggregate counts. The output states the finite orders, numbers of bijections, exact lines, loops, and fundamental black bridges. Bridge detection retains parallel edges and ignores loops correctly. It does not treat a degree-one calculation as a bridge classification.

The port test independently compares direct nonpath complement with restoration from maximal runs for all directed Hamilton paths of every selected substituted tournament: quotient orders 2 and 3, every labelled quotient orientation, module sizes 1,2,3, fixed transitive singleton/pair and cyclic triangle modules, total order at most 7. This is a deliberately bounded test of an all-order elementary identity, not an exhaustive substituted-tournament census.

The agreement graph hostile, regular S3 carrier, and two-triangle sign-synchronization test are independent additional controls. No tests claim that all odd H other than 7 and 21 occur, that a proposed invariant tuple is sufficient, or that all non-leaf bridges have been ruled out.

The remaining structural problem is to compress this exact certificate/transfer while preserving enough port information to prove nontrivial all-order connectivity, bridge, or adjacency families. An exact state representation is necessary infrastructure, not by itself the desired classification.

## 5. A structural consequence: labelled odd girth and quotient holonomy

This section is NOT a claim about the unmarked merged metagraph. Let L_n have **every labelled tournament on one fixed vertex set** as a vertex, with an edge T--K_P(T) for every directed Hamilton path P. For n>=3 this graph has no loops. Two distinct paths cannot give the same labelled partner, because their edge sets would agree, and a directed Hamilton path is determined by its edge set.

Encode a tournament by an orientation bit vector t on the N=binom(n,2) unordered pairs, and write p for the characteristic vector of its path edges. Each move is

    t -> t + 1 + p       over F2.

For a sequence of k valid moves with path vectors p_1,...,p_k,

    t_k + t_0 = (k mod 2)*1 + p_1+...+p_k.                 (5.1)

This elementary conservation law retains information discarded by the isomorphism quotient.

**Theorem 5 (labelled odd-cycle obstruction and sharp family).**

1. Every odd closed walk of length k in L_n has k>=ceil(n/2).
2. L_n is bipartite whenever n is 0 or 3 modulo4.
3. For n=4r+2>=6 its odd girth is exactly n/2=2r+1. In particular its shortest odd cycles have unbounded length with n.
4. L_5 has odd girth3. No general exact odd-girth formula for n=1 modulo4 is asserted here.

*Proof of 1.* For an odd closed walk, (5.1) gives p_1+...+p_k=1. Every pair therefore occurs an odd, positive number of times among the paths. Since each path uses n-1 pairs, k(n-1)>=N, giving k>=n/2.

*Proof of 2.* Each move toggles exactly m=N-(n-1)=(n-1)(n-2)/2 bits. This integer is odd precisely for n=0,3 modulo4. Parity of the number of 1-bits is then a bipartition.

*Proof of 3, explicit construction.* Put n=2k with k odd and k>=3; use labels in Z/(2k). For i=0,...,k-1 take the undirected Hamilton path

    P_i=(i, i-1, i+1, i-2, i+2, ..., i-(k-1), i+(k-1), i-k),

with entries reduced modulo2k. Its edges have endpoint sum either 2i-1 or 2i modulo2k. The odd-sum pairs form a perfect matching of k edges; the even-sum pairs form a matching of k-1 edges (the two formal self-pairs are omitted). These are exactly the edges of P_i. As i varies the sums cover every residue once, so the paths partition all edges of K_(2k).

Orient each pair initially so that its occurrence in path P_i is forward in the displayed order after i prior moves: use that displayed direction for even i and the reverse for odd i. A pair in P_i has been toggled in each of the preceding i moves, so P_i is indeed a directed path at step i. After all k moves each pair was retained once and toggled k-1 times, an even number; the tournament returns to its start. No intermediate repetition is possible: a proper segment of r moves toggles the nonempty union of its disjoint path sets if r is even, or their nonempty complement if r is odd. Thus this is a simple k-cycle, matching the lower bound.

*Proof of 4.* On labels {0,1,2,3,4}, take three path words

    (2,0,1,3,4), (3,0,1,4,2), (4,0,1,2,3).

Their common edge is 0->1; every other unordered pair occurs exactly once. Give that common edge the displayed direction, and orient every other edge so its displayed direction is forward at its step, as in the previous construction. The common edge is never toggled; every other edge is toggled twice. This gives a labelled triangle. `holonomy.out` includes its explicit bit-vector states. ∎

### Exact holonomy correction for quotient walks

A walk closed only after relabeling need not obey the labelled bound. If the end is g.T_0, define Delta_g=t(g.T_0)+t(T_0), the unordered-pair disagreement vector. Equation (5.1) says exactly

    p_1+...+p_k = (k mod2)*1 + Delta_g.

For odd k this forces

    |Delta_g| >= binom(n,2)-k(n-1).

If the walk closes after both relabeling and converse, replace k mod2 by (k+1) mod2 on the right. More uniformly, for signed holonomy g with converse bit delta, the all-pairs term is (k+delta) mod2. Hence short ordinary-class odd cycles at large n necessarily carry substantial nonidentity relabeling disagreement; merged cycles can additionally use the converse sign. They do not contradict the labelled odd-girth theorem.

This is the missing coordinate in attempts to deduce a quotient graph's global cycle shape from a path-state operation: one must track the vertex identification accumulated around a walk, not only its endpoint class, degree, or S2/S3 label. The statements establish neither connectedness nor non-leaf-bridge exclusion for the merged metagraph.

Reproduction: `python3 holonomy.py`. The implementation constructs and checks actual directed cycles, not only path-decomposition counts, at n=6,10,14,18,22; checks the explicit n5 triangle; and exhausts labelled path moves at n3,4 for the parity control. The proofs cover the full stated families independently of these finite tests.

### Equality anatomy in the unresolved odd-order residue

There is also an exact restriction on a hypothetical lower-bound cycle at n=4r+1. Put k=(n+1)/2=2r+1. If an odd closed walk has length k, each edge is used an odd positive number of times. At vertex v let e_v be the number of these paths for which v is an endpoint. Its total incident path-edge count is 2k-e_v=n+1-e_v, whereas covering every edge at least once requires n-1. Thus e_v<=2 and the excess incident count is 2-e_v. Since every excess edge multiplicity is even, e_v is even. Consequently exactly k vertices have e_v=2, and the other k-1 vertices have e_v=0. At every latter vertex exactly one incident edge is used three times; no edge is used five times. The triply used edges therefore form a matching of size (k-1)/2=r, and all other edges occur once.

So attainment of the lower bound at n=4r+1 necessarily gives a Hamilton-path decomposition of the multigraph K_n+2M for a matching M of size r, together with orientation compatibility on every tripled edge. This condition is also sufficient when that compatibility is supplied: for an edge occurring at step i (zero-based), with j earlier occurrences, its initial direction must be its displayed path direction reversed i-j times. These requirements must agree at its three occurrences. If they do, they define an initial tournament, every proposed step is valid, and the odd cover returns to the start. The lower bound excludes any shorter odd subcycle, so the resulting closed walk of length k is a simple cycle.

This reduces the remaining n=1 modulo4 extremal question to an explicit path-decomposition and sign-compatibility problem. The n5 construction has M={01}. It is not asserted here that compatible decompositions exist at every larger order.

**FINITE-EXACT additional boundary closure: n9.** The lower bound is attained at n9 as well, so L9 has odd girth5. An exact directed cycle is certified by the following successive path words:

    (5,0,8,6,1,2,3,7,4)
    (4,2,3,6,0,1,5,7,8)
    (5,4,8,3,0,1,7,2,6)
    (6,5,2,8,1,3,4,0,7)
    (7,6,4,1,0,2,3,5,8)

The repeated matching is {01,23}. With bits ordered by lexicographic unordered pairs and a 1 denoting the lower label pointing to the higher, the tournament states are

    5927068074 -> 45360910277 -> 59874146331 -> 9281284065
    -> 58367449686 -> 5927068074.

`holonomy.py` checks every directed arc and every complement transition using only Python integers. The optional `probe_n9.py` documents discovery by SciPy/HiGHS MILP followed by exhaustive ordering/orientation search; its output is not trusted without the independent integer verifier. This closes one additional finite boundary, not n=1 modulo4 in general.
