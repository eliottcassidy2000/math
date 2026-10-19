# The forest primitive floor and the two parity functions of the eight-vertex double-path

2026-10-19. **PROVED in the stated fixed-order and source-stack scopes; exact finite witnesses supplied; independently reviewed.** No external-priority claim. This is about the eight-vertex tree D itself. It is separate from the already proved nine-vertex U_c/U_a leaf selectors.

## 1. Inherited inputs and the actual new conclusion

Use the reviewed [nine-tree packet](../twin_signature_parity_20261019/NINE_VERTEX_GOOD_TREES.md), the [twin-signature packet](../twin_signature_parity_20261019/TWIN_SIGNATURE_PARITY.md), and THM-4533, `01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md`. Their relevant proved inputs are:

- Every orientation of every tree through six vertices has constant spanning embedding parity on tournaments of the same order.
- Every orientation of P4 has odd spanning embedding parity.
- A symmetric seven-tree is good, since all its proper subtrees have at most six vertices and hereditary symmetry suffices.
- Deletion–reversal gives phi_S+phi_(reverse_e S)=phi_(S-e), retaining the entire vertex set.
- Disjoint unions of parity-rigid components have parity equal to the product of their constants times the multinomial shuffle coefficient.

The earlier nine-tree packet proves that D plus one isolated vertex is odd in every orientation on every nine-tournament. It does not give the following classification of D's own 128 labelled orientations on eight-tournaments.

Let D have vertices 0,...,7 and ordered edge list

    01, 02, 03, 14, 15, 26, 47.

An orientation mask has seven bits: bit one reverses the corresponding listed edge, and bit zero retains the listed direction. The critical bits are 2 and 4, namely the leaf edges 03 and 15. Let A be D with leaf5 deleted, so A is the asymmetric seven-vertex spider with leg lengths 1,2,3. Let F(T) be the embedding parity of A plus one isolated vertex, in an eight-tournament T. As proved below, F does not depend on the chosen orientation of A.

**Theorem.** For every orientation mask m and every eight-tournament T,

    phi_(D,m)(T)=1+epsilon(m) F(T) in F_2,
    epsilon(m)=bit_2(m)+bit_4(m) mod2.               (1)

The function F is nonconstant. Consequently:

- Exactly 64 of D's 128 labelled orientations have odd embedding count in every eight-tournament.
- The other 64 have one and the same nonconstant parity function 1+F.
- There are exactly 32 Rédei oriented isomorphism classes and 32 nonconstant oriented isomorphism classes of D. Every orientation has trivial oriented automorphism group.

No classification of all eight-trees or all good trees is claimed. Section4 additionally proves a specific all-order primitive-degree floor for forests.

## 2. Proof of the universal selector

### Orientation independence of F

Delete any edge of A. Its two tree components have orders a and 7-a, with 1<=a<=3. Both have at most six vertices, hence their oriented embedding parities are constant. Retaining the isolated eighth vertex, the disconnected forest has shuffle coefficient

    8!/(a!(7-a)!)=8 binom(7,a),

which is even. Its parity is therefore zero for every orientation. Deletion–reversal shows that flipping any single edge of A does not change F. All orientations are connected by such flips, proving the claimed independence.

### Every first derivative of D

- Deleting 03 or 15 leaves a copy of A plus one isolate, hence gives F.
- Deleting 01 leaves two P4s. Their shuffle coefficient binom(8,4)=70 is even, so this derivative is zero.
- Deleting 02 or 14 gives components of orders two and six, both good, with even shuffle coefficient binom(8,2)=28.
- Deleting 26 leaves a seven-tree in which vertices2 and3 are twin leaves at0; deleting 47 gives the reflected case with twin leaves at1. Each seven-tree is symmetric and therefore good. With the retained isolate, the shuffle coefficient is eight, so both derivatives are zero.

All these statements are pointwise on every eight-tournament and hold for every orientation of the remaining edges. Thus only bits2 and4 affect the parity function, and either flip adds exactly F.

### The odd anchor

At mask0 the two rooted P4 halves, 3-0-2-6 and 5-1-4-7, are oriented compatibly, and 01 is directed from the first root to the second. Partition any eight-host into an unordered pair of four-sets. Each P4 has odd embedding count on either four-set. Swapping the two labelled halves pairs candidate embeddings into pairs for which exactly one direction of the root-to-root edge is present. Therefore the parity contribution of each unordered partition is one. There are binom(8,4)/2=35 partitions, an odd number. Hence mask0 has constant odd parity.

Starting at mask0 and applying the previously determined first derivatives proves (1).

### Oriented isomorphism classes

The two rooted P4 halves are rigid and isomorphic. Each has branches of lengths one and two from its root, so no nonidentity rooted automorphism exists. The undirected automorphism group of D is therefore exactly C2, generated by the swap

    (0 1)(2 4)(3 5)(6 7).

That swap reverses the central edge, so no orientation is fixed by it. All 128 labelled masks split into pairs. The critical-direction agreement condition is invariant under swapping halves, giving 32+32 oriented isomorphism classes, all with trivial oriented automorphism group.

## 3. An explicit nonconstant defect

Use the eight-tournament whose out-neighbor bit masks are

    [210,188,169,81,164,201,22,72].

This is the D8 witness already frozen in `twin_signature_parity_20261019/minimality_witnesses.json`; the current calculation retains its lineage and adds the common defect count. Let F use D with edge15 deleted and vertex5 retained as the isolate, with every remaining listed edge directed forward. Exact counts are:

    transitive host:          emb(F)=480,  emb(D,mask5)=165;
    displayed hostile host:   emb(F)=617,  emb(D,mask5)=302.

Thus F is zero on the transitive host and one on the displayed host. This proves its nonconstancy and proves that the epsilon-one branch in (1) is genuinely nonconstant, rather than merely uncertified.

The supplied checker computes all 128 D orientation counts simultaneously by literal permutation enumeration, verifies every first derivative and selector value on both hosts, and independently checks selected counts by an injective backtracking algorithm. It also enumerates the complete undirected automorphism group and all oriented orbits, and replays the sharp seven-core witness in Section4. It preserves the isolate in F. Its 2,428 checks pass identically normally and with Python optimization. The two-host calculation is a witness and implementation control; the universal host quantifier in (1) comes from Section2.

## 4. The integral-basis consumers

### 4.1 A sharp all-order forest primitive floor

Work modulo two. For every oriented forest S on n vertices, its expansion in the source-free induced-subtournament basis has zero coefficients at every core order

    1<=r<=min(6,n).                                 (2a)

Equivalently, all nonconstant primitive data of an arbitrary forest embedding parity begin at core order seven.

**Proof.** By iterating the source-evaluation identity, h^(n-r) phi_S is a sum, with its actual multiplicities, of spanning embedding-parity functions of oriented r-vertex forests. For r<=6, every tree component is good. The disjoint-component shuffle formula therefore makes every such forest parity constant. Hence h^(n-r) phi_S is constant. The integral Taylor projector P, reduced modulo two, kills every constant in positive degree: on C_r it gives sum_j(-1)^j binom(r,j)=(1-1)^r=0. Thus P h^(n-r)phi_S=0, exactly the primitive coefficient at degree r. This proves (2a) for every n; no bound on the original forest size is used.

A direct host-level consequence, retaining the fixed spanning order, is

    phi_S(R with n-r successive universal sources)
       = phi_S(transitive n-tournament),     0<=r<=min(6,n),   (2b)

for every r-vertex tournament R. Indeed h^(n-r)phi_S is constant on all r-hosts, and its value at a transitive r-host is the right side. Thus replacing the bottom at most six vertices of a transitive source stack by an arbitrary tournament is invisible to every spanning oriented-forest parity of that order. This is a source-stack statement; it is not invariance under arbitrary subtournament substitutions or graph minors.

The threshold seven is sharp. The inherited seven-tree witness has

    edges 01,02,03,14,25,46; orientation mask52,
    source-free hostile out-neighbor masks [118,44,96,85,102,8,34].

Its exact embedding counts are 66 on the transitive seven-tournament and 81 on that hostile host. The hostile has no universal source, so its core order is exactly seven. These counts are independently replayed by both permutation enumeration and injective backtracking. Their provenance is the n=7 record of the same inherited `minimality_witnesses.json` file. They show that (2a)--(2b) cannot be extended through core order seven.

The forest hypothesis is substantive. A directed triangle already has three spanning embeddings into the cyclic three-host and none into the transitive three-host. The claim is not an all-oriented-graph primitive floor.

### 4.2 The double-path has one degree-seven primitive defect

Use the source-free induced-subtournament basis and source operator h from the reviewed [Hasse deletion packet](../tournament_hasse_frobenius_20261019/PADDING_FROBENIUS_PROOF.md). Work modulo two. Put

    c=phi_A(transitive seven-tournament),
    kappa_A=phi_A+c·1_7.

The source-evaluation identity h phi_A sums the spanning parities of the oriented vertex deletions A-u over source-eligible vertices u. Every such deletion is a forest on six vertices; each component is good, so its parity is constant. Hence h phi_A is constant. Evaluating it at the transitive six-tournament shows that constant is c. Therefore

    h kappa_A=0,             kappa_A belongs to K_7.

Ordinary padding d=D_1 satisfies d1_7=8·1_8=0 modulo two, so

    F=d phi_A=d kappa_A,
    phi_(D,m)=1_8+epsilon(m) D_1 kappa_A.            (2)

This localizes every nonconstant term to source-free core order seven. All nonempty primitive coordinates of orders below seven and the same-order primitive coordinate of order eight vanish. The primitive is nonzero because F is nonconstant; D_1 is injective on the primitive K_7 summand in the integral reconstruction basis.

This is an actual consumer of the new basis: D's 128 orientation functions take only the two affine values 1 and 1+D_1 kappa_A, and its nonconstant information is erased by the next padding since d squared is zero. It explains the inherited all-nine-host padding phenomenon without claiming a general classification.

## 5. A credited Hamiltonian-path consumer of the same basis

This paragraph identifies an inherited theorem in the new coordinates; it is not a new forbidden-value theorem. Let H(T) be the number of directed Hamiltonian paths. Adjoining a universal source forces that vertex to be first, so hH=H across degrees. Therefore the source-free coefficient at a core A is

    a_A(H)=sum_(S subset V(A)) (-1)^|S| H(A-S).

The odd-cycle-collection formula in [THM-002](../../../01-canon/theorems/THM-002-ocf.md), credited there to the Rédei–Berge work of Grinberg–Stanley, writes H(B) as the sum of 2^(number of cycles) over vertex-disjoint collections of directed odd cycles in B, including the empty collection. Substituting that formula above, inclusion-exclusion cancels every collection with any unused vertex. Hence

    a_A(H)=sum_(spanning odd-cycle covers of A) 2^(number of cycles).

This is a nonnegative integral coefficient and is even for nonempty A. It identifies the inherited cycle-cover expansion with the primitive reconstruction basis. It supplies no new exclusion of Hamiltonian-path counts7 or21 and no new extremal theorem.

## 6. Reproduction

From this packet directory:

    python3 -B verify_double_path_selector.py
    python3 -B -O verify_double_path_selector.py

Both reproduce the supplied selector receipts.
