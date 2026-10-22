# Labelled attachment responses: repaired one-probe formula and finite equivalence

Status: PROVED formula and finite-equivalence reduction; FINITE-EXACT two-probe separation through core order six. Universal two-probe reconstruction is OPEN here. No novelty claim for weighted-automaton equivalence. No repository files or Git state were changed by this worker.

## Objects, inheritance, and loss ledger

A labelled tournament Q has vertex set [k]. Attach an ordered transitive base x_1,...,x_m with x_i→x_j iff i<j. A mask a_i has bit u equal to 1 iff core vertex u→x_i. Let F_Q(a_1...a_m) be the number of directed Hamiltonian vertex words. The empty response is H(Q). Labels and masks remain fixed when comparing two cores; this is not an unlabelled-isomorphism question.

Closest inherited mechanism: Session26's exact Hamiltonian slot automaton, `/workspace/shared/session26-counter/proof.md` and `counter.py`. Its final-state dependence on Q, not just its state bound, is essential. The ordered path-cover mechanism traces to THM-3181, `THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`. The hostile is opposite labelled C3 orientations: identical zero/one-probe responses but different two-probe responses. The corrected near miss is the insertion-only one-probe formula below. The least-used sidecar is the ordered endpoint pair of a spanning two-path cover; its unordered sum is insufficient in general.

Concept board: labelled ports; ordered path covers; symmetric versus antisymmetric endpoint data; common transfer realization; reachable linear quotient. Scalar averaging and fixed-depth measurements destroy coordinates, whereas a common finite transfer realization makes equality decidable without promising reconstruction.

## 1. Correct one-base polynomial

For k≥1 let S_u and E_u count Hamiltonian paths of Q starting and ending at u. For distinct u,v let C_uv count ordered pairs (P_1,P_2) of nonempty vertex-disjoint directed paths covering all core vertices, with P_1 ending at u and P_2 starting at v. Set C_uu=0. Single-vertex paths are allowed.

Then for a Boolean mask b,

F_Q(b) = H(Q) + Σ_u [E_u−S_u+Σ_v C_uv] b_u − Σ_{u<v}(C_uv+C_vu)b_u b_v.

Proof: split every extended Hamiltonian path according to the position of the new vertex x. If x is first the contribution is Σ_u S_u(1−b_u). If x is last the contribution is Σ_u E_u b_u. If x is internal, deletion of x produces exactly an ordered spanning two-path cover, and the required two new arcs contribute C_uv b_u(1−b_v). These cases are disjoint and exhaustive. Expand and collect the Boolean multilinear coefficients. This also proves degree at most two and gives a unique coefficient representation on {0,1}^k.

This formula does not assert that the separate S, E, or ordered C entries are recoverable from one-base responses. Precisely the displayed linear combinations and unordered C sums are recovered. This distinction is important for census refinement: grouping by a finer unrelated statistic can lose genuine one-response collisions.

### Refuted insertion-only shortcut

The tentative formula H(Q)+Σ_v(H(Q)−S_v)b_v−Σ_{u<v}E^{HP}_{uv}b_u b_v is false. Here E^{HP}_{uv} means the number of core Hamiltonian paths using the actual oriented edge between u and v. For Q:0→1 and b=(0,1), the extension is C3 and F_Q(b)=3, while that shortcut gives 2. The first failed implication is that deletion of an internal inserted vertex leaves a Hamiltonian path of Q. It may instead leave two paths whose connecting edge is oriented backwards. The correct endpoint-cover data include such bypasses.

For this hostile H=1, S=(1,0), E=(0,1), C_01=C_10=1, hence the correct polynomial is 1+2b_1−2b_0b_1.

## 2. Finite-exact refinement census

`census.cpp` enumerates every labelled tournament by increasing unordered edge pairs. A set mask bit gives the smaller→larger orientation. It computes all Hamiltonian paths and all ordered spanning two-path covers by core permutations: a valid cover corresponds to one cut with all other adjacent edges present. Equivalently a permutation has zero bad adjacent edges and permits every cut, or exactly one bad adjacent edge and must be cut there.

It groups cores by the correct polynomial coefficients, then only refines nonsingleton groups by exact two-base responses. The latter are computed independently by full-vertex subset endpoint DP, not the polynomial or slot automaton. Refinement continues until all groups are singletons or every ordered mask pair is exhausted. Singleton groups never need more probes. The printed number of remaining collision groups can increase during refinement because one group may split into multiple nonsingleton groups; that is not a failure.

Results from `census.out`:

- k=3: 8 cores, 7 one-response classes, 1 collision group of size 2; zero two-response collisions.
- k=4: 64 cores, 56 one-response classes, 8 collision groups of size 2; zero two-response collisions.
- k=5: 1024 cores, 932 one-response classes, 92 collision groups of size 2; zero two-response collisions.
- k=6: 32768 cores, 31244 one-response classes, 1484 collision groups of maximum size 4; zero two-response collisions.

Core orders 0,1 have only one tournament, and k=2 has its two orientations already separated at one base vertex. Thus zero/one/two-base data distinguish all labelled tournaments through order six. No all-order implication follows from this census.

The earlier exploratory groups 911 for k=5 and 30992 for k=6 used Hamiltonian edge usage instead of spanning two-path covers. They were not one-response census counts and are superseded, not evidence for the final statement.

## 3. Finite witness theorem for unrestricted response equality

For any two tournaments Q,R on the same labelled k-vertex set, k≥1, set D=2 k! 3^k. Then the following are equivalent:

1. F_Q(w)=F_R(w) for every finite attachment word w.
2. F_Q(w)=F_R(w) for every word with |w|≤D−1.
3. The common reachable row space annihilates f_Q−f_R in the inherited slot realization.

This is an equality-decision theorem, not a proof that full response determines Q.

Proof: Session26 gives a D-state representation F_Q(w)=α M_w f_Q. For every permutation of the core, the slot transitions depend only on the permutation, attachment mask and empty/occupied endpoint states. Q enters only through the final test of empty interior slots. Therefore α and all M_a can be shared between Q and R. Their difference is represented in dimension D, not 2D:

F_Q(w)−F_R(w)=α M_w δ, δ=f_Q−f_R.

Let V_m=span_Q{α M_w:|w|≤m}. V_0 has dimension one, since α is nonzero. Also V_{m+1}=V_m+Σ_a V_m M_a. Once V_m=V_{m+1}, this space is invariant under every M_a and equals the full reachable space. Until stabilization, dimension strictly increases. Since ambient dimension is D, stabilization occurs by m=D−1. Consequently all differences vanish iff δ annihilates V_{D−1}. If the series differ, some word of length at most D−1 detects the difference. More sharply, D can be replaced by the reachable dimension r, giving r−1. k=0 has only the empty core and constant response one.

The number D is an inherited upper bound, not a minimal realization or practical complexity claim. Alphabet size is 2^k, and D grows factorially/exponentially with k. Explicitly testing every word through D−1 is unnecessary: construct a basis by linear reachability.

## 4. Executable exact equality test

`equivalence.py` is self-contained Python standard library. It builds the common slot matrices with repeated destinations retained as integer transition multiplicity, constructs δ, and explores reachable rows. A candidate row is kept only if independent of previously stored rows, tested by exact Fraction elimination. Only independent rows are expanded. At most D rows are stored for expansion and at most 2^k D candidates are tested, apart from initialization. If a row has nonzero dot product with δ, its retained actual input word is returned as a witness. Otherwise completion gives an invariant spanning space annihilating δ and certifies equality.

Rational elimination does not alter the actual row-word pairs used for expansion: reduced vectors are kept separately for membership tests. There are no floating-point or modular-zero shortcuts. Arithmetic and memory costs include rational coefficient growth; no unit-cost estimate is presented as bit complexity. Intended use is small k.

Positive controls exercise the equality branch rather than returning immediately when Q==R. For k=2 the common reachable dimension is 28 and the longest retained basis word has length 4, below the crude D−1=35. Opposite C3 masks 2 and 5 yield word (1,2) and response difference 2, below D−1=323. This witness can be checked independently by full-tournament endpoint DP.

## 5. Classical attribution and scope

The linear-dependence-pruned word tree is classical weighted-automaton equivalence, not a new general algorithm:

- Wen-Guey Tzeng, “A Polynomial-Time Algorithm for the Equivalence of Probabilistic Automata,” SIAM Journal on Computing 21(2), 1992, 216–227. Official record: https://epubs.siam.org/doi/10.1137/0221017 . The official abstract states a polynomial-time equality test.
- Stefan Kiefer, Andrzej S. Murawski, Joël Ouaknine, Björn Wachter, James Worrell, “Language Equivalence for Probabilistic Automata,” author-hosted paper: https://people.mpi-sws.org/~joel/publications/langeqprobaut11.pdf . Section 2 explicitly treats real-weighted automata and describes pruning by linearly dependent state distributions; Section 2.1 gives forward/backward space algorithms. The paper also attributes the broader mechanism to Schützenberger.

The finite witness D−1 bound above is proved directly; the cited paper prints a slightly looser n bound for two arbitrary automata. Our useful specialization is the shared slot realization: only the final vector changes, so there is no need to double dimension. No priority claim is made for that elementary specialization or for the one-probe cover identity.

## Reproduction and boundaries

From this directory:

    g++ -O3 -std=c++17 census.cpp -o census
    for k in 3 4 5 6; do ./census "$k"; done
    python3 audit.py
    python3 -O audit.py
    python3 equivalence.py 3 2 5
    python3 equivalence.py 2 0 0

`audit.py` independently verifies the one-base polynomial against endpoint DP for all labelled cores through k=5 and all one-base masks, checks the common automaton against endpoint DP through k=3 and all words through length 3, and exercises positive and negative exact-equivalence controls. Explicit RuntimeErrors, not assertions, preserve checks under optimization.

Open here: universal two-base reconstruction; full-response injectivity on arbitrary labelled tournaments; structural characterization of the one-response collision classes; minimal observable state dimension; and improved uniform witness length. The finite-equivalence theorem answers none of those injectivity questions by itself.

## 6. Converse boundary current and invisible circulation

Status: PROVED. The terminology “converse current” below is introduced only as a descriptive label here, not attributed as established graph terminology.

Let Q^op reverse all core arcs while retaining every label. Reversing a Hamiltonian word exchanges S and E. Reversing both paths of an ordered two-path cover and then exchanging their order gives C^op_uv=C_vu. Thus every quadratic coefficient of the one-probe response is automatically converse-invariant. Define the skew endpoint-cover current

A_uv=C_uv−C_vu,

its divergence div(A)_u=Σ_v A_uv, and the converse boundary vector

g_u=2(E_u−S_u)+div(A)_u.

Then exactly

F_Q(b)−F_Qop(b)=Σ_u g_u b_u.

In particular Q and Q^op have identical responses for every one-base attachment iff g=0. This is a full criterion, not just a sufficient symmetry condition. Since Σ_uE_u=Σ_uS_u=H and a skew matrix has total row sum zero, Σ_ug_u=0.

The lost coordinate has a precise linear description. With endpoint boundary E−S fixed, this comparison sees only divergence of A, not its circulation component. On the complete labelled graph, the formal real vector space of skew currents has dimension k(k−1)/2; divergence has image the zero-sum vectors, dimension k−1 (k≥1), and kernel dimension (k−1)(k−2)/2. The kernel consists of cycle currents and is spanned, for example, by oriented triangles based at a fixed vertex. This is a statement about the linear observation map on formal currents. It DOES NOT say that arbitrary currents, endpoint counts, or perturbations in that kernel are realized by tournament path-cover counts. Also g=0 may involve cancellation between nonzero endpoint boundary and nonzero divergence; it need not imply either separately vanishes.

### Vertex-transitive cores are always one-probe converse-blind

If Q is vertex-transitive, then S_u and E_u are each constant in u, hence both equal H/k. The row sums of C are constant by transitivity, as are its column sums; their total sums agree, so the constants agree. Consequently div(A)=0 and g=0. Arc-transitivity is NOT required for this converse-invisibility conclusion.

If Q is additionally arc-transitive, C_uv has a constant value c_+ on oriented arcs u→v and a constant c_- on reversed pairs. Thus C_uv+C_vu=c_++c_- is constant over unordered pairs, and the linear coefficient is constant by vertex-transitivity. Its one-base response depends only on the mask's Hamming weight. This latter stronger conclusion is not asserted from vertex-transitivity alone.

For an explicit infinite family without any external classification theorem, the cyclic tournament on Z/(2r+1) with i→j iff j−i mod (2r+1) is in {1,...,r} is vertex-transitive by translations. For every r≥1 it is distinct from its converse on the fixed labels, but the two have identical zero/one-base responses. This proves that one probe does not reconstruct arbitrary labelled cores in every odd order at least three. It does not establish or refute two-probe reconstruction.

### Explicit nonzero invisible C3 current

Take Q on labels 0,1,2 with arcs 2→1,1→0,0→2 (mask 2 in the supplied edge convention). Then H=3, S=E=(1,1,1), and

C = [[0,0,2],[2,0,0],[0,2,0]].

Hence A=C−C^T is twice the directed three-cycle current: A is nonzero, but every row sum is zero. Thus g=0. The response is

F_Q(b)=3+2Σ_u b_u−2Σ_{u<v}b_ub_v=3+3s−s², s=|b|,

and is identical for the converse. A genuinely ordered two-base experiment separates them: masks (1,2), with first base→second base, give F_Q(1,2)=13 and F_Qop(1,2)=11. These numbers are independently checked by the full-vertex Hamiltonian DP and the common slot realization. This illustrates an invisible circulation made observable by a richer probe, not a theorem transferring lattice higher-genus or electrical-network results into tournament reconstruction.

### Exact low-order collision structure

`structure.py` classifies every one-response collision pair through order five by checking actual core adjacency matrices. At k=4 the eight collision pairs are exactly reversals of a cyclic three-vertex module (all outside vertices have a uniform relation to its members). At k=5 the 92 pairs divide into 80 such cyclic-module flips and 12 converse pairs of regular five-vertex tournaments. There are no unclassified pairs in these finite universes. The counts are consistent with the direct constructions: C(5,3) choices of module, two orientations for the outside edge, and four independent uniform outside/module relations give 80; the 24 labelled regular five-vertex cores give 12 converse pairs. These are FINITE-EXACT statements, not a proof that cyclic-module reversal preserves one-response data at arbitrary ambient order, nor an all-order collision classification.

Reproduce with `python3 structure.py`. The script also checks the converse criterion on all 1,099 labelled tournaments of orders one through five.

## 7. Two ordered ports recover the full ordered-cover sidecar

Status: PROVED, with independent exact controls. This section extends the previous freeze. Use J_uv=C_uv−C_vu for the skew current (the A_uv of Section 6), and use P_u for a vertex potential to keep types distinct. Let e_u denote the singleton attachment mask at u, and define

D_uv=F_Q(e_u,e_v)−F_Q(e_v,e_u),

P_u=S_u+E_u+Σ_w C_uw.

Then the exact two-port contrast is

D_uv=P_u−P_v−J_uv.

### Adjacent-block proof

Name the new vertices a,b. Relabel the new vertices in the second experiment so that in both experiments core→a occurs only at u and core→b only at v. The resulting tournaments differ only by the arc a→b versus b→a. Hamiltonian words not using that arc have identical contributions and cancel. Count those using the consecutive block ab in the first tournament:

- Block at the start: H−S_v, since b may exit to every core vertex except v.
- Block at the end: E_u, since a may be entered only from u.
- Block internal: Σ_{w≠v}C_uw, obtained by deleting the block to leave an ordered spanning two-path cover.

The corresponding ba contributions are H−S_u, E_v, and Σ_{w≠u}C_vw. Their difference is exactly the displayed identity. No inference that deletion leaves a single Hamiltonian core path is used.

### Triangle circulation

For any labels u,v,w, define the triangle curl of a skew array K by curl(K)_{uvw}=K_uv+K_vw+K_wu. Vertex-potential differences telescope, so

curl(D)_{uvw}=−curl(J)_{uvw}.

Thus ordered two-port contrasts reveal the complete triangle-circulation coordinates of the endpoint-cover current. Since triangles span the complete graph's cycle space, they recover its formal cycle-space component. For a vertex-transitive tournament P is constant, and the sharper identity D=−J holds entrywise. These are exact finite combinatorial identities; the Hodge/current language is an interpretation of the proved gradient/cycle decomposition, not a borrowed theorem about physical electrical networks or an implication of lattice higher-genus results.

### Zero-mask calibration recovers the potential too

There is a stronger observation: zero-mask ports recover P, not merely its differences. Namely

P_u=F_Q(e_u,0)−F_Q(0,e_u)=F_Q(e_u,0)−F_Q(e_u).

For the first equality use the same arc-reversal/block proof, now with core→a only at u and no core vertex pointing to b. The ab block contributes H at the start, E_u at the end, and Σ_w C_uw internally. The ba block can occur only at the start and contributes H−S_u. Their difference is P_u. For the second equality, the first base vertex in the word (0,e_u) is a universal source and therefore must start every Hamiltonian word; deleting it gives precisely the one-base experiment e_u.

Consequently zero/one/two-base measurements determine the full ordered-cover matrix and both endpoint vectors. Here are explicit reconstruction formulas. Let

L_u=F_Q(e_u)−H,

B_uv=H+L_u+L_v−F_Q(e_u+e_v), for u≠v.

The one-base polynomial gives B_uv=C_uv+C_vu and L_u=E_u−S_u+Σ_vC_uv. Obtain P from the calibrated two-port responses, then set

J_uv=P_u−P_v−D_uv,
C_uv=(B_uv+J_uv)/2,
S_u=(P_u−L_u)/2,
E_u=L_u+S_u−Σ_v C_uv.

The formulas prove uniqueness and exact recoverability of (H,S,E,C) from this probe battery. They do NOT prove that those data determine Q or that they determine the entire two-base response function. General two-base responses may contain additional ordered-three-path-cover information.

A sufficient fixed battery consists of the empty response, k singleton one-base responses, C(k,2) pair-union one-base responses, k calibrated words (e_u,0), and both orders of every distinct singleton pair. This uses 1+2k+3C(k,2) response entries, each with at most two base vertices. F_Q(0,e_u) need not be queried separately because it equals F_Q(e_u).

### Hostile to a tempting overstatement

The stronger proposed shortcut D_uv=±2H(Q−{u,v}) is false. For the two-vertex core 0→1, C_01=C_10=1 and P_0=P_1=2, so D_01=0. But H(empty)=1, making the proposed right-hand side ±2. The endpoint-cover current, not a plain deleted-core Hamiltonian count, is the correct object. This is the smallest possible distinct-port hostile.

### Controls and freeze

`two_port.py` independently reconstructs all S,E,C values from actual extended-tournament endpoint-DP responses, then compares them to direct core-permutation enumeration. It checks all 75 labelled tournaments of orders one through four, 22 deterministic labelled examples of orders five/six, and vertex-transitive cyclic cores of orders three/five/seven. It also verifies every triangle-curl identity in those cases and the two-vertex hostile. There are 97 base audit cores plus 3 cyclic controls. Normal and optimized outputs agree. The cyclic controls here are not claimed to be Paley tournaments: in particular the order-seven cyclic control has H=175, distinct from the Paley-seven example.

The opposite-C3 witness from Section 6 remains consistent: for its chosen orientation, P=(4,4,4), J_01=−2 and D_01=2, namely 13−11. With the opposite cyclic orientation used by `two_port.py`, D_01=−2 instead. The sign is an orientation convention, not a discrepancy.

## 8. Hamiltonian-polynomial derivative bridge and orientation hostile

Status: PROVED. Introduce independent formal directed-edge variables x_uv for every ordered pair of distinct labels. Define the Hamiltonian-path polynomial

P_k(x)=Σ_{σ∈S_k} Π_{j=1}^{k−1} x_{σ_j,σ_{j+1}}.

This is multiaffine in each directed-edge variable. Evaluate it at the 0/1 adjacency matrix of Q. For u≠v,

C_uv = (∂P_k/∂x_uv)(Q).

Proof: a monomial contributes to this derivative exactly when its permutation contains u immediately followed by v. Deleting that factor cuts the vertex word into two nonempty paths, the first ending at u and the second starting at v. All other factors evaluate to one exactly when these two paths are directed paths in Q. This is a bijection with the covers counted by C_uv. The derivative can be nonzero even when Q_uv=0: the removed variable is no longer part of the evaluation. This is precisely the absent-arc bypass coordinate missed by the refuted insertion-only formula.

Thus the quadratic-size two-port battery recovers the entire first derivative (gradient) of the Hamiltonian-path polynomial at the tournament's adjacency point, together with the Hamiltonian value and endpoint counts. This does not claim that the polynomial's gradient map is injective on tournament adjacency points.

No Hamiltonian-path monomial contains both x_uv and x_vu: using both would revisit u or v and create a two-cycle. Therefore the polynomial is affine under reversal of one core edge. If u→v is the actual edge and Q' reverses only that edge, then

H(Q')−H(Q)=C_vu−C_uv=−J_uv.

This gives an exact finite-difference interpretation of the current: it is an edge-reversal sensitivity of the Hamiltonian count, not the edge's orientation sign.

### Current need not point along the tournament edge

Take the transitive three-vertex core 0→1→2 with 0→2. Its cover matrix is

C=[[0,1,0],[1,0,1],[2,1,0]].

Although the actual nonadjacent edge is 0→2, J_02=C_02−C_20=−2. The current points against that edge in the sign convention J_uv>0 for flow u→v. For adjacent edges J_01=J_12=0 despite their definite orientations. Reversing 0→2 creates C3, increasing H from 1 to 3, exactly −J_02=2. Neither the sign nor support of J reconstructs the core's arcs by itself.

This hostile blocks an unwarranted jump from response-current tomography to tournament orientation reconstruction. The strongest established result here is exact tomography of the endpoint and ordered-cover sidecar, plus finite labelled separation through six; general reconstruction remains open.
