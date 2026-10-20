# Merged-metagraph foundations and a closed transitive-pipe gap

Status: recovered theorem map; new elementary mixed-landing proof; exact n=5 hostile replay. Source snapshot: math commit129154b04024. For the stronger pure-blue classification reduction proved during this investigation, see PURE_BLUE_STRONG_CYCLE_REDUCTION.md.

## 1. Exact objects and operation firewall

Use vertices 0,...,n-1 and fixed Hamiltonian path n-1 -> ... -> 0. Tile positions are S={(x,y):x-y>=2}, m=binom(n-1,2), with t(x,y)=1 meaning x->y. The staircase cube X=F_2^S contains all tournaments with that fixed path.

* All-arcs converse reverses all binom(n,2) arcs, including the path; it does not preserve this fixed-path gauge.
* Path-preserving antipode a(t)=t+1 reverses exactly the m nonpath arcs. For n>=3 its free two-orbits are the literal complement-line instances L=X/<a>.
* Grid reflection r(t)(x,y)=t(n-1-y,n-1-x) permutes bits, with no bit negation. Its tournament is the all-arcs converse relabelled by vertex reversal. Thus it preserves the converse-merged node q(t).
* M consists of tournament isomorphism classes modulo all-arcs converse. Its fibre F(u) consists of the normalizations of all Hamiltonian paths in its one or two tournament classes, modulo automorphism.
* Blue means r(t)=t. Since a and r commute, both ends of one line are blue or both black. A line boundary is the unordered multiset {q(t),q(at)}. Parallel lines remain distinct. A loop consumes two tilings/half-edges.
* The local one-tile/Hamming-neighbor graph and the arbitrary-one-arc isomorphism metagraph are NOT this antipodal line graph. Even claims identifying their supports require their own theorem; THM-1430 labels the all-order avoidable-arc identification conditional.
* 'Points of symmetry' in this session means pure-blue nodes, exactly as the user specifies. Historical uses for ambiguous score classes are not silently substituted.

The faithful object is the span M <- X -> L, not just the simple projected graph (THM-796-three-sorted-recursive-tiling-line-node-incidence).

## 2. Strong inherited laws

Write H(T)=number of directed Hamiltonian paths, A(T)=|Aut(T)|, N_u=|F(u)| and g_u=number of blue tilings in F(u).

1. THM-781-hamiltonian-path-inverse-metagraph-fibre gives the exact inverse F(u)=disjoint union HP(T)/Aut(T), so N=H/A on SC nodes and 2H/A on a non-SC converse pair. H and A are odd. Consequently N is odd on SC, and N=2 mod4 on non-SC (THM-345-merged-bucket-parity).
2. THM-643-blue-black-line-parity-structure proves r acts on each odd SC fibre and hence g is odd and positive there; non-SC fibres have g=0. Pure-black is exactly non-SC. The black half-edge mass N-g is even in every merged node.
3. Let f=floor((n-1)/2), h=(m+f)/2=floor((n-1)^2/4). Global masses are |X|=2^m, |L|=2^(m-1), #blue tilings=2^h, #blue lines=2^(h-1).
4. More strongly, r fixes every blue line and pairs black lines without fixed points, preserving each unordered node pair AND exact reflection defect delta=t+rt. Thus every black node-pair line multiplicity, including loop multiplicity, is even, separately at every defect (THM-796 §4). For the oriented half-edge matrix A^c(u,v), diagonal entries count loops twice, so black diagonal entries are multiples of four.
5. THM-346-tiling-quotient-bucket-balance: 2 internal_lines + cross_halfedges = fibre_size * number_of_masks for every quotient and mask family. It does not assert lumpability. THM-796 §8 exhibits recursive non-lumpability from n=5.
6. THM-646-line-metagraph-score-complement-skeleton gives pointwise s_v(t)+s_v(at)=c_v, with c=(n-2,n-1,...,n-1,n) in the descending convention.
7. THM-649-flip-parity-trilogy-affine-D-binomial-law / THM-790-blue-parity-law-proved give Delta c3=e_(n-1)-e_0, where e are nonpath outdegrees at the two endpoints. Only the two staircase legs and apex carry current. On blue tilings e_0+e_(n-1)=n-2, hence Delta c3=2e_(n-1)-(n-2), odd at odd n, even at even n, with absolute value <=n-2. Its full blue distribution is binomial, multiplied by the silent interior cube size. No blue loops occur at odd n.
8. THM-785-cyclic-triangle-flow-on-blue-black-metagraph retains endpoint defect epsilon=s_0+s_(n-1)-(n-1). Reflection reverses epsilon while preserving merged endpoints and current orientation. Its absolute magnitude is a genuine sidecar; its signed bins cancel even conditionally on one node pair.

Most structurally powerful inherited package: THM-796's three-sort incidence, exact defect short exact sequence, face phase torsor and coupling tensor; THM-830-b3-deletion-deck-mirror-current-calculus's gap-tournament/pair-groupoid dictionary; THM-851-node-coloured-defect-coherent-closure's factorization deck. THM-851 proves its closure lemma for all finite coloured groupoids and proves the binary composition tensor at all sizes, but full reconstruction is only FINITE-EXACT n=3,5,6,7 (n=4 is a genuine failure; n>=8 remains open).

## 3. Concrete count-preserving impossible rewiring

The independent script probe_n5_rewiring.py enumerates all 64 tilings, canonicalizes by all 120 vertex permutations and converse, and reports the 10 merged nodes. Tile order:

(2,0),(3,0),(3,1),(4,0),(4,1),(4,2).

Two actual blue lines are {0,63} and {4,59}. Their four distinct endpoint nodes have c3 values (3,0) and (4,1), respectively. Replace their incidences by (3,1) and (0,4), retaining those same four nodes; keep every other blue and black line unchanged.

This 2-switch preserves every individual node's blue degree, black degree, fibre size, colour type, blue eligibility, loop count, all global counts and all parity/handshake constraints. Yet it is impossible: at n=5 every blue line has odd c3 jump of magnitude <=3. The new jumps 2 and 4 violate this, and 4 violates the maximum too. Thus count compatibility is strictly weaker than realizable arrangement. The witness is not merely an abstract graph; it is a two-edge surgery on an exact actual n=5 atlas.

## 4. New theorem: the unique transitive pipe always lands mixed

**Theorem.** For every n>=4, the unique complement line incident to the transitive merged node is blue and has a mixed other endpoint. At n=3 its other endpoint is the pure-blue 3-cycle. Its c3 jump is n-2.

**Proof.** The transitive class has exactly one fixed-path tiling, all ones: H=1 and |Aut|=1. Its antipode is U_n, the all-zero tiling, with all adjacent arcs descending and all nonadjacent arcs ascending. Both constant tilings are reflection fixed, so their line is blue. U_n is not transitive for n>=3; its cyclic triangles are precisely the consecutive triples, giving c3(U_n)=n-2 (or use the inherited leg law).

To prove mixed landing, take the Hamiltonian path

P=(0,n-1,n-2,...,1)

in U_n. Its first arc is nonadjacent ascending and all subsequent arcs are adjacent descending. Both endpoints have outdegree n-2. A reflection-fixed normalization of P would make positional reversal an anti-automorphism, forcing the endpoint scores to sum to n-1. But their sum is 2n-4, unequal to n-1 exactly when n!=3. Therefore P normalizes to a black tiling in the same class for every n>=4, while all zero is a blue tiling there.

There is an explicit mask witness with no isomorphism search: relabel old j to j-1 modulo n. This sends P to the fixed descending path, and the normalized mask has support

{(n-1,y):1<=y<=n-3}.

Its reflection is the disjoint other-leg support {(x,0):2<=x<=n-2}; these nonempty masks are distinct for n>=4. Hence the landing fibre contains at least one blue and two black tilings. At n=3 U_3=C3 and its fibre is a singleton. QED.

This closes the named general mixed-landing gap in THM-790 §§5–6 / THM-785 §3, which stated only certified n=4..8 evidence. verify_transitive_mixed_landing.py directly checks the path, normalization support, reflection, scores, and c3 for n=3..100. The proof is all-order and does not depend on those finite checks.

## 4a. Strong core and nonstrong periphery

Every complement partner of a nonstrong tournament is strong for n>=3. A Hamiltonian path enters strong components in their unique linear order, so its first vertex beats its last. Complementing the nonpath endpoint arc closes the preserved path to a Hamiltonian cycle. Thus the nonstrong merged nodes form an independent set with no loops. This is an explicit corollary of the inherited apex-closing mechanism THM-333-apex-tile-strong-connectivity. It does not depend on blue/black colour or node counts. The companion cycle-reduction note combines it with even purity=rigidity to show every even-order point of symmetry at n>=4 is a blue leaf attached to a mixed node. STRONG_CORE_COUNT_RECURSIONS.md gives all-order recurrences for strong and strong-blue tiling counts and the resulting exact colour-resolved cut masses.

## 5. Correction lineage and current purity frontier

* MISTAKE-031/033: complement tiling is not all-arcs converse; MISTAKE-065: reflection is a coordinate permutation, not bit negation. MISTAKE-081: self-converse counts are not ordinary tournament counts at n-1. Avoid stale conflicting OEIS IDs in HYP-3809; the mathematical Burnside formula |M_n|=(#iso classes+#SC classes)/2 is unambiguous.
* HYP-3809-merged-metagraph-bucket-conjecture-atlas's parity and blue-fold targets are superseded by THM-643/647/796. Its even multiplicities of blue-count values are finite data, not a canonical involution. HYP-3811-twin-separator-and-the-spectral-twins-platonic explicitly says flip*r does not descend to the claimed node twinning. Its 'ONLY spectra distinguish' wording is also too strong: committed merged_metagraph_twin_separator_kps.out separates both named mixed twins with c5/cycle-polynomial coordinates too.
* Pure-blue does NOT mean rigid. Rigid SC (H/Aut=1) implies pure-blue, but regular cyclic T5 has H=15, |Aut|=5 and fibre/blue mass 3. Consequently proving H/Aut>1 does not prove mixed. THM-2477's claim that its non-rigidity corollary settles a rotational purity question must be narrowed.
* T5 is not a Paley tournament or a quadratic-residue tournament modulo5. Call it the regular cyclic T5; the standard Paley tournament construction requires q=3 mod4.
* At odd orders, pure-blue does NOT mean leaf/pendant. In our exact n5 coordinates, the regular pure-blue node (canonical code76, representative12, fibre/blue mass3, c3=5) has three distinct blue neighbors: codes40,41,12. This refutes the pendant prose in HYP-3808, THM-643 and THM-2444 §5. The singleton-fibre rigid points are leaves; the nonrigid T5 point is not. The companion cycle-reduction theorem proves that every even-order pure-blue point is rigid, and hence a leaf attached to a mixed node, for all even n>=4.
* THM-2444 refutes floor((n+1)/2)-[n even] at n9 (6, not5). THM-2454 also refutes '+one nonrigid at every odd size' at n11: two nonrigid pure-blue classes occur, including C3*T5*C3.
* THM-2454-pure-blue-alphabet-and-center-law proves: for a palindromic SCC stack, blue mass = product(half fibre sizes)*blue(center), while total fibre size = product(half fibre sizes)^2*fibre(center). Purity is equivalent to rigid halves and pure center. Its complete atom alphabet {1,C3,T5} is FINITE-EXACT through n13, not proved all-order.
* THM-2453 proves the palindrome/Narayana formula within the {1,C3} grammar. THM-2477-rigid-strong-classification-through-sixty-two supplies an inherited claimed proof of rigid strong=C3 through62 and modular reduction all-order; that source leaves unrestricted strong rigidity beyond62 open. This is the recovered source boundary, not a claim about results outside this packet. This does not classify nonrigid pure atoms.

## 6. Least-used useful sidecar

Preserve the defect-labelled factor-incidence deck, not just a new scalar invariant. THM-851 gives explicit n5 masks8 and24 with the same exact defect01 and same merged node but different factorization counts. Its full ordered deck separates every tiling at n5..7; its unordered reversal-compatible deck separates reflection orbits. Forgetting factor defect labels leaves residual collisions even after retaining endpoint nodes and unlabelled midpoint node pairs. The surviving field is (factor defect) <-> (pair of factor nodes).

## 7. Five worthwhile next targets

1. Pure-blue strong atoms via the Hamiltonian endpoint/cycle reflection obstruction, not another unstructured census. The companion proof PURE_BLUE_STRONG_CYCLE_REDUCTION.md closes the even-order obstruction using Lu's published crossing theorem: every nontrivial strong pure-blue atom is odd, regular and circulant, and every even pure-blue node is rigid. Remaining target: classify the odd circulant atoms, potentially using the injection of Hamiltonian cycles into n-cycle automorphisms.
2. Seek a genuine explanation for the odd-n blue-count spectrum's equality on the two c3-parity sides, newly observed independently at n3,5,7. Anti-Redei alone gives odd values but does not force sidewise equality. First distinguish a distribution identity from a node involution.
3. Defect-refined black quotient Q=A_black/2, with loops normalized correctly: derive cut/charge constraints and test whether bridges or components obey sharper restrictions than Eulerian parity in the doubled graph.
4. Continue THM-851 at the first untested size with controlled-forgetting probes: does one defect-labelled factor deck still suffice at n8, and what is the minimal collision if not? This is a faithful arrangement test, not a scalar-count exercise.
5. Use the node/line/face coupling tensor, especially the straight/crossed half-edge phase. Bare node deletion rows are not Markov/lumpable even when they distinguish all nodes; characterize the smallest continuation-equivalent refinement rather than declaring a complete static fingerprint a recursive state.

## Reproduction

python3 probe_n5_rewiring.py
python3 verify_transitive_mixed_landing.py

Matching outputs are included beside their sources. The separate strong-core count note and verify_strong_core_counts.py add an independent literal reachability/mirror audit through n7.
