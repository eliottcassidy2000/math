# Session42 synthesis: reconstructing after symmetry, and when shared information survives

23 October2026. All-size theorems below have explicit carriers. Numerical controls are finite certificates, not replacements for proof. No literature-priority claim is made. The session investigates both sides of proposed discrepancies and seeks a larger exact structure rather than retaining only the favorable examples.

## 1. A complementary graph pair is not one tournament class

The [correspondence theorem](correspondence/CORRESPONDENCE_AND_OBSTRUCTIONS.md) preserves the genuine ordered encoding: relative to a supplied vertex order, one bit per undirected pair specifies whether its tournament arc agrees with that order. Graph complementation reverses every tournament arc. Keeping the order therefore gives an exact correspondence.

Forgetting it is a different operation. The same one-edge graph on three vertices can encode a transitive tournament or C3. More strongly, **only empty and complete graphs have observer-independent tournament type, even allowing converse**. Four explicit observer orders expose every other graph's change in cyclic-triangle count.

Any scalar invariant required to be invariant under both ordinary graph relabeling and the transported tournament relabeling must be constant. The generated actions include individual cube-bit flips. Thus a theorem about arbitrary unlabelled graphs cannot simply be transferred to a single tournament class through this chart.

The meaningful replacement has two layers:

1. A connected, weighted many-to-many incidence correspondence between graph classes and tournament classes, with complement/converse pairing if desired.
2. The inherited **signed** Walsh basis of automorphism-even graphs, which gives a faithful linear basis for tournament observables. This is not an object-level bijection from ordinary complementary graph pairs. The signs retain information erased by ordinary graph-class averaging.

Full Fibonacci addresses can supply a real order. A finite color, residue tag or depth cannot replace that unrestricted address merely because it comes from the same numeration system.

## 2. The replacement correspondence has an exact all-size spectral law

The [coupling operator](coupling/COUPLING_OPERATOR.md) alternates conditional expectation onto ordinary graph-invariant and tournament-invariant functions on the uniform edge-bit cube. On graph classes it is a reversible Markov operator whose stationary weights are actual labelled orbit sizes, not uniform weights on isomorphism classes. The paired construction also retains complement/converse.

For a simple edge-pattern graph F, let

mu(F)=E_sigma product_{ij in E(F), i<j} sign(rank_sigma(j)−rank_sigma(i)),

where sigma is a uniformly random total order. The graph-orbit Walsh sum for F is an eigenvector with eigenvalue **mu(F)^2**. This exactly measures the overlap of unsigned graph orbit sums with the inherited signed tournament basis.

The [signed-moment theorem](sign-moments/THEOREM.md) proves

|mu(F)| <= 1/3 for every nonempty F,

with equality precisely for **P3 or C4 plus isolated vertices**. The proof is short and structural: condition on the first vertex. If every deletion leaves an edge, induction applies; otherwise the graph is a star plus isolates, whose moment is an elementary integral. Equality forces the stated two shapes.

Consequently the coupling has **exact spectral gap 8/9 at every n>=3**, paired or unpaired. The n=2 boundary is handled separately. The slow modes are a quadratic wedge statistic and, when n>=4, a quartic four-cycle statistic; n=3 has only the wedge mode. Define the complete signed adjacency matrix by X_ii=0 and X_ij=(−1)^(edge bit ij) for i!=j, not by a sparse adjacency convention. Then the latter statistic is

(sum of signed C4 products) = (tr(X^4)−n(n−1)(2n−3))/8.

This is a genuine second-/fourth-order connection. It is not an identification with the Gaussian quartic invariant below.

The theorem does not imply a size-independent worst-state total-variation time: the initial stationary-mass factor is retained. It is also **not** a spectral theorem for the native tournament metagraph, nor a demonstrated A000568 enumeration speedup.

Information loss remains visible despite rapid mixing. At n=4 the unsigned coupling has rank3 whereas tournament observables have dimension4. An explicit cubic-star observable distinguishes dual source/sink extensions of C3 but is erased by ordinary graph averaging. Converse pairing removes that odd mode. Small-order full paired visibility is not asserted at every size.

## 3. Two congruences per attachment: exactly when the claim is true

For genuine cyclic chart groups C_m and C_n, a nonempty block transporter is an offset plus a subgroup. Its integer lift has unique basis

(h,0), (s,v), with h,v>0 and 0<=s<h.

Thus [each attachment](lattices/README.md) has the exact conditions

y−y0 = 0 mod v,
v(x−x0)−s(y−y0) = 0 mod hv.

Its size is mn/(hv); nontrivial one-sided kernels are retained. These are not necessarily functional attachments.

Across a network, the SAME chart variable must occur at every incidence of a cell. The simultaneous solution set is empty or one lattice coset modulo the chart-period lattice. Its exact count follows from the lattice index, and Smith normal form gives consistency and the composite-order alignment group W/image(D).

This is an all-size rule for a supplied cyclic-chart model. An arbitrary numbering of noncyclic charts does not satisfy it: an actual R9/C3² example breaks the fictitious C9-coordinate description. Full groups instead use subgroup cosets and Goursat's quotient relation. Kernel fibres remain necessary for shared ports.

After authenticating chart choices, the complete bounded integer model can use Session41's coupled primitive-move criterion. Independent residual cycles cannot replace the moves of the full shared-incidence matrix. Actual path positions, converse branches and nonlinear electrical guards remain separately retained unless a defined action encodes them.

## 4. Larger cells: hierarchical visibility does not make ownership cheap

The [hierarchical family](hierarchical-kernels/HIERARCHICAL_KERNELS.md) uses actual regular tournaments on ternary words and prefix-equality cross blocks. Its explicit nested chart groups give a compact description of exponentially many hidden symmetries. Increasing the observed prefix depth reduces the kernel order, but need not reduce the number of actual owners required to eliminate it.

For a 3^k-vertex cell, its full internal group has order3^((3^k−1)/2). A nontrivial prefix depth l<k leaves a one-sided kernel of order3^((3^k−3^l)/2), yet its vertex-owner base remains3^(k−1). With m cells at uniform partial resolution, the exact named-cell owner count is m·3^(k−1); at full resolution it drops to3^(k−1), not zero, because a common diagonal frame survives.

For [mixed attachment depths](hierarchical-kernels/MIXED_DEPTH_AND_SHARED_FLAGS.md) l_e, let G_d retain edges with l_e>=d+1, counting isolated cells, and let c_d be its number of components. The exact named-cell chart group has order3^(sum_{d=0}^{k−1}3^d c_d). Its vertex-owner base is3^(k−1)c_(k−1). This nested-component rule retains the actual embedded equality coordinates, rather than just each edge's number of choices.

For the specified actual concatenated Hamiltonian-path orbit, indexed arc flags must cover the required bottom-triple owner slots. Their minimum is exactly t−nu, where t is the number of slots and nu is the maximum matching size of their nonloop shared-incidence graph; an internal flag supplies a singleton cover. The simpler ceil(t/2) formula holds in the uniform family but can fail after mixed gluing.

An actual **regular order63** example has15 required owner slots but needs **nine indexed arc flags**, not eight. Its full-resolution cell triangle is {1,3,5} inside a cyclic seven-cell host with nine vertices per cell. The slot graph has a matching of size6 and a structural obstruction to a larger one. This is a native shared-port witness, not a relaxed-margin artifact or a Collatz result.

Global regularity requires weighted balance; arbitrary mixed depths need not supply it. Uniform depth, balanced internal cycles and appropriate circulant depth assignments do. These are all-size obstructions to rules based only on coarse connectivity or scalar kernel sizes, alongside compact coordinates for the obstruction. The exact upper bounds retain named cells and the specified path orbit; the unlabelled whole-tournament group is not silently identified.

## 5. Square clocks: unique at odd orders, two choices at even orders

A square-sum path's **square clock** is its ordered list of adjacent square roots. The [reconstruction theorem](SQUARE_WORD_RECONSTRUCTION.md) holds even for arbitrary ordered adjacent sums of a Hamiltonian permutation of {1,...,N}:

- Odd N: at most one realization.
- Even N: at most two realizations.

Every vertex has the affine form x_i=(−1)^i t+c_i. The first moment fixes t for odd N; the second moment gives a quadratic for even N. If two realizations exist with start difference d>0, then 2d divides N, and their alternating position classes are exactly alternating d-long blocks of the integer interval. Their twin map shifts one class by+d and the other by−d.

For unrestricted N=2m permutations, the [independent enumeration](audit-square-words/REVIEW.md) gives

number of double-fibre words = tau(m)(m!)²,
number of distinct words = (2m)!−tau(m)(m!)²,
proportion of permutations with a twin = 2 tau(m)/binom(2m,m).

The normalized double-fibre counts have divisor Lambert series sum_{d>=1} z^d/(1−z^d). These are not formulas for square-restricted Hamiltonian counts.

For square-sum graphs, the first twin occurs **exactly at36**. Fifty-two smaller divisor-block candidates have explicit degree/connectivity obstructions. Both order36 paths are retained as witnesses.

The [palindromic extension](PALINDROMIC_SQUARE_CLOCKS.md) has a different sharp threshold: **48**. A palindromic clock makes reversal equal the twin. Quotient by the twin pairs, but retain the marked central endpoint. At36 a two-vertex separator produces four components; at40 the possible central mark set is empty. The order48 witness passes every edge and reflection check.

Thus an exact two-coloring becomes meaningful through reconstruction, not through a guess that all square-sum edges should be properly two-colored.

## 6. The three islands are ancestry sets, not permanent components

Direct square-sum edges toward the smaller endpoint. The only founders are1,2,4. The [ancestry theorem](square-colors/RESULTS.md) proves that the complete pure-founder islands are

{4,5,11,12},
{2,7,9,16,18},
{1,3,6,8,10,15,17,19,30}.

The user's first two proposed islands are complete; the third gains15,19,30. They are not the later graph's connected components. Vertex13 already has mixed ancestry,22 has all three, and **every n>=71 has all three founders**;70 is the last exception. A finite base71..99 plus a largest-predecessor induction proves the all-size tail.

A coherent single allegiance is obtained by choosing one earlier parent at each birth. This gives a three-root forest. Its exact joint partition function retains multiplicity: through25 there are1024 parent-selector forests but192 colorings. Uniform forests do not give uniform colorings, and positive individual marginal probabilities need not be jointly possible.

The birth DAG has merges and unbounded indegree, unlike a Berggren tree. Therefore no bijection preserving all its predecessor incidences identifies the two. The inherited Fibonacci-to-Berggren three-ray map remains a genuine positive connection, with its normalization and exact address retained.

## 7. A mass coloring experiment, with a test against overinterpretation

The same packet fixes procedural vertex/edge feature grammars using parity, square roots, digit sums and Zeckendorf features. It tests3104 formulas, or1584 distinct colorings up to color names. Of these,1549 decode the supplied complete Q32 carrier of392 paths when actual endpoints and ordered color words are retained.

That high success rate is evidence against singling out one coloring as profound merely because it passes a small test. An exact equality-constraint rank computes each random-color collision probability, and gives sufficient two-color decoding bounds for this finite path family. An inclusion-exclusion partition function retains shared equality constraints.

The [independent complete holdout](color-holdouts/RESULTS.md) then enumerates all unoriented paths at orders33,34,35,36:669,4041,17175,12960 respectively. Without retuning the feature grammar, only1263 of the1549 training successes survive every held-out order. Survivors are4/56 vertex-two-color,15/60 edge-two-color,281/359 vertex-three-color and963/1074 edge-three-color candidates. Every family has failures already at33, despite passing the earlier two selected hostile pairs. Independent enumeration and rescoring reproduce these losses.

There is also a positive all-size palette map: edge color e(a,b)=c(a)+c(b), together with an actual starting vertex, makes the ordered edge-color and vertex-color words mutually reconstructible. It preserves their endpoint-conditioned collision partitions. For connected graphs an odd cycle kills the global lifting ambiguity over F3, but a constant-flip ambiguity remains over F2. This is a proved incidence-algebra distinction, not a proposed allegiance identity.

Proper coloring, allegiance and path observation remain different tasks. Proper bipartiteness already fails at15, while the first triangle occurs at30. No eventually periodic coloring properly colors the entire infinite square-sum graph: m²t² and3m²t² have the same residue modulo m and sum to(2mt)².

## 8. Fibonacci colors give a direct additive bridge

Put alpha=phi^(-2). The [additive packet](fibonacci-addition/RESULTS.md) identifies the canonical colors with rotation intervals of lengths alpha,1−2alpha,alpha. For every ordered color pair c,d and target N, its representation count equals (N−1) times an explicit circular interval convolution at phase {Nalpha}, with a uniform logarithmic error. The Fourier series has literal sine factors for interval widths and cosine factors for the retained target phase.

For canonical A, the proportion of square roots k for which k² has two distinct A summands is exactly **2alpha**. Its normalized unordered representation count has limiting mean alpha²/2. The [full distribution](fibonacci-addition/SQUARE_COUNT_DISTRIBUTION.md) has an atom1−2alpha at zero and density4 on(0,alpha/2). All nine color-pair counts share one phase; factorization of their average does not make them independent.

Every permitted branch retains the forced canonical A positions, so it inherits lower bounds, not necessarily the same limits. At the opposite end, the legal selector B→A and C unchanged produces A_max, which represents **every integer N>=7 by two distinct summands**, with last failure6. This fixed selector uses only A and C and does not inspect primality or the target; its role is additive representation, not proper graph coloring.

For one target, the exact branch representation polynomial is

2^u z^r_AA (1+z)^r_AB (3+z)^r_BB,

where the r values count canonical unordered distinct pair types and u counts unused B choices. It describes every intermediate attainable representation count. Several targets share selector variables and cannot be combined by multiplying separate marginal polynomials.

Thus both branches have genuine, different additive roles. None is a proper square-sum coloring, since forced A already contains2 and7. Prime-weighted representations require new estimates; ordinary additive abundance is not Goldbach.

## 9. Gaussian one-eighth compression and its limits

The [Gaussian packet](gaussian/RESULTS.md) gives the exact one-octant representative count, including origin, axis and diagonal stabilizers. The generic orbit size8 is not an exact N/8 count. No larger origin-fixing Euclidean lattice-isometry group exists beyond D4.

There are valid stronger decoders: primitive rays with dilation, arithmetic norm-shell generation, and the ring generator i. They reduce stored descriptions by changing the operations supplied to the decoder. For arbitrary marked subsets, the discarded memberships still have to be recorded.

The complete single-point orbit coordinates are u=x²+y² and v=x²y². The first norm shell with multiple octant orbits is25, containing(5,0) and(4,3). Its missing quartic angular coordinate satisfies cos(4theta)=1−8v/u².

Even these complete separate orbit coordinates do not determine addition of two inputs. The relative frame controls Re(z conjugate(w)); independently replacing1 by−1 changes the norm of the sum. This is an exact group-action example of the same shared-alignment issue, not a transfer of an arithmetic conjecture.

Bott's real period8 is a stable homotopy/Clifford-algebra theorem, not the order of D4. The packet preserves that distinction and the previously proved prime-local golden/Gaussian ring bridges.

## 10. New incoming work is credited and changes the frontier

Incoming43938b527ec, af23ceb861d and f7abf1477ea were preserved during publication. The [bounded integration review](late-incoming/INTEGRATION.md) distinguishes what was rechecked from inherited finite censuses. Their [cyclic attachment and charge packet](../cyclic_attachment_lattices_20261023/README.md) independently develops exact coset gluing and integrates the Session41 branch language and Session42 shared-frame results.

Its restricted-character theorem tests constancy of supplied local observation tables on the SAME solution coset. Equal restricted characters must be combined before testing cancellation. A generator-line hostile proves that integer representatives cannot be treated as additive group characters: modular carries can remain invisible on each separate generator line. The number of Fourier terms is not by itself a bit-complexity bound for large cyclotomic fields.

The nonabelian prefix bank uses complete certified coordinate images, not sampled witnesses. Generic authentic finite-group coset feasibility is covered by the cited Mal'tsev theory in explicit domains; the packet does not implement a general solver or turn implicit large groups into small inputs. Its full-branch Fibonacci charge ledger is conservation, not a debt bound or descent proof.

A controlled UNION of cosets for native path observations is a worthwhile next hypothesis. No universal bound on its number of pieces is proved.

## 11. The underlying structure and remaining obligations

The repeated useful operation is **retain common incidence data before taking a quotient or marginal**. It appears here in group transporters, Gaussian addition, ancestry choices, graph/tournament Fourier signs and path-clock reconstruction. Those are explicit parallel constructions, not a claim that all their objects are isomorphic.

Concrete remaining obligations are:

- A compact complete native-port model at unrestricted tournament orders, beyond supplied named chart systems.
- Exact compatible observers when path choices, electrical guards or larger-cell kernels introduce non-coset filters.
- Enumeration gains established by an actual algorithm and measured against existing A000568 methods; a fast quotient chain alone does not count all classes.
- Color tests on independent complete larger carriers, not just selected hostile pairs.
- The stronger globally reused noncanonical Fibonacci expansion dictionary, distinct from occurrence-dependent branch choices.
- Prime-weighted additive estimates if pursuing Goldbach, and a genuinely uniform scale argument for the separate Collatz obligations.

All new packets preserve explicit source pins, normal/optimized controls and independent audits. No open conjecture is declared solved by resemblance or finite testing.
