# Session43 synthesis: observers, restricted gaps, arithmetic lifts and exact compiler boundaries

23 October2026. This session treats each proposed connection as a map with a declared source, action, measure and retained coordinates. A failed identification is replaced, where possible, by a precise positive construction. All-size proofs, finite computations and inherited claims remain distinct.

## 1. Every tournament supplies an observer, but not a canonical one

The [native observer packet](native-ports/OBSERVER_BUNDLES_AND_NATIVE_PORTS.md) formalizes (T,P), with P an actual ordered Hamiltonian path. Automorphisms act freely on such paths, so the observer fibre over a tournament class has size H(T)/|Aut(T)|. Normalizing P gives the cube of (n−1)(n−2)/2 free nonpath arcs.

Thus existence supplies observer order, exactly as suggested. It does not select one actual order naturally when Aut(T) is nontrivial, and it changes the sampling weights: uniform observer masks weight classes by H/|Aut|, not uniformly.

Normalized converse reflects free bits; it does not complement them, because full arc reversal and reversal of observer order each contribute a complement. Free-arc flip and reflection commute, giving exact blue/black budgets. These recover the established observer dictionary rather than claim to determine all cross-class wiring.

The crucial native distinction is:

- Complete tournament isomorphisms automatically transport genuine paths.
- Partial indexed-owner constraints on a complete transporter remain coset constraints.
- Assembling paths from internal cell charts before whole-tournament compatibility is known can require unions of cosets.

An actual almost-regular six-vertex example has six paths in a specified two-cell family, arranged in two diagonal cosets, while the host has45 Hamiltonian paths in total. The family is not silently identified with the whole path set.

## 2. Which restrictions preserve the 8/9 gap?

The [restricted-coupling theorem](coupling/THEOREM.md) concerns Session42's graph/tournament averaging correspondence, NOT the native metagraph walk. Its graph labels use full backward graphs on n vertices; complement/converse pairing means the full graph-complement and tournament-converse quotient labels. Compressed backward graphs on n−1 vertices form a different incoming carrier, discussed below.

Reducing invariant subspaces retain their eigenvalues. Coarsening either observation while keeping the same uniform edge-cube measure also preserves gap at least8/9. Equality requires the surviving P3/C4 mode and its corresponding projection to survive. Observing H modulo8 as a coarse statistic is such a valid operation.

Conditioning on a residue or a Hamiltonian observer is different:

- Fixed Hamiltonian-path conditioning isolates the transitive/empty-graph tile for every n>=3. The whole conditioned chain has gap0, even after complement/converse pairing.
- Removing that tile at order4 leaves a seven-tile component with exact gap7/15.
- At order5, conditioning on H≡7 modulo8 gives gap226/315, below8/9. Conditioning on residue5 instead gives gap1657897/1801800, above8/9. These are exact rational calculations; smaller nonempty residue carriers have no such failure.
- Refining both sides with a common named datum can create a nonconstant common observable and hence disconnect mixing.

There is a positive robust replacement. If a new joint weight obeys 0<a<=w<=b, the new conditional-expectation chain has gap at least8a/(9b). The proof compares conditional variances and retains the normalization. Hard filters have a=0 and lose this guarantee. A finite rational positive-semidefinite criterion certifies any explicitly supplied restricted incidence matrix.

## 3. Native branch structure: exponential coset lists, small transfer states

For a contiguous cyclic-cell template, the native packet gives an exact all-size path count p times the product of allowed boundary-palette sizes. Interleaved templates require shared phase differences whose sums vanish around every incidence cycle. For each fixed consistent boundary-difference choice, a partial forest retains p^c lifts; a single spanning template is connected and retains p per consistent choice. Indexed owners can fix those lifts. A label-sensitive mod8 filter must be checked on the actual lifts unless its constancy on a quotient has been proved.

The [branch-width theorem](native-branch-width/THEOREM.md) sharpens the representation boundary. There are actual regular order-3t hosts, t odd, whose declared path families have3·2^(t−1) members and require exactly2^(t−1) affine cosets in ANY coset cover. Projecting to phase differences produces the binary box {0,1}^(t−1), which contains no nontrivial subgroup image over F3; every contained coset therefore fits in one three-element diagonal fibre.

Nevertheless exact counting and endpoint observations require only a three-state transfer matrix. The endpoint correction is period six; fixing only the endpoint difference still leaves three common-phase lifts, while fixing the initial owner removes them. For odd prime p and a nonempty proper palette E⊊F_p, the corresponding version retains p states and exact coset width |E|^(t−1).

This answers the larger-universe question constructively: coset-list size, transfer width and native completeness are different resources. Exponential coset width alone is not a computational-hardness proof.

### The harder boundary is now a proved native decision problem

The [native complexity theorem](native-complexity/THEOREM.md) proves NP-completeness of a precisely constrained problem: given an actual tournament, its supplied directed-C3 cells, and a spanning template, choose one shared cyclic shift per cell so that the template becomes a Hamiltonian path. The partition can be equitable and every cross block circulant or uniform.

The reduction retains actual vertices and arcs. Equality chains share each graph-color variable, boundary palettes encode equality or inequality, and an Euler trail plus local cyclic padding produces a literal spanning template. This is not a claim that finding an unrestricted Hamiltonian path in a tournament is hard.

There is also a pinned regular-host version. Double the actual host using two copies with the same internal tournament, a center with A→center→B, and cross arcs A_i→B_j exactly when j→i in T for i!=j, with B_i→A_i on the diagonal. The resulting order2N+1 tournament is regular. Append a genuine fixed suffix and retain its actual owner pins. Feasibility in the allowed template family remains exactly the original problem.

Consequently a universal polynomial-time solver for this supplied family would imply P=NP. This does not prohibit exact all-size rules: bounded incidence treewidth gives a three-state-per-variable dynamic program, and fixed branch parameters give other exact methods. It also does not transfer hardness to complete tournament-isomorphism transporters or to an unspecified native metagraph task. Additional historical electrical/mod8 guards require their own preservation proof.

## 4. Finite mod8/mod9 frames and Gaussian moments

The [Gaussian phase packet](gaussian-jacobian/THEOREMS.md) proves a precise obstruction and repair. For a finite group of size N acting on a reduced commutative Q-algebra, vanishing of the first N Reynolds-averaged powers forces the polynomial to vanish. Newton identities prove this, and cyclic examples make the bound sharp.

By contrast a nonzero circular complex Gaussian variable Z has E[Z^m]=0 for every positive m. Finite orbit cancellation therefore cannot stand in for infinite angular cancellation at every degree.

The positive repair retains integer charge. Cyclic averaging P_L keeps charges divisible by L; continuous circle averaging keeps only charge zero. They agree whenever the authenticated nonzero charge bandwidth lies strictly below L. Exact support propagation can give a sharper test.

P8P9=P72. The first alias for the powers of Z occurs at72; finite-angle radial sampling then gives36!, whereas the true pure Gaussian moment remains zero. Shared total charge must be combined before averaging: Z_1 conjugate(Z_2) survives a common phase but is lost under independent phases.

The [bounded-degree planar repair](FINITE_GUARDS_AND_DEGREE.md) is equally positive and sharp. Exact C_L equivariance for weights(1,−1), with degree D and L>D+1, forces full C* equivariance. A Keller map under that full action is diagonal linear. Consequently simultaneous exact C8 and C9 symmetries certify linearity through degree70. The degree71 automorphism (x+y^71,y) passes both finite symmetries but is nonlinear. It is not a Jacobian counterexample.

Provenance matters: the primary GMC3 counterexample and its dimension-two limitation were checked. The repository's GMC2/Lean endpoints were inspected, not rebuilt or independently certified here. None of this session's theorems establishes general planar JC or Gaussian moment conjectures.

## 5. The 4→2→5→3→7→13 cycle has a literal arithmetic repair

The [cycle and geometry packet](cycle-geometry/PROOF.md) recovers the Monster/Baby Monster exponent-difference source. A modulo9 interpretation fails as a full quotient:11 and2 share a residue but have different outputs.

Instead extend prime valuation to composite bases:

h_N(b)=max{k:b^k divides N}, b>=2.

For the specified integer N=2^5·3^7·5^3·7^4·11·13², the cycle

4→2→5→3→7→4

is literal, and13 also maps to2. Folding4 and13 preserves their forward transition, but not incoming histories, prime status or arithmetic labels.

The exact divisor-height polynomial is

5616z+120z²+16z³+4z⁴+2z⁵+z⁷.

Among the5759 divisor bases greater than1,143 enter the cycle and5616 terminate at1. The map at1 is not assigned a finite valuation.

There is a complete finite compatibility theorem. For prescribed bases b and heights k_b, set L_p=max_b k_b v_p(b) and N0=product p^L_p. A realization exists exactly when h_N0(b)=k_b for every prescribed base. Then N0 is the least realization in divisibility order, allowing N0=1 in the all-zero case. Shared prime coordinates, rather than independently chosen heights, determine compatibility. This yields entire families of cycles with an exact construction test.

## 6. A sphere, a logarithmic cover and a hyperbolic orbifold are different structures

The five-cycle can be realized as the monodromy of w^5=z, but that preserves the cyclic sheet action, not integer valuation arithmetic. A separate reversible completion of the side arrow11→1 gives a seven-state passport(5,2). Adding a joining branch leads to the explicit double-star polynomial

P_(a,b)(z)=(a+b)^(a+b)/(a^a b^b) · z^a(1−z)^b.

For(a,b)=(5,2) it has degree7 and monodromy S7. The completion and joining branch are stated choices, not consequences of the Monster.

Its source is a Riemann sphere, yet its branch-order orbifold(10,2,7) is hyperbolic. These statements concern different retained structures. The inherited affine triangle-group family(2,3,2p) separately gives a genuine spherical/Euclidean/hyperbolic transition at p=2,3,>=5, with its bad-prime and lattice guards retained.

The complex plane itself is not intrinsically a spiraling multivalued surface. The exponential cover makes logarithmic winding single-valued upstairs; a negative-real cut chooses one branch downstairs. Keeping the sheet coordinate is the precise analogy with retaining a kernel lift or integer phase, rather than identifying the objects.

## 7. Tournament path covers give a sharp all-order binary law

The [cyclotomic packet](cyclotomic/NOTE.md) starts from the inherited directed path-cover identity

F_T(x)=sum_m m! u_m(T)(x−1)^(n−m),

where F counts forward steps over all vertex orders. For a prime p and pi=zeta_p−1, every term has pi-valuation at least n−s_p(m), where s_p is the base-p digit sum. This proves an all-size bound using n−max_(1<=m<=n)s_p(m). Palindromicity forces even valuation for odd p, strengthening the rounded bound. At p=3 the cyclotomic value is a root-of-unity multiple of an integer, giving an ordinary power-of-three divisibility statement and recovering the mod9 law.

At x=−1, even orders vanish. For every odd n>=3 the exact largest UNIVERSAL binary divisibility exponent is

B(n)=min(n−v2(n+1), n+1−floor(log2(n−1))).

Every tournament's value is divisible by2^B(n), and an explicit tournament attains that valuation at every odd order. The exceptional formula branch is the binary-repunit orders. The transitive tournament is an extremizer exactly at orders2^r−1 and3·2^j−1 within this range; other orders have explicit reversed-extreme-arc witnesses. Classical tangent valuations and the classical extremizer family are credited.

This concerns F_T(−1), not H(T) and not the native electrical mod8 guard. The Gaussian radial identity E|Z|^(2m)=m! supplies an exact finite moment packaging of the path-cover formula; it does not turn it into a pure-power Gaussian nullness statement.

## 8. All universal linear congruences are computable without enumerating tournaments

The [affine polynomial compiler](cyclotomic-lattice/NOTE.md) proves that the integer affine hull of all order-n forward-step polynomials is

A_n + sum_(s=1)^floor((n−1)/2) Z · [2(x−1)^(2s) A_(n−2s)],

where A_n is Eulerian. The inherited odd-cycle expansion gives one inclusion; a unitriangular family of actual single-backward-arc tournaments gives the other.

Therefore every universal linear congruence, and every universal cyclotomic valuation minimum, can be tested on only O(n) explicit tournaments. The odd-prime closed-form equality suggested by505 exact order/prime tests remains a conjecture; the finite-generator minimum algorithm is proved for all orders.

This is not a realizability theorem. At order5,

A5+3C1 = [7,38,30,38,7]

has nonnegative integer coefficients, palindrome symmetry, total120, nonnegative cycle-excess coordinates, and satisfies every universal affine-linear congruence. Yet it requires three directed triangles and no directed5cycle, which is structurally impossible. The all-order excess-coordinate obstruction and a recovered proof of the known H≠7 gap are retained separately. This does not settle the complete Hamiltonian-count spectrum or explain21 by a numerical analogy.

## 9. Incoming work supplies a different positive restriction theorem

Incoming e858f56b332 and8a73412b983 were preserved. The [bounded integration review](late-incoming/INTEGRATION.md) records a useful distinction rather than a contradiction.

The incoming [observer-restriction packet](../hamiltonian_observer_restrictions_20261023/README.md) uses the COMPRESSED graph on n−1 vertices, mapping its edge(i,j) to a backward arc(i,j+1). Our conditioned coupling uses the FULL backward graph on n vertices. The labelled bit sets are bijective, but their unlabelled partitions are incomparable from order5. Their complement pairings already differ at order3.

The incoming all-size two-component descent theorem is therefore a theorem for its compressed incidence carrier, not an automatic proof for our full-backward carrier. Its boundary repair and the distinction between arbitrary and reference-pointed selectors are retained. Our full-backward nontransitive component is connected in the exact tests through5; the corresponding all-size extension remains open here.

Likewise, the incoming shared-coordinate frame cube counts complete automorphisms of a supplied actual path. It does not enumerate all independent-chart path assemblies, so it is compatible with our union-width and constrained-complexity results.

A separate incoming trace construction preserves a gap bound at least8/9 by observing RETURNS to an accepted set through an auxiliary ambient chain. This is not the kernel rebuilt by hard conditioning. Its stationary mean duration is1/pi(S), and intermediate excursions need not preserve the native source and ports. The positive theorem retains return clocks and still requires a native lift/cost proof. It cannot be substituted for the disconnected conditioned chain or its7/15 component.

Thus five operations must remain separate: observer normalization, choice of graph quotient, coarse observation, hard conditioning, and guarded return sampling with its clock. The incoming results were scope-reviewed and supplemented by small independent carrier witnesses, not fully reenumerated here.

## 10. Remaining boundaries

- A native family compiler must prove its state inventory complete for the family actually claimed; complete isomorphism transport and independent path assembly are not interchangeable.
- Conditioning and guard observation require different spectral arguments and weights.
- All-power Gaussian questions require an unbounded charge ledger or an independently proved finite reduction.
- Integer affine congruences constrain possible polynomials but do not enforce cycle incidence or full metagraph realizability.
- A divisor-height fold or a chosen branched cover must not be promoted to an arithmetic isomorphism.

The recurring positive method is to retain the observer, shared coordinate, integer lift or incidence relation that a tempting quotient erased. The resulting theorems are useful precisely because they state where the transfer is exact and where additional structure is still required.
