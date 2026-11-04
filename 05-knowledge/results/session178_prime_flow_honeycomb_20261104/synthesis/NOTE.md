# Session 178: exact arrows, geometric wires, and arithmetic balance

2026-11-04. **Status:** synthesis of the sibling constructive proofs and their scoped computational checks, with inherited results explicitly credited. No claim of historical priority, Collatz convergence, or odd-perfect-number nonexistence is made. The new work provides usable representations and a support-specific arithmetic decision procedure; these are not equivalent to settling the global conjectures.

## 1. Three ways to make arithmetic visible

The session answers the geometric question constructively in three different senses. Throughout, C(n)=n/2 for even n and C(n)=3n+1 for odd n. The degree-six triangular lattice of hexagon centers is distinguished from the degree-three honeycomb graph.

**All sites, exact time, long arrows.** The [detector packet](../detector/detector.md) specifies a bijection between positive integers and triangular-lattice sites, shell by shell. Shell r has 6r sites; the labels through that shell end at 1+3r(r+1). Decode endpoint labels a,b and test both possible orientations using the exact Collatz equations. This labels every site and preserves every ordinary step as an arrow, but most arrows are nonlocal jumps. “Space-filling” means a site bijection, not a single Collatz orbit visiting all integers, a continuous space-filling curve, or a nearest-neighbor arithmetic embedding.

For this precise spiral convention, exactly eight Collatz arrows connect neighboring sites:

    1→4, 2→1, 3→10, 28→14, 30→15, 32→16, 34→17, 36→18.

The completeness proof bounds the larger endpoint's shell by six and then exhausts that finite region. It is not a conjecture inferred from a large sample, nor a statement about every possible lattice labeling. The cycle 1→4→2→1 is retained arithmetically, although its edge 4→2 is not neighboring in this spiral.

**All integers, local halving, long odd arrows.** The [honeycomb packet](../honeycomb/RESULTS.md) uses the unique decomposition n=2^a q with q odd and assigns Phi(n)=(a,q−1). This is a bijection onto the half-plane subset x≥0, y≥0 even, not onto all honeycomb sites. Every even step is one horizontal honeycomb edge. At a boundary point (0,y), the odd arrow lands at (a,z) determined by

    a≥1, z≥0 even, 2^a(z+1)=3y+4.

This gives an exact global arrow chart. Following the odd jump and its subsequent horizontal halvings returns to x=0 and realizes the accelerated odd map. The contracting y≡0 mod4 and expanding y≡2 mod4 cases are exact formulas, not a proof of average contraction. A valuation computation or equivalent arithmetic equation remains part of the detector.

**Finite trees, local wires, variable time.** The same packet gives explicit coordinates for every finite binary tree of height d inside the honeycomb after subdividing edges. The complete tree uses

    (d+2)2^(d+1)−4

inserted wire vertices, in addition to its original branch vertices. Wire lengths at a given level are D+4 and D+2, where D=2^(d−t). Edge interiors are disjoint; extra lattice edges are omitted. Deleting unused branches embeds any finite inverse tree rooted at a nonperiodic Collatz vertex, including 5. Arithmetic labels belong to the marked branch vertices; routing sites are auxiliary. One Collatz step becomes one entire directed wire, not one lattice move.

These finite layouts depend on d and are not nested. They therefore do not supply a global infinite layout by taking their union. This limitation is compatible with the separate inherited global existence theorem below.

## 2. What the arrow oracle costs

The detector packet supplies a fixed finite alphabet for a binary arithmetic strip. For least-significant-first bits x_i,y_i, each cell checks

    A x_i+c_i=y_i+2c_(i+1), c_0=B, c_L=0.

Use (A,B)=(2,0), with carries 0,1, for doubling, or (3,1), with carries 0,1,2, for the odd branch. The latter additionally requires x_0=1. Neighboring cells must agree on branch tags and adjacent carries; the multiplier cannot change midway. Summing the cell equations proves the whole affine equality.

This is a constant-state local verifier using O(log max(n,m)) cells and sequential reading time once the endpoint digits are supplied. It is not a constant-time verdict at one geometric site. Endpoint bits must be tied to the actual selected endpoint labels, and copying them across a drawing costs wires or additional certificate data. The strip certifies the integers it encodes; arbitrary unverified labels cannot be silently substituted for those endpoints.

A fixed residue pair cannot replace these data: for every M, the true doubling pair (n,2n) and false doubling pair (n,2n+M) have identical residues modulo M. Conversely, exact arithmetic from unbounded coordinate labels is allowed and is precisely what the site oracle uses. Unbounded valuations can be emitted or counted along a wire; they cannot all be stored in one fixed finite internal state.

## 3. Inherited geometric results and the actual additions

[Session 140's geometry proof](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/collatz_hex/PROOF.md) already established a global faithful directed honeycomb subdivision, using a cited infinite-grid minor theorem. Original vertices and edge paths are retained, including any unknown components. Its route lengths must be unbounded. That existence theorem is not replaced by the present explicit finite construction, and neither implies a canonical computable global layout.

The same source proves a uniform five-ordinary-step binary predecessor supply for every root prime to three, giving at least 2^r distinct predecessors after 5r steps. Consequently polynomial-growth lattices cannot support an injective uniformly bounded-step representation of an entire Collatz component. [Session 140's principles](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session140_collatz_hexagonal_geometry_20261102/principles/PRINCIPLES.md) already supplied triangular-site conjugacy and exact first-return suspension. Those broad mechanisms are inherited. The present additions are the fully specified spiral with its complete eight-edge list, explicit carry strips and their interface costs, the dyadic-ray chart, and concrete finite honeycomb wire coordinates.

[Session 167's honeycomb note](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session167_recurring_structures_20261103/honeycomb/NOTE.md) already counts 704 vertices in the depth-21 inverse ball at 5, against 694 honeycomb sites at radius 21. The difference is 10. The present honeycomb check reproduces it; it is not new. For a depth-d tree routed with maximum edge length L, the capacity condition becomes

    |B_tree(d)|≤1+3(Ld)(Ld+1)/2.

Thus L=1 fails for this finite example; the count alone does not rule out L=2 there. The all-depth predecessor growth supplies the stronger unbounded-stretch conclusion.

Session 167 also already gave the common colored honeycomb carrier for three arithmetic reflections and a pinned three-flip tournament cube. With target differences p=b−a, q=c−a, g=gcd(p,q), the arithmetic deck subgroup is Z(q/g,−p/g), whereas the cube kernel is 2Z². Primitive arithmetic winding closes the cube only after two laps. Those exact return laws are distinct from a Collatz dynamical conjugacy and remain useful inherited test cases for observer information.

## 4. An exact theorem for what observation forgets

The [unification packet](../unification/NOTE.md) considers a finite weakly connected directed multigraph G and a surjective vertex observation pi:V→Y. The quotient Q retains every individual edge and its direction, even if parallel arrows or loops appear. For signed integral edge chains, its theorem is

    0→Z1(G)→Z1(Q)→ker(pi_*)→0,

where the last map is the original boundary D(z)=partial z. Therefore |V|−|Y| independent integer-linear counters are necessary and sufficient to distinguish every observed circulation from an actual circulation. This is a universal additive-counter rank, not an entropy bound or a lower bound for arbitrary nonlinear encodings. For disconnected G, the terminal group is ker(pi_*) intersect im(partial), so the full connected rank formula must not be imported.

Potential labels h give accumulated holonomy H(D(z)). A true circulation has zero holonomy, but a compressed family of potentials may fail to detect nonzero boundary. For the path 2→4→8, subtracting its two edges has boundary −e_2+2e_4−e_8, yet multiplicative ratio (4/2)/(8/4)=1. All prime-valued potential counters vanish on this nonzero defect.

Likewise parity makes the genuine Collatz edge 4→2 appear to be a loop, with hidden boundary e_2−e_4. The full cycle 1→4→2→1 has zero boundary. These are finite exact examples, not statements about every orbit. Even zero universal boundary certifies only a signed circulation, not a prescribed ordered trajectory. A nonnegative circulation may decompose into several cycles.

The useful synthesis is operational: specify whether the consumer needs exact arrows, marked first returns, signed-flow closure, ordered trajectories, or arithmetic balance before choosing which information to discard.

## 5. Prime flow now has a finite exponent interface

[Session 177](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session177_odd_perfect_20261104/synthesis/NOTE.md) supplied Euler-aware allowed order graphs, LCM neighborhood closure, strong connectivity of actual perfect-number support, and exact cut accounting. It also proved broad induced tournament realizability among selected primes. Such local graph realizability leaves exterior factors and multiplicities uncontrolled; it is not perfection.

The [current prime-flow proof](../primeflow/PROOF.md) makes the exponent branches finite once a support S and Euler prime q are fixed. If d=ord_r(p)>1, a legal edge contributes

    v_r(sigma(p^a))=v_r(p^d−1)+v_r(a+1)

when d divides a+1. At d=1 the contribution is v_r(a+1), with divisibility by r required for a nonzero edge. This denominator exception is retained. Euler-compatible exponent parities filter permitted orders.

Writing s=|S| and C for the maximum sum of the constant LTE terms entering one target, any B satisfying

    C+(s−1)(floor(log_3(B+1))+1)<B

bounds every exponent strictly below B. Column capacities, role rounding, LCM conjunctions, and connectivity then shrink the domains soundly. On S={3,5,7,11,13}, both candidate Euler graphs initially pass connectivity and the abundancy ceiling. The q=5 box B=16 forces a_13≤1, incompatible with its positive even role. The q=13 box B=32 reduces the bounds to (4,2,4,4,1), after which source 5 has no outgoing edge. These are independently reproduced support certificates, not improvements on known general prime-count bounds.

For a proper prime-power block A_T, with outgoing and incoming factors O_T,I_T, exact cuts give

    I_T<O_T<2I_T if q is outside T,
    O_T<I_T<2O_T if q is inside T.

They retain more arithmetic information than connectivity. Bounds valid on one common branch can contradict these strict intervals; independently achievable extremes need not be simultaneously achievable, so they are safe bounds rather than exact feasible descriptions.

Enumerating the finite role-compatible exponent domains, rejecting exterior odd factors, and imposing every odd valuation balance gives an exact finite decision procedure for that fixed support and Euler prime. The Euler condition supplies exactly one factor of 2. The proof gives no uniform bound on the support size or largest prime, so it does not finish the global existence problem.

The unification packet explains an important distinction: sigma-component reactions and a balancing edge 2→1 encode perfection as zero prime-valuation holonomy, not as literal state-boundary closure or a Collatz cycle. Even N=6 has vanishing arithmetic holonomy but nonzero universal vertex debt in that representation. This does not weaken the exact arithmetic equation; it identifies which consumer the valuation counters actually certify.

## 6. Constructive next targets

These are proposals, not completed results of the session.

1. **Endpoint-bound arithmetic wires.** Combine the explicit digit strips with selected finite honeycomb routes and prove a complete certificate including bit copying, branch markers, and wire lengths. Acceptance should certify the exact chosen arithmetic endpoints, not merely some encoded pair.
2. **Adaptive inverse-tree layout.** Exploit actual Collatz branch sparsity to reduce the finite wire cost below the complete-binary template while retaining disjoint routes and ordinary-step markers. Separately, an effective nested global construction would need a compatibility theorem; the current layouts cannot simply be united.
3. **Finite-support arithmetic solver.** Turn the proved exponent box into a transparent row-table solver retaining all exterior factors and exact balances. Measure the incremental effect of LCM, capacity, and cut filters on a declared support universe. The universal fixed-support proof and an implemented benchmark are different deliverables.
4. **Restricted complete observer counters.** Characterize when prime valuations or another smaller potential family are injective on the boundary debts actually generated by a specified arithmetic language. The 2→4→8 defect is a mandatory control; universal state injectivity alone is insufficient.
5. **Fourth-target square transport.** Extend the inherited three-reflection honeycomb carrier to a finite square-sum graph with additional targets and positivity guards. The 15-vertex path's extra chord 1–3, of target sum 4, is a concrete first test. Preserve winding and layer information rather than treating the old cylinder as the entire graph.

The proofs and reproducible computational scopes are in the sibling packets. This synthesis introduces no additional experiment or unproved global implication.
