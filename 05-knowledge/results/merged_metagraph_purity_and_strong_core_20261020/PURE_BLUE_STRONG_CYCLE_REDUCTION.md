# Pure-blue strong atoms: odd circulants and the even-rigidity corollary

Status: PROVED elementary reduction, using CITED classical Hamiltonicity and Xiaoyun Lu's published crossing theorem. The unrestricted pure-blue atom classification remains OPEN. This note uses the user's term 'points of symmetry' for pure-blue merged nodes.

## Definitions and cited input

A tournament T is pure-blue when every normalization of every directed Hamiltonian path has a reflection-fixed fixed-path tiling. Equivalently, for every Hamiltonian word P=(p_0,...,p_(n-1)), its positional reversal permutation p_i -> p_(n-1-i) is an anti-automorphism of T. This is THM-644(a) plus the exact Hamiltonian-path fibre inversion THM-781.

Classical input: every strong tournament has a directed Hamiltonian cycle (Camion; also implied by Moon's vertex-pancyclicity imported in 05-knowledge/reference/CORE-PAPERS-TOURNAMENT-STRONG-EARS-2026-08-25.md).

Published crossing input: Lu's Theorem7, in 'On king–serf pair in tournaments', Discrete Mathematics 343(7) (2020), 111903, DOI https://doi.org/10.1016/j.disc.2020.111903. Its precise statement and primary-source access scope are recorded in LU_CROSSING_SOURCE.md. The needed consequence is that no strong tournament of order>=6 is Hamiltonian-crossing for any bipartition.

## 1. Cycle-cut reflection lemma

Let T be pure-blue and strong, n>=3. Choose a directed Hamiltonian cycle indexed by Z/nZ, with arcs v_j -> v_(j+1). Cutting immediately before v_i gives Hamiltonian path

P_i=(v_i,v_(i+1),...,v_(i-1)).

Purity says that the permutation

r_i(v_j)=v_(2i-1-j)

is an anti-automorphism. Hence the composite r_(i+1) r_i is an automorphism, and direct calculation gives

r_(i+1)r_i(v_j)=v_(j+2).

Thus rotation by2 along EVERY Hamiltonian cycle is an automorphism.

Tournament automorphism groups have odd order: an even-order element would power to an involution, whose transposed vertex pair would reverse its arc. Therefore, if 4 divides n, rotation2 has forbidden even order n/2. No pure-blue strong tournament has order divisible4.

If n is odd, rotation2 generates rotation1. The Hamiltonian cycle labels therefore identify T with a circulant tournament on Z/nZ. In particular T is vertex-transitive and regular, with every score (n-1)/2. This implication is necessary, not sufficient: regular circulant examples at n7 are mixed.

## 2. Excluding the remaining even orders

Suppose n=2k with k odd, n>=6. The rotation2 automorphism has exactly two vertex orbits:

U={v_0,v_2,...}, W={v_1,v_3,...}.

Scores are constant a on U and b on W. The anti-automorphism r_0 exchanges U,W, so a+b=n-1. Since the right side is odd, a!=b.

Every Hamiltonian path P has an anti-automorphic positional reversal, so its endpoint scores sum to n-1. They cannot both lie in U (sum2a), or both in W (sum2b). Thus T is U-Hamiltonian crossing. Lu's Theorem7 contradicts strong connectivity at order>=6.

Orders divisible4 were already excluded. No order-two tournament is strong. Therefore:

**Every nontrivial pure-blue strong tournament has odd order and admits a regular circulant presentation along every Hamiltonian cycle.**

No large census is used. The two ingredients test genuinely different loss coordinates: all cycle cuts give the automorphism, while all Hamiltonian paths exclude the two-orbit survivor.

## 3. Every even pure-blue node is rigid

Let tc(T)=H(T)/|Aut(T)|. THM-2454's center law follows from the canonical ordered strong-component decomposition: a self-converse stack has converse-paired outer components, and possibly one self-converse central component. The grid reflection exchanges paired component fibres. All independent path choices can be reflection-fixed exactly when every outer component fibre is a singleton and the central component, if present, is pure-blue.

If T has even order, an unpaired central component would have even order, since the outer pairs contribute an even number of vertices. Section2 excludes a pure-blue strong central component of even order. Hence there is no center. Every paired component is rigid, so stack multiplicativity gives tc(T)=1.

Conversely, every self-converse rigid class is pure-blue by the odd-fibre involution argument (indeed its fibre is a singleton). Therefore, for every even n,

**T is pure-blue iff T is self-converse and H(T)=|Aut(T)|.**

This does NOT hold at odd order: cyclic regular T5 has H=15, |Aut|=5 and three blue tilings, so is pure but nonrigid.

Consequences:

* Every even-order pure-blue node has exactly one tiling. For even n>=4 it has exactly one incident complement line: these points of symmetry really are blue leaves. At n2 the antipode is not free and this line convention does not apply.
* For even n>=4 that edge cannot be a loop, since a loop requires two tilings in the same fibre.
* Every even pure-blue leaf at n>=4 is attached to a mixed node. To see this, its tournament cannot be strong. Along any Hamiltonian path of a nonstrong tournament, the first and last vertices belong to distinct extremal strong components, so the first beats the last. That arc is nonpath for n>=3. Flipping all nonpath arcs reverses it and therefore closes the unchanged Hamiltonian path into a directed Hamiltonian cycle. The complement partner is strong. At even order it cannot be pure-blue, and it is blue-touched, so it is mixed.
* The old unrestricted 'all pure-blue points are leaves' claim remains false at n5 and n7. The repaired theorem is universal at even order; at odd order singleton fibre suffices, but pure does not.
* Conditional on the inherited PROVED-through62 strong-rigidity classification in THM-2477 and Narayana-stack count THM-2453, whose lengthy proofs were not reaudited here, this yields pure-blue(n)=c(n/2) for every even n<=62, where c(0)=c(1)=c(2)=1 and c(j)=c(j-1)+c(j-3). This packet does not extend that inherited finite-range rigid-atom classification; the all-even purity/rigidity equivalence itself is unconditional under the cited crossing theorem.

## 4. The nonstrong periphery is independent

For every n>=3, every complement-line partner of a nonstrong tournament is strong. Indeed, a directed Hamiltonian path visits the canonically ordered strong components in order. Its first vertex dominates its last vertex because these lie in distinct extremal components. Flipping all nonpath arcs reverses that endpoint arc while preserving every path arc, closing a directed Hamiltonian cycle. Thus every line touching the nonstrong locus immediately enters the strong locus. In particular the nonstrong merged nodes form an independent set, have no loops, and admit no within-periphery complement lines.

This is an explicit metagraph corollary of the existing THM-333-apex-tile-strong-connectivity mechanism, not a new proof of the apex-closing fact. That older file uses the opposite bit convention and uses SC for strongly connected; here SC elsewhere means self-converse, so the words strong/nonstrong are written out. The conclusion is independent of blue/black colour, purity, parity, and all node counts. The converse is false: both endpoints of a complement line can be strong.

## 5. Hamiltonian-cycle automorphism injection

For odd pure-blue strong T, map each directed Hamiltonian cycle, counted up to cyclic rotation, to its successor permutation on V(T). Section1 proves that this n-cycle is an automorphism. The successor permutation determines the entire directed cycle, so the map is injective. Therefore

HC(T) <= #{g in Aut(T):g is one n-cycle} <= |Aut(T)|-1.

The -1 uses n>=3, so the identity is not an n-cycle. This does not assert that every n-cycle automorphism traces directed arcs of T. Controls: C3 has one Hamiltonian cycle and |Aut|=3; cyclic regular T5 has two Hamiltonian cycles and |Aut|=5.

Potential next obligation: combine a valid lower bound for HC in regular tournaments with an upper bound on odd-order automorphism groups to bound possible new pure-blue strong atoms, then attack the residual circulant cases by a single non-anti-reversible Hamiltonian path. Keep the exact cycle-count convention and the reference's constants; no such all-order exclusion is claimed here.

## Scope and failure boundary

This note reduces, but does not solve, the odd strong alphabet question {C3,T5}. It does not identify purity with regularity, circulancy, transitivity of automorphisms, or tc=1. It does not replace all Hamiltonian paths by only cuts of one Hamiltonian cycle. It also does not reprove or extend the all-order rigid classification problem left open beyond62 by the cited THM-2477 snapshot.
