# Session 178: Boundary debt of an observed cycle

## Status and point

This packet gives a precise common construction for local arrows, observer holonomy, functional graphs, and arithmetic prime-factor accounting. The graph exact sequence is an elementary instance of relative graph homology, not a historical novelty claim. Applying it is useful because it separates three different claims: an observed cycle, an actual closed flow, and a trajectory that closes. None implies the next without additional conditions. No Collatz convergence or odd-perfect-number existence claim follows.

## 1. The universal observer defect

Let G be a finite connected directed multigraph; connected always means connected after forgetting directions. Loops and parallel edges are allowed. Let π:V→Y be a surjection. Form Q by identifying vertices according to π, **retaining every edge with its identity and direction**. In particular, parallel observed arrows must not be merged.

Use free abelian groups C1=Z^E and C0(G)=Z^V. The boundary is ∂e=e_target−e_source; let π_*:Z^V→Z^Y send each vertex basis vector to its observed basis vector. Then ∂Q=π_*∂. Write Z1(G)=ker ∂ and Z1(Q)=ker ∂Q. These are signed integral circulations, not necessarily nonnegative flows or individual directed walks.

Define the boundary debt of an observed circulation z by D(z)=∂z. Since ∂Qz=0, D(z) lies in K=ker π_*.

**Theorem 1 (exact observer sequence).** There is an exact sequence of free abelian groups

    0 → Z1(G) → Z1(Q) --D→ K → 0.

In particular, Z1(Q)/Z1(G) is canonically K, of rank |V|−|Y|. For each fibre choose an anchor a_y. The differences e_v−e_a_y, v≠a_y, form a basis of K. Thus precisely |V|−|Y| integral counters suffice to detect every failure of an observed circulation to be an actual circulation; fewer integral linear counters cannot do so.

**Proof.** The common edge group identifies Z1(G) with a subgroup of Z1(Q); D's kernel is exactly Z1(G). Every d∈K has total coefficient zero. In a connected graph, every total-zero vertex chain is a boundary: choose a root, join each vertex to it by an undirected path, and use signed oriented edges along the path to realize e_v−e_root. Therefore d=∂z for some integer edge chain z. Since π_*d=0, this z belongs to Z1(Q). This proves surjectivity. In each fibre, coefficient sum zero gives the displayed basis. An exact detector with r integral counters is an injective homomorphism K→Z^r. Tensoring with Q forces r≥rank K; anchor coordinates attain equality. ∎

The sequence splits as groups after choices of paths, but there is no asserted canonical splitting. If G is disconnected, replace the codomain by K∩im ∂; the connected theorem's full rank formula need not hold.

**Trajectory caveat.** D(z)=0 means that the signed edge chain is a circulation. It does not mean a specified ordering of its edges is a composable closed trajectory. Even a nonnegative circulation can split into several directed cycles. The theorem deliberately certifies chains, not order-sensitive walks.

## 2. Arithmetic counters and holonomy

Let h:V→A take values in an abelian group, and extend it linearly to H:Z^V→A. Label each edge with dh(e)=h(target)−h(source). Its accumulated label on an observed circulation is

    Hol_h(z) = Σ_e z_e dh(e) = H(D(z)).

This is an exact identity. Every true circulation has zero holonomy. All observed circulations have zero holonomy if and only if h is constant on every observer fibre: necessity follows from Theorem 1's surjectivity and the anchor basis; sufficiency follows by factoring h through Y. Thus a quantity that is a pure potential difference before observation can appear to have holonomy after incompatible states are identified.

If vertices are positive rational numbers, take h(n)=(v_p(n))_p, with finite-support prime valuations. The multiplicative counterpart is

    R(z) = ∏_e (target(e)/source(e))^(z_e),
    v_p(R(z)) = Hol_p(z).

These formulas are algebraic identities, not stochastic or information-theoretic analogies. Nevertheless, the prime counters can lose information about D: the linear extension of an injective state encoding need not be injective on vertex chains.

**Counterexample.** Let G be the path 2→4→8, and let the observer identify all three vertices. The two edges become distinct loops in Q. For z=edge(2,4)−edge(4,8),

    D(z)=−e_2+2e_4−e_8 ≠ 0,
    R(z)=(4/2)/(8/4)=1.

All prime-valued holonomies vanish, but the observed circulation is not a circulation in G. Both individual states and each arrow were retained correctly; multiplicative compression still erased a second-difference defect. Two universal counters are necessary here, while the sole nonzero prime coordinate supplies only one.

## 3. Collatz example, with no global conjecture

Use the ordinary Collatz map T(n)=n/2 for even n, and 3n+1 for odd n. On the genuine cycle 1→4→2→1 take the parity observer. The quotient has two vertices and three separately retained edges. Its signed cycle rank is 2, whereas the original graph has cycle rank 1; K is generated by e_2−e_4.

The observed even→even loop coming from 4→2 has debt e_2−e_4 and multiplicative ratio 1/2. It is visibly a loop to the parity observer but is not an actual closed orbit. The full three-edge Collatz cycle has zero debt and ratio 4·(1/2)·(1/2)=1. The theorem explains precisely which extra apparent cycle parity identification has introduced. It makes no statement about trajectories outside these three vertices.

More generally, take any finite set of actual Collatz arrows, including both endpoints. If the resulting graph is connected, Theorem 1 applies. A numerical window must retain outgoing endpoints or explicitly record boundary edges; discarding arrows is a different operation from observing vertices.

## 4. Exact connection to sigma-factor flow

Let N=∏_{p∈S}p^{a_p}, with distinct primes and positive integers a_p. Define the finite prime reaction vector for each component by

    s_p = v(σ(p^{a_p})) − a_p e_p.

Multiplicativity of σ on coprime factors proves the exact accounting identity

    Σ_{p∈S}s_p = v(σ(N)/N).

Equivalently, form a directed multigraph whose individual edges are the rational transitions p^{a_p}→σ(p^{a_p}). The prime label on each such edge is precisely s_p. Different components may be disconnected; Theorem 1's connected rank assertion is not invoked for this graph. The additive potential identity still applies without connectedness.

For a candidate equation σ(N)=kN, k a positive integer, introduce a balancing edge k→1. The product of all component ratios and the balancing ratio is 1 exactly when the candidate equation holds. Observe all vertices as one state, so every edge is an observed loop. Then the candidate equation is **zero arithmetic holonomy**, not universal zero boundary debt and not a directed prime-flow cycle. The distinction is essential, even for genuine perfect numbers.

For N=6 the edges are 2→3, 3→4, and balancing edge 2→1. Their ratios multiply to (3/2)(4/3)(1/2)=1, but their total boundary is e_1+e_4−2e_2, which is nonzero. Arithmetic zero holonomy therefore cannot be substituted for literal graph closure. This is a constructive bridge, not an assertion that perfect numbers are Collatz cycles.

Let F_{pq}=v_q(σ(p^{a_p})). The same identity in coordinates is

    Σ_p F_{pq} − 1_{q∈S}a_q = v_q(k), for every prime q.

This recovers exact sigma-factor supply/demand, including exterior primes. Projecting to q∈S is an observer operation on the accounting group and can hide exterior debt. For fixed exact components, if every inside coordinate is balanced, the residual

    E = σ(N)/(kN)

has zero valuations on S but can still have nonzero exterior valuations (and for arbitrary k can be rational). Consequently inside balance alone does not certify perfection. If k has no prime factors outside S, then E is an integer supported outside S, and full equality holds precisely when E=1. The common odd-perfect setup k=2 has 2 outside S, so that special integer conclusion cannot be invoked until the 2-adic balance is checked as well.

This is why local tournament realization in Session 177 is insufficient: orientations on selected primes discard multiplicity, exponents, absent exterior vertices, and coordinate balance. Our construction supplies an exact accounting language, not a new arithmetic obstruction to all exponent choices.

## 5. What the unification does and does not buy

- **Local arrow detection:** accurate local arrows can become spurious cycles solely through vertex identification; edge identities matter.
- **Observer information:** |V|−|Y| is the exact number of independent universal additive debt counters for connected finite graphs. This is an algebraic rank, not an entropy claim.
- **Holonomy:** observed closed chains can carry potential differences because their original endpoints do not cancel.
- **Collatz:** parity exposes an explicit hidden endpoint debt in a concrete finite functional graph.
- **Sigma networks:** prime valuation labels are exact reaction balances; their vanishing is weaker than universal state-chain closure.
- **Tournament interfaces:** any finite tournament is a directed graph to which the theorem applies under vertex observations. No transfer of tournament Hamiltonicity to arithmetic trajectories is claimed.

The next useful mathematical problem is not to assert a single universal dynamical model. It is to determine which reduced counters remain complete on a specified admissible class of arithmetic networks. The 2→4→8 example shows why completeness needs proof.

## Reproduction

Run `python3 verify.py` and `python3 -O verify.py` in this directory. Checks use explicit exceptions, exact integer/rational arithmetic, all connected simple directed orientations through four vertices with all set partitions, and concrete Collatz/sigma examples. Finite checks corroborate the proof; they are not its justification.
