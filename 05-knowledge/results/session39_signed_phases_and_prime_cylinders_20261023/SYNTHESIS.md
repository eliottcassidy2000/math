# Session39: hidden attachment phases, topological data, and signed Collatz cylinders

## Overview

The useful common idea is **transport with a forgotten coordinate**, not an identification of tournaments, quantum Hall fluids and Collatz dynamics. This session gives actual maps and actual obstructions:

- unresolved cyclic-three tournament cells carry affine F3 connections;
- long attachment cycles hide their phase from bounded local motif observers;
- linking forms and quadratic refinements retain information lost by the dyadic group tower;
- signed Collatz valuation words define arithmetic cylinders, with exact prime and target-basin statements;
- a formal graph-to-Chern–Simons construction has explicit parity and time-reversal limitations.

No proof of Collatz convergence, a microscopic Hall phase, or an all-size tournament-metagraph classification is claimed. “Hidden” here means invisible to a specified mathematical observer; it is not a claim about astrophysical dark matter.

## 1. Incoming work and the actual new attachment carrier

Incoming `e41a2c82b34fea4dffc6c1ee5dab07bfb8c56ded` proves a [single cyclic-three-cell decoder](../incidence_wiring_20261023/MARKED-CELLS.md): retain common-owner selector overlaps, and one chirality bit when at least two owners occur. Its unresolved obligation is joint alignment of multiple cells.

Incoming `bc3c88ec033d809dee0d665314005308012f0bd6` adds [rooted-five port repairs and marked electrical guards](../metagraph_five_port_phase_frontier_20261023.md). In particular it already proves generic fixed-motif-cutoff failures in prime-cyclic families, obtains a native maximal linking squareclass, and distinguishes marked potentials from unmarked rank data. Those results have incoming priority and remain intact.

Our exact carrier applies to **every equitable tournament partition into singleton and triple cells**. Internal regularity forces each triple to be a directed C3. Each mixed 3-by-3 block has constant row and column sum. It is therefore uniform, a permutation matrix, or the complement of a permutation matrix. Normalizing a mixed block to a permutation gives

    x ↦ εx + δ,  ε in {+1,−1}, δ in F3.

Retain block polarity as well. Rotating the coordinates of each oriented triple changes δ by a gauge transformation but does not discard ε. Spanning-tree normalization leaves the jointly based cycle holonomies. Local chirality signs alone do not determine these phases.

There is an explicit topological representation: with ω=exp(2πi/3), the affine transport acts on the unit circle as z↦ω^δ z^ε. Positive ε is a rotation; negative ε includes complex conjugation. Composition reproduces the affine semidirect-product law, giving the dihedral subgroup of O(2). In the all-positive case, cycle translations become genuine U(1) holonomies. This preserves the connection data, but supplies neither a Chern–Simons action nor an anyon-braiding interpretation. A graph loop phase and a quasiparticle's framing twist are different observables.

This is not merely an artificial block model: the incoming order-31 two-mark cases with seven singleton cells and eight triple cells meet the theorem's hypotheses. Applying the decoder still requires their actual retained partition and alignment; it does not authorize naming unresolved vertices independently.

The general affine fixed-section principles already occur in THM-2622 and THM-2089. The new content is the authenticated tournament carrier and its constructive hostile families, not a claim to invent holonomy.

Proof and checks: [affine connections](marked-cells/AFFINE_CONNECTIONS.md).

## 2. A complete small family, then arbitrarily long hidden phases

For three intrinsically distinguished C3 cells with permutation cross-blocks, all 6³=216 presentations fall into exactly 16 tournament isomorphism classes under allowed cell rotations. Twelve have automorphism group C3 and four have trivial automorphism group. The cell scores 3,4,5 prevent hidden cell permutations.

Already with all matching signs positive, the three cycle phases give Hamiltonian-path counts

    2697, 2625, 2721.

They share the local chirality signs. Thus an explicit phase, not an additional count of those signs, is necessary.

For every g≥5, put g directed triples around a support cycle. Use matching blocks on the cycle, uniform forward blocks elsewhere, zero phase on an open support path, and phase δ on its closing edge. The three values δ=0,+1,−1 give three distinct strong self-converse tournaments Tδ on 3g vertices with the **same exact triangle-edge incidence signature E**.

The intrinsic cell scores are

    s0=3g−6;
    si=3g−2−3i for 1≤i≤g−2;
    s(g−1)=5.

They are pairwise distinct. Any isomorphism must preserve cells, and connected matching transport leaves exactly the global C3 rotations as automorphisms. Cycle phase distinguishes the three classes. The explicit converse map is

    (i,x) ↦ (g−1−i,−x).

It preserves δ after reversing the tournament, so the three phases are not merely converse aliases. Bidirectional arcs between matching-adjacent strong cells prove strong connectivity.

At g=5 the three Hamiltonian counts are

    45,132,777; 45,150,825; 45,102,765.

The exact E-blindness has a structural cause: girth at least five excludes both triangular support loops and competing length-two transport routes. Each local triangle-edge histogram sees at most one transported offset, which its local summation erases. The proof works for the declared signed matching-support class, not arbitrary unresolved tournament partitions of larger cell sizes.

## 3. Detection depth versus an actual mark

Every motif using fewer than g cells sees a forest of matching constraints; its phases can be removed by cell rotations. Hence all induced motif counts using fewer than g vertices agree across the three phases.

The two nonzero phases remain indistinguishable even on full g-vertex transversals: each uses one vertex per cell, so simultaneous coordinate inversion exchanges +1 and −1 without encountering an internal triangle arc.

With the intrinsic cells marked, the detection thresholds are exact:

- a g-vertex matching transversal distinguishes zero from nonzero phase;
- a (g+1)-vertex observation with two vertices in one cell and their directed internal arc distinguishes the two nonzero phases.

This boundary is stated for the named-cell observer. It must not be silently promoted to a minimal unrooted motif size in every instance.

An **actual distinguished vertex in a triple cell**, however, repairs this family. It splits its cyclic triple; a resolved triple separates every matching-adjacent cell. With intrinsic cell colors retained, directed refinement becomes discrete in at most 1+eccentricity of the marked cell rounds. This is a scoped constructive repair, not a universal one-mark theorem: the incoming singly marked doubly regular controls still fail.

The static singly-rooted strengthening is sharper: the two nonzero phases have identical complete rooted motif-count vectors at every corresponding vertex through cutoff g. Their intrinsic cell partition is equitable, so ordinary **unindividualized** refinement initialized by those vectors remains blind at every depth. This is not the actual-mark procedure above.

The inherited lift L(T)=C3[T,{a},{b}] adds two vertices while preserving strong connectivity, self-converseness and equality of the static rooted profiles. Its unique degree extremes recover the original tournament, so nonisomorphism survives. Applying it zero, one or two times covers every order n≥11. Writing n=3g+2r with r∈{0,1,2}, a separator in this observer format must have

    k(n)≥g+1=(n−2r)/3+1,
    liminf k(n)/n ≥ 1/3.

This is a necessary observer cutoff, not an algorithmic complexity lower bound or a theorem about arbitrary reconstruction data. The same lifts give three distinct strongly connected, self-converse exact-E collisions at every n≥17, without requiring unequal Hamiltonian counts at every lifted order.

Details: [rooted cutoff](marked-cells/ROOTED_CUTOFF_BOUNDARY.md), [all-order lift](marked-cells/ALL_ORDER_ROOTED_BOUNDARY.md).

Classical Stockmeyer–Kocay deck counterexamples already give much stronger unrooted-count failures, so no broad reconstruction-novelty claim is made. The additional distinction here is pointwise rooted-count equality combined with ordinary-refinement blindness and a one-vertex individualization repair. Counts are not indexed hypomorphy: even one specified three-vertex set differs between the phase twins. The [primary-literature reconciliation](RECONSTRUCTION_PRECEDENTS.md) separates these notions and records the relevant Lopez and Egrot–Hirsch results.

Proofs and complete finite receipts: [girth-phase hostiles](marked-cells/GIRTH_PHASE_HOSTILES.md). The original packet checks all 216 three-cell presentations, 7,776 two-cycle affine connections, all 7,776 signed C5 phase assignments, mixed block polarities, and the marked reconstruction controls.

## 4. What seven flux quanta and two electrons do, and do not, specify

The primary Nature Physics paper using this example employs **four vortices per electron** at filling 2/7. Its different parton description has two fermionic partons sharing three flux quanta. It does not assign three flux quanta to each physical electron. See [Yang et al., Methods and Extended Data Fig. 1](https://www.nature.com/articles/s41567-026-03338-9).

Three-flux attachment is nevertheless a coherent alternative conditional construction: odd flux attachment changes fermions into composite bosons. If those bosons occupy a bosonic integer-Hall sector with

    K_B = [[0,1],[1,0]],  t=(1,1),

then the flux-attached matrix is

    K = K_B + 3ttᵀ = [[3,4],[4,3]],
    tᵀK⁻¹t = 2/7,   det K = −7.

This does not prove that the specified microscopic system forms that phase. The finite pair of numbers N=2 and Nφ=7 does not supply geometry, shift, interaction or spin-vector data. In particular, mean-field kN flux bookkeeping and spherical k(N−1) Jastrow degree bookkeeping are different statements.

Odd K requires spin/local-fermion qualifications. A characteristic refinement is not automatically the same observable as a representative quasiparticle exchange phase. The report gives explicit conventions, the charge vector, and the relevant primary literature, including [Senthil–Levin's bosonic integer-Hall construction](https://arxiv.org/abs/1206.1604).

Source review and exact matrix controls: [flux and twists](flux/FLUX_TWISTS_AND_TOURNAMENTS.md).

## 5. Exact tournament-to-topological-data constructions and their limits

There are at least two different constructions, and they do not have the same output.

First, use a grounded native metagraph Laplacian A as a *formal* Abelian Chern–Simons matrix. Its characteristic covector diag(A) modulo two is exactly the retained self-converse-node indicator. Ordinary tiling fibres are odd, non-self-converse merging doubles them, and loops contribute even degree. In that declared basis, a local excitation Az has exchange parity

    (−1)^(sum of z_i over retained self-converse nodes).

This is a precise parity dictionary for the chosen formal construction, not an identification of graph vertices with physical electrons. See [native characteristic parity](flux-characteristic/NATIVE_CHARACTERISTIC_ADDENDUM.md).

Second, for an even-order tournament take its skew adjacency S and the doubled BF matrix

    K_S = [[0,S],[Sᵀ,0]].

Modulo two, S=I+J is invertible; its integer determinant is an odd square. Thus K_S is an even nonsingular integral matrix of signature zero, with odd fourth-power discriminant order. This formal bosonic construction cannot reproduce the native metagraph's dyadic sector. The two constructions are not interchangeable.

Nor is the graph critical group torsion in the ordinary first homology of a graph: that homology is free. The critical group is a cut-lattice discriminant. Adding two-cells or a physical realization needs additional specified data.

## 6. Beyond the complete dyadic tower

Incoming work and an independent Session39 reconstruction agree that the maximal native order-eight 2-primary linking summand is ⟨1/256⟩. The coordinate coefficients 105/256 and 153/256 are compatible: 9²·153≡105 modulo 256. The numerator is not an invariant without its generator.

Scaling all conductances by three preserves every labelled dyadic lift space and every binary divided-energy form, yet changes the highest linking unit class from 1 to 3 modulo eight. This is a weighted-network control, not a same-fibre native port counterfeit.

At order five, four of the seventeen switches invisible to the full dyadic tower are distinguished by the same declared rooted characteristic quadratic refinement. Thirteen remain invisible, including the earlier 41096/36680 tree-count hostile. Thus a quadratic sidecar can help, but does not complete reconstruction.

The native linking form has no antiisometry: every order-256 element has self-linking numerator 1 modulo eight, whereas sign reversal would require 7. Converse fixes the merged vertex divisors and hence their critical group, although it exchanges literal black twins. Thus tournament converse cannot implement formal K→−K in this construction. The [converse versus reversal obstruction](linking/CONVERSE_IS_NOT_TIME_REVERSAL.md) proves this precisely; it is not a claim about every possible physical construction from tournaments.

Loops are invisible to A alone, but their counts are recoverable from supplied native degrees. Even A plus these degree sidecars does not authenticate individual path-labelled partners inside a fixed aggregate edge cell.

Proofs, exact positive and hostile controls: [linking packet](linking/README.md).

## 7. The signed Collatz graph and its genuine prime connection

For positive odd n define

    Uε(n)=(3n+ε)/2^v2(3n+ε),  ε=+1 or −1.

The positive 3x−1 map is conjugate by sign to the negative-integer 3x+1 map. Its three known odd cycles are

    1;
    5→7→5;
    17→25→37→55→41→61→91→17.

No claim that these exhaust all positive starting behavior is made.

For a finite prescribed valuation word (a1,…,ar), put A=sum ai. Each signed map realizes that word on exactly one odd residue class modulo 2^(A+1). The two classes are negatives. Therefore every finite word occurs for infinitely many prime starts and infinitely many composite starts. For a fixed word, its asymptotic proportion among primes is 2^(−A).

This is a real arithmetic connection, but it also proves an observer limitation: no procedure using only finitely many valuations at each input can always decide primality. It is not a limitation on primality algorithms that read the integer itself, and no uniformity as the word length grows is supplied.

The classical signed cycle equation and existing repository prime-power work are credited in [signed cylinders and primes](negative-collatz/SIGNED_CYLINDERS_AND_PRIMES.md).

## 8. One extra step reaches any admissible target, with the clock retained

Let B be the usual carry of a length-r valuation word, so its affine map is (3^r n+εB)/2^A. For any positive odd target y coprime to three, choose a final exponent b satisfying

    2^(A+b)y ≡ ε(2^A+3B) modulo 3^(r+1).

The fact that two generates the units modulo powers of three gives exactly one class of b modulo 2·3^r. Sufficiently large members give positive odd starts

    n = [2^(A+b)y−ε(2^A+3B)]/3^(r+1)

following the prescribed prefix and then one more odd step to y. Exact valuation admissibility is proved, not assumed from a formal affine equation.

For 3x−1, the common prefix (1,2) occurs at starts 38837,6069,21, which reach 1,5,17 after final exponents 17,12,2 respectively. Every finite valuation cylinder therefore meets all three known cycle basins. Even with eventual entry into one of those basins promised, an always-finite valuation-only observer cannot decide which basin.

This inverse-orbit ubiquity has classical precedent: [Zarnowski, Theorem 3](https://www.fq.math.ca/Papers1/46_47-2/Zarnowski_1-09.pdf). The report claims no priority for that general density principle. Its explicit signed exponent clock and observer formulation are the useful connection here.

Do not intersect two unrelated infinitude statements: primes occur in every cylinder, and each known basin occurs in every cylinder, but this does not prove infinitely many primes in every basin-cylinder intersection.

There is also a quantitative clock boundary. Reaching a fixed target in at most r actual odd steps covers only polylogarithmically many starts below X for fixed r; even r=o(log X/log log X) covers o(X). This is a bound on native arrival depth, not proof-complexity or an obstruction to an unbounded-depth 153N argument. No missing uniform source budget has been supplied.

Details: [closed-target cylinders](negative-collatz/CLOSED_TARGET_CYLINDERS.md), [coverage clocks](negative-collatz/COVERAGE_CLOCK_BOUNDARY.md).

## 9. The next precise frontier

The attachment theorem now supplies a real joint decoder for the triple-cell cases left by marked refinement. The next task is to combine that decoder with the actual source/path/target flags and to handle unresolved cells of larger size, without flattening valid nonzero holonomy.

The electrical task remains to choose a sufficient marked phase bank for an entire feasible switch lattice, rather than one family. The signed arithmetic task remains an unbounded-depth coverage theorem with its actual clocks and conserved source budget. None is replaced by a shared number or by the word “topological.”

Independent audit packets: [core audit](audit/), [extension audit](audit-extensions/AUDIT_EXTENSIONS.md), [final observer and reversal audit](audit-final/AUDIT.md), [all-order lift audit](audit-lift/AUDIT.md). The [research map](RESEARCH_MAP.md) records the exact source/target and information-loss distinctions.
