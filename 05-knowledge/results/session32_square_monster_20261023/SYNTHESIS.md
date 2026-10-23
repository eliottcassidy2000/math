# Session32 synthesis: exact bridges, preserved ports and arithmetic clocks

**Final research synthesis.** Started2026-10-23 at00:48 UTC. Proofs, cited prior theorems, finite certificates and open connections are distinguished throughout. These are verified connections and repaired formulations, not a claim that all results are new to the literature.

## 1. The supplied integer and prime partition are exact group data

The integer is

N=39,046,339,332,000=|Monster|/(|Baby Monster|·29·41·59·71).

The three sets are precisely:

- {17,19,23,31,47}: equal nonzero prime valuations in both orders;
- {29,41,59,71}: present in Monster, absent from Baby Monster;
- {2,3,5,7,11,13}: present in both, with larger Monster valuations.

The same three sets are also the prime-support shells of Fi22, Baby Monster and Monster; see the [independently source-checked Fischer support note](FISCHER_SUPPORT_SHELLS.md). This does not import the richer pre-audit modular-filtration claims in a concurrent packet.

There is a stronger action interpretation. Monster acts on its2A involutions with stabilizer2.B. The first set consists of primes whose full Sylow subgroups fix a2A point; the third has some order-p subgroup fixing a point but no full Sylow subgroup doing so; the absent set has no order-p fixed subgroup at all. The actual conjugacy-class index uses2.B, so its small-prime core is N/2, not N.

The exponent arrows of N are exactly13→2→5→3→7→4 and11→1. Identifying13 and4 as graph states closes a five-cycle. This is not a well-defined modulo9 quotient of the whole exponent map:11 and2 coincide modulo9 but have different output residues. Even the restricted cycle cannot be represented by an integer polynomial modulo9. Replacing N by the centralizer-normalized N/2 changes the quotient cycle to a two-cycle. The normalization and lost arithmetic data therefore matter.

Sources, proof and independent audit: [Monster note](monster/MONSTER_BABY_PRIME_PARTITION.md), [fixed-point interpretation](SYLOW_FIXED_POINT_PARTITION.md), [valuation quotient](VALUATION_CYCLE_AND_QUOTIENT.md), [independent audit](monster-audit/AUDIT.md).

## 2. The four absent primes have a geometric edge defect

Aricheta's cited modular-curve theorem supplies a genuine centralizer/supersingular-rationality bridge, with good-level hypotheses retained. Exact finite-field calculations show that each of29,41,59,71 has one nonrational Frobenius pair of coarse points on X_0^+(2), despite rational underlying supersingular j-values.

The two points represent parallel unordered2-isogeny edges between the same rational vertex pair:

- p29: {2,25};
- p41: {3,32};
- p59: {28,48};
- p71: {41,66}.

Frobenius fixes the vertices and swaps those edges. Forgetting level/edge structure hides exactly the observed defect. The resulting full Atkin–Lehner quotients X_0*(2p) have genus1. The computations were checked by two producer routes and a separate independent finite-field enumeration. Coarse-orbit rationality is not automatically an F_p-defined oriented isogeny; p2 is excluded by the level condition.

This distinguishes the Monster-only primes from the Baby primes. It does not distinguish the two different nonzero valuation regimes by itself, nor supply a Hamiltonian-gap theorem.

## 3. A direct Monster–Paley action, not only an analogy

The Monster7A normalizer has a quotient by Held isomorphic to AGL(1,7). Its seven-point affine coset action is the full orientation/converse symmetry group of the Paley tournament P7. The subgroup7:3 of order21 preserves orientation, while the other coset reverses it.

This is an action of a quotient of a proper normalizer, not a seven-point action of the whole Monster. The coset action is essential: conjugation on the cyclic7-subgroup alone fixes the identity and does not supply translations. Likewise global converse is not the operation of complementing only the nonpath arcs used for metagraph wiring.

Exact permutation checks give H(P7)=189, |Aut(P7)|=21, nine fixed-path presentations and three blue presentations. Thus the verified order7 Paley center has deficit6. This explains a real occurrence of7 and21 in a shared group action and a native metagraph fiber. It does not explain why7 and21 themselves are absent as Hamiltonian counts.

The CFT source establishes a parafermion embedding but leaves some Held-action/lifting statements at the level of expectation or evidence; those are not silently promoted here. See the [Paley bridge](paley/PALEY_MONSTER_BRIDGE.md) and its [independent audit](paley-audit/AUDIT.md).

## 4. The square-sum arrival picture is completely explicit

For the graph G_n on1..n with distinct a,b joined when a+b is a square, the arriving vertex has degree

d(n)=floor√(2n−1)−floor√n.

The component counts at12,13,14 are3,2,1. Connectivity persists thereafter. The three components at12 are exactly the user's paths of sizes3,4,5. Arrivals13 and14 merge components;15 creates the first cycle. Component and cycle-rank updates distinguish these operations.

The degree-two births are globally, not merely through36:

13–15,19–32,36–40,49–50.

Degree-one births are3,5–12,16–18. Thus the phrase 'each node up to12 adds one edge' needs the exceptions1,2,4. Every fixed birth degree j has an explicit finite generating polynomial obtained by intersecting two square-root intervals. The first total counts by degree0..9 are3,12,24,34,47,59,68,84,92,105; the tempting12j rule already fails at3.

A further all-scale law explains silent square births: the downward floor jump cancels an upward jump exactly at the two Pell families m²−2k²=−1 or−2, with n=k². First examples are9,25,289,841,9801,28561. At25 the square49 channel arrives while square25 leaves the birth window, so degree2 persists even though its neighbor ports change. This arithmetic clock alone is not a Hamiltonian criterion, as n9 shows.

Details: [square-sum proof packet](square-sum/README.md), [Pell clocks and generating polynomials](PELL_ARRIVAL_CLOCKS.md), [independent graph audit](graph-audit/AUDIT.md), [Pell addendum](graph-audit/PELL_RESET_ADDENDUM.md).

## 5. Why15,23,25 are special, and why3–4–5 totals are insufficient

The exact path counts up to reversal are one each at15,16,17; three at23; ten at25; zero at the other orders2..25. Small cut certificates handle18–22; a corrected endpoint-aware forced-cycle argument handles24. At25 the ear11–25–24 is compulsory and repairs the obstruction. Degree2 alone does not force both edges in a path: G23 has a Hamiltonian path ending at the degree-two vertex22.

For the user's25-vertex blocks, internal edge counts are4,3,5, but the first block is branched. Contracting the three trees as a multigraph preserves cycle rank8 yet destroys Hamiltonian feasibility: deleting the three resulting supervertices leaves six core components, exceeding the allowed four. This is a direct obstruction in the proposed arrangement, not a hypothetical unrelated counterexample.

An exact port-state construction repairs the loss. The blocks have5,101,149 endpoint-labelled states. The endpoint budget leaves1,543 triples, which produce19 degree-feasible assemblies; nine contain disconnected cycles, leaving precisely ten Hamiltonian paths. Endpoint allocation and component pairing are needed, not only degree or cycle-rank totals. Endpoint identities may be anonymized only when the target does not prescribe endpoints and the required pairing/multiplicity data is retained.

At31 a separate forced-chain kernel has exactly two spanning2-factors, of component sizes11+20 and14+17, hence no Hamiltonian cycle. At32 there is exactly one spanning cycle up to reversal/rotation; removing32 gives the unique G31 Hamiltonian path with prescribed endpoints4 and17. These are finite exact threshold refinements. Gerbicz's all-n path existence from25 and cycle existence from32 are cited prior results; their full base/gluing tables were not reaudited here.

See [port transfer](ports/PORT_TRANSFER.md) and the [cycle-threshold proof and certificate](cycle-threshold/README.md).

## 6. Two further exact maps, with their losses identified

Path rooted-forest recurrences have transfer matrix [[3,−1],[1,0]], integrally conjugate to the cat map [[2,1],[1,1]]. Consequently the path blocks have a real common operator with the modular cat-map thread. But tr(A²)=7 and the rooted-forest count on P4 is21: they are different observables. The3/4/5 components of G12 have rooted-forest counts8,21,55, product9240. At25 the branched block has52 forests rather than the55 of a path with the same edge count.

Applying square-sum adjacency directly to the15 Monster primes produces a connected14-edge tree with five leaves, not the interval graph G15. Removing the four Monster-only primes leaves two Baby-support components; they restore a connecting path and a pendant path. Same5/4/6 block cardinalities do not preserve within-block adjacency. This finite arithmetic observation is not claimed to be forced by moonshine.

See [cat-map path forests](CAT_MAP_PATH_FORESTS.md) and [prime-support square graph](PRIME_BLOCK_CARDINALITY_HOSTILE.md).

## 7. Withheld research lane

The separate attachment-derived audit, proofs and test outputs are retained privately pending explicit publication approval. They are not included in this public packet.

## What remains open

The session supplies explicit maps and explicit obstructions, not a universal equivalence among all these subjects. It does not explain the forbidden Hamiltonian values7,21, prove the missing CM/MRB defect mechanism, settle Collatz, or classify complete tournament metagraph wiring. The strongest transferable principle is specific: preserve the port, orientation, level or collision coordinate that carries the target property before taking a quotient.

Specific further proof obligations are in [Next proof targets](NEXT_PROOF_TARGETS.md). They are not tasks still running at closeout.
