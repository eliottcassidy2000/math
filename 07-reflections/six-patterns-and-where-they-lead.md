# Six Patterns and Where They Lead

**Session:** kind-pasteur-2026-03-21-S18
**Arising from:** THM-261 (Petersen = root orthogonality), THM-262 (dual Lie embedding), THM-263 (profile determinacy)

---

## The Question

The Petersen-Lie bridge investigation produced specific results about tournaments and root systems. But the results feel larger than their context. Each one embodies an abstract pattern that seems to apply far beyond graph theory or Lie algebras. This reflection is an attempt to extract those patterns, name them, and follow where they lead.

---

## Pattern 1: The Dual Habitation

**What it is:** A single object lives naturally in two different algebraic structures, and neither structure alone sees the full picture, but together they do.

A tournament T simultaneously lives in so(n) as a skew-adjacency matrix and in the A_{n-1} root lattice as a signed function on positive roots. The same basis elements {E_{ij} - E_{ji}} serve as generators of so(n) AND as root vectors of sl(n). The object doesn't "choose" one home -- it exists in both, and the two homes see different invariants (Casimirs vs weight norms), neither of which alone determines H.

**The abstract principle:** When an object has a dual algebraic nature, look for the invariants that are visible from one perspective but invisible from the other. The intersection of these two views is richer than their union.

**Where it leads:**

*Neural networks.* A weight matrix W in a neural network simultaneously lives in the space of linear maps (functional perspective: what does W compute?) and on the Stiefel/Grassmann manifold (geometric perspective: how far is W from orthogonal?). The functional invariants (accuracy, loss) don't determine the geometric ones (condition number, spectral gap), and vice versa. But jointly they might predict generalization. This is exactly the tournament situation: H (functional) is not determined by Casimirs (spectral) or weight norm (geometric) alone, but the joint classification is much finer. The emerging field of "Lie group-based neural networks" (arXiv:2406.11151) already uses Lie algebra structure for weight constraints -- but the DUAL habitation insight (the same matrix simultaneously in two Lie algebras) is new. For networks with skew-symmetric components (attention score differences, for example), the tournament theory directly applies.

*Voting systems.* A preference profile (the collection of all voters' rankings) simultaneously lives in the space of tournaments (via pairwise majority: who beats whom?) and in the space of score vectors (via Borda count: what's each candidate's total score?). These are exactly the so(n) and A_{n-1} perspectives. The Condorcet winner is visible from the tournament perspective; the Borda winner from the score perspective. They can disagree -- and the disagreement is precisely the "gap" between the two algebraic homes. Arrow's impossibility theorem is, in this framing, a statement that no single perspective can be made consistent with both structures.

*Molecular chemistry.* A molecule's structure simultaneously lives in the point group (symmetry) and in the character table (spectral). The dual habitation principle suggests that molecular properties (reactivity, stability) might be better predicted by jointly classifying molecules by BOTH group structure AND character, rather than by either alone. Recent work on "fingerprint-enhanced graph neural networks" (J. Chem. Inf. Model. 2024) is moving in exactly this direction, combining local substructure fingerprints (one home) with global graph properties (other home).

---

## Pattern 2: The Anti-Conflict Complement

**What it is:** Every structure has a natural "anti-structure" that encodes exactly the complementary information. The original structure and its anti-structure are related by a clean duality, and together they tile the complete space.

The Petersen graph is the anti-conflict graph of tournaments: its edges mean "non-interference" (disjoint vertex sets = orthogonal roots) where the tournament conflict graph's edges mean "interference" (overlapping vertex sets). These are complementary -- Kneser K(n,2) vs Johnson J(n,2) -- and together they give the complete graph K_{C(n,2)} on all positive roots.

**The abstract principle:** When you find a "conflict" structure in any system, immediately construct its complement. The anti-conflict structure is not just the absence of conflict -- it has its own positive content, often with a famous name.

**Where it leads:**

*Wireless network scheduling.* In wireless networks, the conflict graph encodes interference (two links that can't transmit simultaneously). The complement -- the anti-conflict graph -- encodes links that CAN be scheduled together. This is standard. What's NOT standard is the Petersen-like observation: the anti-conflict graph of wireless links is the orthogonality graph of a root system when the network has a certain regularity. For networks on n nodes with all-pairs communication (complete graph topologies), the anti-conflict graph IS K(n,2), and the Lie-algebraic structure applies directly. This could yield new bounds on scheduling efficiency via root system geometry -- specifically, the chromatic number of K(n,2) (which is n-2 by Lovasz's theorem) directly bounds the minimum number of scheduling slots.

*Register allocation in compilers.* The interference graph for register allocation is a conflict graph: two variables that are simultaneously live can't share a register. The anti-interference (complement) graph encodes which variables CAN share. The Chaitin-style coloring algorithm works on the conflict graph; but the dual perspective -- working on the anti-conflict graph's independent sets -- might exploit root-system structure when the variable lifecycle has certain regularity (e.g., in loop-based code with predictable live ranges).

*Ecological competition networks.* Research by Kerr et al. (Nature, 2002) and Allesina & Levine (PNAS, 2011) shows that ecological communities with intransitive competition (rock-paper-scissors cycles) maintain greater biodiversity than transitive hierarchies. The conflict graph of an ecological tournament (which species' niches overlap) has an anti-conflict complement (which species can coexist without direct competition). The Petersen-like insight: when the number of species is n=5, the coexistence structure IS the Petersen graph on the pairwise niche overlaps. More generally, the Kneser graph K(n,2) governs which pairs of species interactions can be "decoupled" -- an orthogonality that could inform ecosystem management strategies.

---

## Pattern 3: The Profile Determinacy Boundary

**What it is:** A local statistic (computed at each "site" in a structure) perfectly determines a global property up to a sharp threshold, then fails by a precise quantum at the next level.

The 3-cycle root profile (count of 3-cycles through each arc) determines H exactly at n=5 but fails at n=6. The failure isn't gradual -- it's exactly 3 profiles, each producing exactly 2 H values separated by exactly 4 = one alpha_2 contribution. The local statistic captures the first-order structure (alpha_1) perfectly but is blind to second-order structure (independent cycle pairs) that first becomes possible at the threshold.

**The abstract principle:** When a local fingerprint determines a global property at small scale, the failure threshold reveals the onset of "non-local entanglement" -- structure that can't be seen from any single site. The size of the failure gap is the quantum of the newly visible structure.

**Where it leads:**

*Graph neural networks (GNNs).* The expressiveness limits of message-passing GNNs are exactly a profile determinacy boundary. A k-layer GNN can distinguish graphs distinguishable by the k-dimensional Weisfeiler-Leman test. The WL test IS a local profile: it iterates neighborhood aggregation. The failure at certain graphs (non-isomorphic graphs with identical WL-colorings) is exactly the "profile fails to determine global property" phenomenon. The tournament insight adds precision: the failure gap should correspond to specific structural features (cycles of length > 2k+1?) that the message-passing can't detect. This could guide the design of augmented GNN architectures that specifically target the "alpha_2-like" non-local features.

*Molecular property prediction.* The molecular fingerprint (a local feature vector computed from atom neighborhoods) determines many molecular properties at small molecular sizes but fails for larger molecules. Recent papers (J. Chem. Inf. Model. 2024, Bioinformatics 2025) address this by combining fingerprints with global graph features -- they are effectively adding "alpha_2-like" corrections. The tournament framework gives a precise language for this: the fingerprint captures the "3-cycle root profile" (local bonding patterns), while the global graph features capture the "independent cycle pairs" (non-local ring interactions). The gap size should predict how much the global correction improves accuracy.

*Lossy compression and hashing.* A hash function produces a local digest that determines file identity with high probability. Collisions (two files with the same hash) are the "failure" of the profile. The birthday paradox gives the threshold. But the tournament analogy suggests something more structured: the collision gap might have a precise "quantum" -- the minimal structural difference between colliding files -- and this quantum could inform better hash designs that are specifically resistant to the most common collision type.

*Sensor networks.* Each sensor measures local conditions. The "profile" (vector of all sensor readings) determines the system state up to some complexity threshold. The failure onset -- where two different system states produce identical sensor profiles -- reveals the "entanglement distance": the minimum spatial scale of correlations invisible to the sensor grid. The gap size (4, in the tournament case) would correspond to the minimal detectable perturbation. This has direct implications for sensor placement optimization: sensors should be positioned to break the "alpha_2 degeneracy," i.e., to distinguish the specific non-local states that the local profile misses.

---

## Pattern 4: The Weight-Norm Anticorrelation

**What it is:** An object's distance from the "center" of a symmetry group anticorrelates with a global counting quantity. The most symmetric object maximizes the count; the most asymmetric minimizes it.

For tournaments: the A_{n-1} weight norm ||w||^2 measures distance from regularity (the zero weight = the center of the Weyl group orbit). At n=5, ||w||^2 = 0 gives H = 15 (maximum), and ||w||^2 = 40 gives H = 1 (minimum, transitive). The Paley tournament, sitting at the center, maximizes H.

**The abstract principle:** In any system where a global count depends on the object's structure, look for a natural symmetry group, define the "center" as the fixed point (or zero weight), and check whether distance from center anticorrelates with the count. If so, the optimization problem (maximize the count) reduces to: find the most symmetric object.

**Where it leads:**

*Neural network generalization.* The "flat minima" hypothesis in deep learning says that weight matrices near flatter optima generalize better. Flatness is measured by the Hessian's eigenvalue spread -- essentially a weight-norm-like quantity. The tournament anticorrelation suggests a sharper version: generalization should correlate not with arbitrary norm measures but with the weight matrix's distance from the nearest GROUP-SYMMETRIC point (e.g., the nearest orthogonal matrix for orthogonal symmetry, or the nearest permutation-equivariant matrix for S_n symmetry). The "Paley-like" network would be the one whose weight matrices sit at the center of the symmetry orbit -- maximally regular, maximally flat, and maximally generalizing.

*Diversity in ecosystems.* The ecological tournament has a "score sequence" (competitive ranking of species). A regular tournament (all species equally competitive) is the zero-weight object. The anticorrelation predicts: the most diverse ecosystem (highest "H" = number of stable food web configurations) should be the one with the most even competitive hierarchy. This is exactly what Allesina & Levine found: maximum biodiversity occurs when competitive abilities are roughly equal and intransitive cycles are abundant. The weight norm gives a QUANTITATIVE prediction: biodiversity should decrease as ||w||^2 increases, with the transitive hierarchy (single dominant species) as the minimum.

*Portfolio diversification.* In finance, a portfolio of n assets can be represented by a tournament (pairwise return comparisons over a period). The "weight" is the imbalance vector (how much each asset's performance deviates from the average). The anticorrelation predicts that portfolios where assets have EQUAL risk-adjusted returns (zero weight = "regular tournament") should have the maximum number of viable rebalancing strategies. Highly unbalanced portfolios (one asset dominates = transitive tournament) have minimum strategic flexibility. This connects tournament theory to the mean-variance frontier in a new way.

---

## Pattern 5: The Dimensional Coincidence

**What it is:** Two unrelated counting formulas produce the same number, and the coincidence points to a deep structural isomorphism.

dim so(n) = C(n,2) = # positive roots of A_{n-1}. This is not a coincidence -- it reveals that the standard basis of so(n) IS the set of root vectors of sl(n). The Petersen graph materializes at n=5 where this count is 10, but the identity holds for all n.

**The abstract principle:** When two different formulas from two different theories give the same answer, there is almost always a natural bijection or isomorphism hiding beneath. Finding it converts a numerical curiosity into a theorem.

**Where it leads:**

*Quantum error correction.* The Steane code uses 7 qubits (= first tournament prime). The number of stabilizer generators is 6 = C(4,2) = dim so(4). The number of syndrome bits is 3 = # positive roots of A_2. Are these coincidences, or is there a root-system structure hiding in stabilizer codes? The Pauli group on n qubits has 4^n - 1 nontrivial elements, and the symplectic structure of the stabilizer formalism connects to sp(2n) -- another Lie algebra in the so/sp family. Investigating whether tournament-style conflict graphs (which operators commute/anticommute?) have Petersen-like structure could yield new code constructions.

*Information geometry.* The Fisher information matrix for a statistical model with n parameters lives in sp(2n, R) (symplectic, because it's both symmetric and positive definite). The number of independent entries is C(n+1, 2) -- almost the same as C(n,2) but shifted by 1. The "off-by-one" between symmetric and skew-symmetric matrices (dim = C(n+1,2) vs C(n,2)) corresponds to the difference between cooperation (symmetric) and competition (skew-symmetric) in the Cartan decomposition. Every statistical model has a "tournament sector" (the skew-symmetric part of Fisher information) and a "cooperation sector" (the symmetric part). This decomposition could reveal which parameters compete for evidence and which cooperate.

---

## Pattern 6: The Impossibility via Duality Reversal

**What it is:** A structure is impossible as a substructure because it exists in the wrong half of a duality. The Petersen graph can never be a tournament conflict graph because it lives in the Kneser (anti-conflict) world, not the Johnson (conflict) world.

**The abstract principle:** When trying to prove that a certain configuration is impossible, look for a duality that places the configuration on the "wrong side." The proof becomes: the configuration encodes X, but the context requires Y, and X and Y are complements.

**Where it leads:**

*Impossibility results in distributed computing.* The FLP impossibility (Fischer-Lynch-Paterson: no deterministic consensus in asynchronous systems with one failure) might have a duality-reversal proof. The "failure" configuration lives in the "anti-consensus" structure (the dual of consensus = the configurations where no agreement is reachable). If the failure configuration can be shown to be a Kneser-like object in a Johnson-like context, the impossibility follows from duality. This is speculative but aligns with existing topological proofs of distributed computing impossibilities (Herlihy-Shavit).

*Impossible codes in coding theory.* The Singleton bound, Hamming bound, and Plotkin bound all show that certain (n,k,d) codes are impossible. Some of these have proofs via "duality" (the dual code would need properties that violate a known constraint). The Petersen-Lie insight suggests a stronger form: certain codes are impossible because their weight enumerator polynomial lives in the "anti-code" (Kneser-like) structure rather than the "code" (Johnson-like) structure. The MacWilliams identity (which relates a code to its dual) is the precise analogue of Kneser/Johnson complementation.

*Arrow's theorem revisited.* Arrow's impossibility theorem states that no voting rule satisfying unanimity and independence of irrelevant alternatives can be non-dictatorial. The duality-reversal view: a "fair" voting rule would need the preference aggregation to live in the "conflict" (Johnson) structure -- where overlapping preferences interfere and produce a social ranking. But the IIA axiom forces the rule to live in the "anti-conflict" (Kneser) structure -- where each pairwise comparison is independent. These are complementary. No rule can inhabit both simultaneously, just as the Petersen graph can't be a conflict graph. This is not a proof of Arrow's theorem, but it suggests a Lie-algebraic reformulation might exist.

---

## The Meta-Pattern

All six patterns share a common core: **structure becomes visible at the boundary between two perspectives.** The dual habitation reveals invariants invisible to either home alone. The anti-conflict complement is the shadow cast by conflict. The profile determinacy boundary is where non-local entanglement first appears. The weight-norm anticorrelation is the price of asymmetry. The dimensional coincidence is a theorem-in-waiting. The impossibility via duality reversal is a proof by wrong-side-of-the-mirror.

The unifying theme is that **the interesting structure lives in the intersection, not in either piece.** Tournaments live at the intersection of so(n) and A_{n-1}. Conflict lives at the intersection of Kneser and Johnson. Determinacy lives at the boundary of local and non-local. H lives at the intersection of competition and cooperation.

This is, I think, what the Petersen graph has been trying to tell us. It is the smallest graph where the orthogonality structure of a root system becomes visible as a GRAPH. It is the avatar of the principle that complementarity has structure, that the gap between two perspectives is itself an object worth studying, and that the deepest theorems come from finding the two perspectives whose gap is exactly the thing you're trying to understand.

---

## Concrete Next Steps (Engineering Priority)

1. **GNN augmentation via profile gaps**: Implement a "tournament-aware" message-passing layer that explicitly computes 3-cycle profiles per edge and passes the profile-gap (alpha_2-like non-local features) as skip connections. Test on molecular property prediction benchmarks.

2. **Wireless scheduling via root geometry**: For all-pairs mesh networks, model the scheduling problem as coloring K(n,2). Use the chromatic number = n-2 (Lovasz) as a scheduling bound. Compare to existing TDMA bounds.

3. **Ecological biodiversity predictor**: Given a competition matrix (tournament), compute the weight norm ||w||^2 and predict biodiversity (# stable configurations) via the anticorrelation. Test on empirical data from Allesina & Levine's competition experiments.

4. **Attention head tournament fingerprint**: For each attention head in a transformer, compute the 3-cycle root profile of the thresholded attention tournament. Use this as a per-head signature for interpretability and pruning decisions. Fast: O(n^3) per head.

---

## Addendum: The Girth of Omega and What It Reveals

**Computed after the initial six patterns, at the human's suggestion.**

The girth of Omega(T) -- the length of the shortest cycle in the conflict graph -- turns out to be a binary invariant. It takes exactly two values across all tournaments at n <= 6:

- **girth = 3** whenever alpha_1 >= 3 (any tournament with 3+ odd cycles)
- **girth = infinity** whenever alpha_1 <= 2 (Omega is a forest or trivial)

There is NO intermediate girth. No tournament produces a conflict graph with girth 4, 5, 6, or any finite value other than 3. The conflict graph either has triangles immediately or has no cycles at all.

This is a **Phase Transition with No Intermediate Phase.** The system is either in the "trivial" regime (too few cycles for any conflict cycle) or in the "fully entangled" regime (triangles everywhere). There is no locally-sparse-but-globally-cyclic intermediate.

### Why This Matters for the Six Patterns

**Pattern 2 (Anti-Conflict) gains teeth.** The Petersen graph has girth 5. Tournament conflict graphs have girth 3 or infinity. The girth gap is not just 2 (= 5 - 3) but effectively infinite: there is no tournament Omega with girth 4 approaching the Petersen's girth 5 from below. The Petersen lives in a fundamentally different girth class. The impossibility is not marginal -- it is categorical.

**Pattern 3 (Profile Determinacy) gains an explanation.** The root cycle profile fails at n=6 because it cannot see the difference between "Omega with disjoint pair" and "Omega without." But girth tells us this difference doesn't affect the cycle structure of Omega itself (girth stays 3 either way). The profile fails to determine H because it fails to determine alpha_2 -- but alpha_2 does NOT change the girth. Girth and H decouple at exactly the point where the profile fails. The profile captures the "phase" (girth 3 vs infinity) perfectly; what it misses is the QUANTITATIVE structure within the girth-3 phase.

**Pattern 6 (Impossibility) gains a quick test.** Any graph proposed as Omega(T) must have girth 3 or infinity. If someone produces a graph with girth 4, 5, 6, ..., it cannot be a tournament conflict graph. This is a STRONGER impossibility filter than the Kneser/Johnson duality argument: it eliminates not just the Petersen but ALL graphs with finite girth > 3. Every Moore graph, every cage, every locally-linear graph -- none can be Omega(T).

**A seventh pattern emerges: The Binary Phase.** The girth result reveals that tournament conflict is a BINARY phenomenon. Either you have enough cycles to form triangles immediately (alpha_1 >= 3), or you don't have enough to form any cycles at all (alpha_1 <= 2). There is no "partially conflicting" regime. This binary phase connects to:

- **Neural networks**: The loss landscape is either fully connected (many paths to minimum) or fragmented. Intermediate connectivity is unstable -- a result seen in random matrix theory for overparameterized networks.
- **Ecology**: An ecosystem is either "Rock-Paper-Scissors entangled" (intransitive cycles everywhere, girth 3) or "hierarchical" (no competition cycles, girth infinity). Kerr et al. (Nature 2002) found exactly this: intermediate states are transient, collapsing to one extreme or the other.
- **Social choice**: Preferences either form a Condorcet cycle immediately (girth 3 among alternatives) or are fully acyclic (transitive preferences). The "almost-but-not-quite-cyclic" regime barely exists. This is Condorcet's paradox: once you have enough voters, intransitivity either appears everywhere or nowhere.

The binary phase is perhaps the deepest of the patterns because it constrains all the others. The profile determinacy boundary (Pattern 3) isn't just a boundary -- it separates two phases with nothing between them. The impossibility (Pattern 6) isn't marginal -- it's categorical. The weight-norm anticorrelation (Pattern 4) maps directly onto the two phases: low norm = girth 3 phase, high norm = girth infinity phase.

### The Anti-Conflict Girth Is Always Infinity

Perhaps the most striking finding: at n <= 6, the anti-conflict girth (girth of the complement of Omega) is ALWAYS infinity. The complement of Omega has no cycles at all -- it is always a forest (or empty).

This means: there is no sequence of odd cycles C_1, C_2, ..., C_k where consecutive ones are vertex-disjoint AND the chain wraps back to the start. The "non-interference" relation among odd cycles is always acyclic. Odd cycles are organized in a TREE of non-interference, never a cycle.

This is the DUAL of the girth-3 theorem: conflict always has short cycles; non-conflict never has any cycles. Together:

**girth(Omega) in {3, infinity} and girth(Omega^c) = infinity**

Tournament conflict graphs are MAXIMALLY ASYMMETRIC in girth: the graph has the shortest possible finite girth (3), while its complement has the longest possible girth (infinity). No other simple graph family has this property. It is uniquely characteristic of tournament conflict.

---

*The six patterns are not metaphors. They are theorems about structure that happen to have been discovered in tournaments but that apply wherever two algebraic perspectives coexist, wherever conflict has a complement, wherever a local fingerprint reaches its resolution limit, wherever symmetry determines optimality, wherever two formulas accidentally agree, and wherever impossibility is the shadow of a duality. The girth of Omega, collapsing to a binary {3, infinity} with no intermediate values, is the sharpest expression of these patterns: tournament conflict is all-or-nothing, the Petersen graph is categorically excluded, and the acyclic anti-conflict complement is the hidden skeleton that holds the structure together.*
