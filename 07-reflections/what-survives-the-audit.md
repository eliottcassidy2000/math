# What Survives the Audit

**Session:** kind-pasteur-2026-03-21-S18z
**Arising from:** The protein folding reflection (S18y) and the need for rigorous self-correction

---

## The Audit

The S18 arc has produced 25 reflections, 5 computations, and dozens of claimed connections. Some are rigorous. Some are creative but ungrounded. Some are numerological. This reflection sorts them.

---

## TIER 1: Rigorous Mathematics (Provable, verified)

These survive any audit. They are theorems or verified computations.

**THM-261: Petersen = A_4 root orthogonality graph.** Proved. Verified computationally. The positive roots of A_{n-1} are in bijection with edges of K_n, and orthogonal roots correspond to disjoint edges = Kneser graph K(n,2). This is textbook-level algebra that was not previously stated in our context.

**THM-263: Root cycle profile determines H at n=5.** Exhaustively verified. Sharp: fails at n=6 with gap exactly 4 = one alpha_2 quantum.

**THM-264: Girth of Omega is {3, infinity}.** Exhaustively verified through n=6. Structural argument for all n.

**The cascade identity: g_{p+1}(rho_p) = -1.** Proved algebraically. Each k-nacci root is one quantum spherical for the next theory.

**g_p(2) = 1 for all p >= 2.** Proved: 2^p - (2^p - 1) = 1. The evaluation point x=2 is universally one quantum hyperbolic.

**The even-alphabet parity law.** Proved: I(G, q) is always odd when q is even, because q^k is even for k >= 1 and alpha_0 = 1.

**H(T_11)/C(11,2) = 95095/55 = 1729.** Verified computation. 1729 = 7 * 13 * 19 is factual. Whether the factorization is structurally deep is a separate question.

## TIER 2: Strong Structural Analogies (Real patterns, not fully proved)

These have genuine content but the connection to the specific algebraic framework is unproven.

**Tournaments as codes.** The partition function bridge (H at x=2 vs distance at x->0) is a real mathematical parallel. The score sequence IS a linear projection losing information, like a syndrome. The OCR IS a fraction of explained variance. Whether this rises to the level of a THEOREM (like Greene's theorem) requires formalizing the matroid structure.

**The binary skeleton.** The 26 binary phenomena are individually verified. Their CLASSIFICATION into types matching the PSL(2,Z) generators (2s and 3s) is a pattern that may be structural or may be an artifact of our choice of categories.

**Cotranslational folding as non-associative.** The biology is confirmed (2024-2025 papers). Domain assembly order matters. The label "octonionic" is creative but ungrounded — non-associativity does not automatically mean octonion algebra.

**The CD tower in the transformer.** The four weight matrices per head is architectural fact. Quaternion transformers with 75% savings are published. Inter-head coupling is an active frontier matching CD level 3. The SPECIFIC claim that the coupling should follow the Cayley-Dickson doubling formula is testable but untested.

**dim(so(4))/dim(p) = 6/9 = 2/3.** This is pure linear algebra applied to the Cartan decomposition of gl(4,R). The 2/3 is exact. Whether it explains the Napolitano phase transition or the T_11 transitivity requires additional evidence.

## TIER 3: Suggestive Connections (Interesting but may be coincidental)

**1729 = j(i) + 1.** True numerically. But j(i) = 1728 is a famous number in its own right, and many numbers are "1728 + something." The factorization 7*13*19 landing on tournament primes IS striking, but the prior probability of any 4-digit number factoring into three primes from a set of 7 is not negligible.

**The Bernoulli-Coxeter chain.** Each tournament prime enters at weight h = p-1, which is a Coxeter number. This is TRIVIALLY true: every integer >= 2 is h(A_k) for some k. The NON-trivial claim is that the exceptional primes {7, 13, 19, 31} correspond to exceptional Lie algebras. This is TRUE (h(G_2)=6, h(E_6)=12, h(E_7)=18, h(E_8)=30 and 6+1=7, 12+1=13, 18+1=19, 30+1=31 are prime). But the INTERPRETATION (exceptional = hard obstructions) is a framing choice, not a theorem.

**The dimension axis.** The k-nacci ladder is real mathematics. The assignment of "tessellation dimensions" to transcendental constants is a definition, not a discovery. The claim that sqrt(pi) at D=1.70 "means" something about the circle is interpretive, not mathematical.

## TIER 4: Likely Numerological (Creative but probably coincidental)

**3.6 = 18/5 = h(E_7)/F_1.** This is 18/5 = 3.6, which is arithmetic. The alpha helix parameter comes from backbone physics (bond angles, van der Waals radii), not from Lie algebra structure. Cherry-picked.

**Ramachandran 2/3.** The generously allowed region is approximately 70%, which rounds to 2/3. But the exact fraction depends on the van der Waals radii used and the definition of "allowed." The match to the Cartan ratio is approximate at best.

**Alpha helix = quaternion via NH/sidechain/CO/sidechain mapping.** The four-component structure of the helix repeat unit is real, but the specific mapping to Q,K,V,O is imposed rather than derived. The i+4 hydrogen bond comes from backbone geometry optimization, not from the Hamilton product.

**Misfolding = sedenion zero divisor.** A creative metaphor (non-zero*non-zero = zero-function). Not derived from any mathematical model of protein energetics.

---

## The High-Leverage Actions

Based on this audit, here is what should be PURSUED vs what should be ACKNOWLEDGED AS SPECULATIVE:

### PURSUE (Tier 1-2, actionable)

**1. Quaternion Evoformer.** Replace the 4 real weight matrices per attention head in a protein structure prediction model with quaternion-valued matrices. The Hamilton product captures Q-K-V-O coupling. This is:
- Directly implementable (quaternion transformer code exists)
- Well-motivated (4 matrices per head is architectural fact)
- Testable (compare protein structure prediction accuracy with real vs quaternion heads)
- Published precedent (quaternion transformers for NLP: ACL 2019, ICLR 2021; sedenion networks: 2025)

Expected outcome: 75% per-head parameter savings. Possible accuracy improvement from inter-component coupling. Low risk experiment.

**2. Order-aware domain assembly for multi-domain proteins.** Current protein prediction treats all residues simultaneously. Cotranslational folding proceeds N-to-C. A model that processes domains sequentially (respecting synthesis order) and allows earlier domains to influence later ones (without the reverse) would capture the non-associative assembly biology.

This does NOT require octonion algebra. It requires a simple architectural change: mask the attention so that domain A can influence domain B but not vice versa (like autoregressive masking, but at the domain level instead of the token level). This is:
- Biologically motivated (cotranslational folding is real)
- Architecturally simple (domain-level causal masking)
- Testable on multi-domain proteins where cotranslational effects are known

**3. Tournament invariants as protein features.** Compute the tournament invariants (H, alpha_k, beta_k) of the directed contact graph (with N-to-C direction) for known protein structures. Correlate with:
- Folding rate (does H predict how fast a protein folds?)
- Stability (does the independence polynomial predict thermodynamic stability?)
- Misfolding propensity (do specific conflict graph properties predict aggregation?)

This is pure data analysis, no new theory needed.

**4. Formalize the tournament-code correspondence.** Prove or disprove: the sorted root cycle profile (SRCP) determines H for all tournaments. This is the central conjecture (verified through n=7 with (c3,c5) profiles). A proof would establish tournaments as codes with the SRCP as the weight distribution.

### ACKNOWLEDGE AS SPECULATIVE

- 3.6 = h(E_7)/F_1 (numerological)
- Ramachandran 2/3 (approximate, definition-dependent)
- Alpha helix = quaternion (structural similarity, not derivation)
- Misfolding = zero divisor (metaphor)
- PSL(2,Z) governs LLMs (ungrounded beyond linear algebra)
- Moonshine connections to biology (far beyond current evidence)
- The dimension axis assigns meaning to transcendentals (interpretive framework, not mathematics)

---

## What This Project Has Actually Achieved

Setting aside the speculative extensions, the S18 arc established:

1. **New theorems** (THM-261, 262, 263, 264) connecting tournament conflict graphs to root systems and Lie algebras.
2. **The cascade identity** g_{p+1}(rho_p) = -1 and g_p(2) = 1, proving the k-nacci hierarchy has universal structure.
3. **The even-alphabet parity law**, reframing Redei's theorem as a consequence of alpha_0 = 1.
4. **26 catalogued binary phenomena**, organized as a structural feature of tournament completeness.
5. **A concrete architecture proposal** (quaternion transformer for protein structure) grounded in published hypercomplex neural network results.
6. **Identification of cotranslational non-associativity** as a biologically real phenomenon matchable to the CD tower's level 3.

The speculative reflections (tessellation, modular group, moonshine, Vitali atoms, dimension axis) are valuable as FRAMEWORKS FOR GENERATING HYPOTHESES. They are not valuable as established mathematics. The discipline is to know which is which.

---

*This reflection does what a morning after a creative night should do: sort the gold from the glitter. The gold is: new theorems about root systems and conflict graphs, the cascade identity, the even-parity law, and a concrete proposal to build quaternion protein-structure models. The glitter is: everything about moonshine, the Ramachandran plot, 3.6 = 18/5, and the dimension axis. The glitter inspired the gold — we would not have found the cascade identity without the creative push through the dimension axis. But the gold stands alone, and the glitter does not. This is what honest mathematics looks like: the creative night generates; the sober morning audits; and the afternoon builds on what survived.*
