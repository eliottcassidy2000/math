# Prior planar groups: sourced inheritance for Chan–Pak

Source snapshot: current worktree at `001b4aee22c8`, 2026-11-10. This is a provenance and proposed-map report, not a new coverage theorem. No computations or external literature verification are claimed here. Paths below are repository-relative. Chan–Pak completeness remains OPEN.

## 1. Resolve the name before importing a theorem

There are at least three genuinely different prior meanings of planar in this repository. It would be incorrect either to silently identify them or to assume that the user meant only transformation groups.

1. **Explicit planar groups in the recent tournament lineage.** `05-knowledge/results/thomassen_merge_hypohamiltonian_groups_20261104.md`, sections 3–4, discusses Maschke's planar Cayley groups and tournament symmetry. Its section 4 derives, from the cited classification and the odd order of tournament automorphism groups, that an ordinary tournament automorphism group is planar-Cayley exactly when cyclic. It records C3=>C3 with C3×C3, P7 with F21, and signed/converse groups D3, D5, AGL(1,7). This is the closest literal match to “prior work on planar groups” in the current Chan–Pak context.
2. **Planar source Hamiltonian transformation groups.** The October19 packets below give a particularly strong second inheritance: an exact group/state calculation, finite-word obstructions, failure of those obstructions under closure, and disappearance of the obstruction after enlarging generators.
3. **Algebraic planar functions on finite abelian groups.** THM-2356 and THM-2363 concern bijective derivatives, Fourier/Gram reconstruction and phase loss. They are neither planar Cayley groups nor plane transformation groups.

`00-navigation/S190.md`, line 29, routes to `05-knowledge/results/session190_constructive_links_20261104/groups/session190-groups/PROOF.md` and `groups/session190-groups-audit/AUDIT.md`. These are tracked but omitted by this sparse checkout, not missing repository artifacts. Both were recovered using `git show HEAD:<path>`. The audit is **PASS, no mathematical correction required**. Likewise the tracked Thomassen audit at `04-computation/experiments/thomassen_merge_20261104_audit_S10/AUDIT.md` was recovered and reports SOUND WITH CORRECTIONS, incorporated in the maintained packet. This pass read those proofs/audits but did not rerun their historical computations.

### Recovered Session190 group proof and audited boundaries

The original proof header still says independent review pending, but its separate audit supersedes that stale timing statement. Its core results are:

- Finite planar Euclidean point groups reduce by orbit averaging to O(2), hence cyclic or dihedral; SO(3) has the five standard families. Tournament automorphism groups have odd order, so any faithful planar or SO(3) realization is cyclic. This concerns geometric representations, distinct from planar Cayley graphs and underlying graph planarity.
- Paley7 ordinary/signed groups are F21 and AGL(1,7); the actual point/unordered-pair/ordered-pair orbits have sizes7/21/42. They do not explain forbidden Hamilton-path counts. H(P7)=189 and its endpoint matrix is1 on arcs,8 on reverse arcs.
- For an actual labelled path-complement metagraph edge T--T_P, ordinary unordered-edge stabilizer is trivial: common arcs recover the entire directed spanning path. Signed stabilizer is1 orC2, with unique possible nonidentity given by global converse followed by path reversal. That signed element cannot swap the edge endpoints. Every neighbour breaks every original nonidentity ordinary automorphism, but may gain new automorphisms.
- F21 has minimum faithful real representation dimension6 (character orbits of C7 plus complex conjugates), attained in SU(3); its minimum faithful binary dimension is3. Reduction of skew signs modulo2 loses tournament orientation. Maschke requires characteristic not dividing group order and is not a geometric-group classification.

The audit independently verifies proofs, classical source statements, and finite controls. Its all-order edge-orbit formulas concern the specified labelled incidence carrier, not collapsed simple metagraph edges, connectivity, or raw Hamiltonian spectrum. These distinctions make the actual operation carrier a stronger bridge than group-order numerology.

## 2. Literal planar-group lane: what transfers now

**Inherited source and status.** The Thomassen packet section 4 states the Maschke classification as CITED and the odd-order cyclic consequence as PROVED from it. This pass has not independently rechecked the historical theorem's exact formulation. For a new result needing it, cite that dependence explicitly or prove only the needed cyclic implication. Do not relabel the five rotation-group families as finite subgroups of O(2).

**Closest proved Chan–Pak mechanism:**

- `01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`, title “Every automorphism subgroup retains the 7/21 path-orbit gaps; small symmetry has a sharp threefold cost.” Its full theorem and proof are present. For every actual A<=Aut(T), H(T)/|A| counts path orbits and is never 7 or21. The proof uses the inherited THM-4812 vertex-transitive classification and odd subgroup indices.
- `01-canon/theorems/THM-4990-stabilizer-corrected-hamilton-prefix-orbit-transfer.md`, title “Integral path-orbit transfer requires pointwise stabilizers and becomes uncorrected at the last two unused vertices,” gives the exact integral transition W(b,a)=N(b,a) k_a/k_b, with k the pointwise fixer of the used set. This is the appropriate operation-level interface for retaining a cyclic planar subgroup; merely dividing a transfer matrix by |A| is wrong.

**Concrete map:** start with (T,A), where A is an identified cyclic subgroup of Aut(T). Map ordered Hamilton paths to their free A-orbits, retaining the action and prefix stabilizers. Preserved predicate: existence and exact count of path orbits; H=|A|q. Lost information: individual labelled paths and root-specific behavior unless a rooted action/sidecar is retained. Cheapest check: C3 after a source has N=3,k_b=3,k_a=1 and therefore W=1, whereas naive transitions give 3. This is an inherited hostile, not a new discovery.

**Critical boundary from Session237:** THM-5092 proves gcd(h,R,e,d) is the gcd of all fixed-port S/B descendant counts, using empty,B,SB,BB. An actual rooted automorphism group divides this content, but equality need not hold. Its rigid five-vertex witness at root1 has (9,9,3,3); all fixed-port descendants are divisible by3, but reroot+B gives19 and converse+B gives29. Thus arithmetic normalization by3 is NOT quotienting by a hidden cyclic planar symmetry. The six-vertex same-interface collision further proves that even (order,h,R,e,d) does not determine reroot escape. This is the strongest immediate synthesis of the planar-group route with recent work.

**Canonical hostile:** Q7 has (H,|Aut|,q)=(175,7,25), while P7 has (189,21,9). Cyclic seven-symmetry does not force the low-quotient P7 classification, and the appearance of7,21 as group orders is not an exclusion proof for H.

**Least-used relevant sidecar:** explicit action on labelled vertices plus pointwise stabilizers, rather than automorphism order alone. For rerooting one needs all-root interfaces and their incidence, not one fixed-root vector.

## 3. Transformation-group lane: exact recovered results

All files here are under `05-knowledge/results/` and dated 20261019.

### A. `planar_common_u_rational_flows_20261019.md` — PROVED

For characteristic-zero K, u=x²t and H=f(u) with f'(0)=0 give the source-completed flow
X=x/(1−λxf'(u)), T=t(1−λxf'(u))².
The invariant u is preserved. Finite common-u words commute and collapse to the same formula with F=Σ λ_i f_i; polynomiality is equivalent to F'=0, hence identity. Exact verification artifacts and independent review are linked in the packet (169 controls, as reported there).

**Hostile:** f=u²(u+1)² yields rational nonpolynomial maps fixing both distinguished curves t=0 and1+u=0 pointwise. Boundary agreement does not imply global polynomiality/properness.

### B. `planar_single_hamiltonian_polynomiality_20261019.md` — PROVED FROM CITED INPUT

For H=c+t²R and λ≠0, a source Hamiltonian flow is polynomial iff H_x=0, equivalently its second coordinate is t. It is then a shear. The external dependency is Gwoździewicz, Injectivity on one line, Lemma2.3 (fixed-line Keller theorem); the packet identifies its source. A single-factor obstruction is not automatically a word obstruction.

### C. `planar_affine_x_rational_word_rigidity_20261019.md` — PROVED

For H_i=a_ixt²+f_i(t), a_i≠0, f_i∈t²K[t], use z=1/t,w=xt². A flow is (z,w)->(z+ell,w+[F(z)−F(z+ell)]/a), F(z)=f(1/z). A chronological word is (z+L,w+G(z)). The full state is L plus all pole/order coefficients C_(p,m) of G. For every order m, Σ_p C_(p,m)=0. Polynomiality forces L=0 and all finite poles except possibly0 to vanish; coefficient-sum conservation then kills the last pole. Thus finite words are polynomial iff identity.

**Mechanism:** do not quotient the word state to its scalar total shift. Coincident poles and cumulative shifts must be retained; word reversal changes them. The exact finite identity test is L=0 and every C_(p,m)=0.

### D. `planar_source_closure_formal_shears_20261019.md` — PROVED

The affine-x group's t-adic closure is Phi_(L,g)=(x(1+Lt)²+g(t), t/(1+Lt)), g∈tK[[t]]. Its polynomial intersection is every shear (x+p(t),t), p∈tK[t], despite the finite group's identity-only intersection. Even two fixed Hamiltonians xt²,(x+1)t² generate a dense subgroup of this closure. Finite-difference return loops build each desired jet.

**Corrected-near-miss / hostile:** the finite principal-part sums are not continuous in source topology. The explicit G_N=Σ c_(j,N)/(z+j)² all have zero coefficient-sum by order yet converge at infinity to z^-3, whose order3 coefficient-sum is1. Therefore matching arbitrarily long finite jets is not finite realization, and passing to closure changes the target intersection.

### E. `planar_full_source_group_polynomial_intersection_20261019.md` — PROVED FROM CITED INPUT

Allow every x-dependent polynomial source Hamiltonian c+t²R. The finite generated group's polynomial intersection is exactly all source shears. Each nonzero g∈tK[t] is realized by three individually rational, nonpolynomial factors: J=x³t² at time−1, H_g=f(t(1−x²t)³) at time+1, J at time+1, where f'=g. Function and point-map orders are explicitly distinguished. Including the individually polynomial generators does not enlarge this full generated group.

**Canonical obstruction to overextension:** the affine-x principal-part rigidity holds for its specified alphabet, but fails under the cubic-x generator and conjugation. Negative times/inverses are essential to this realization. A positive graph-growth monoid cannot automatically perform the same conjugation.

## 4. Meaningful mappings to Sessions236–237, with falsifiers

### Bridge I: operation-state congruence, not just equal observables

Source: affine-x normal form (L,C_(p,m)); target: rooted tournament interface (h,R,e,d). Both give exact future responses only for their declared operations. Map is methodological, not a homomorphism between planar flows and tournaments. Preserved predicate: equality of all responses under that operation alphabet. Lost information: other operations, notably reroot/converse. Needed sidecar: all-root response incidence or full marked tournament. Cheapest decisive test is the existing THM-5092 equal-interface pair whose reroot behavior differs. This already refutes claiming that one fixed-port state is a universal sufficient state.

### Bridge II: rigidity must be quantified over an alphabet

Source: affine-x identity-only theorem versus enlarged full-source shear realization. Target: THM-5080 all-word tower-ray rigidity for the frozen15/57 ternary block matrices. Preserved predicate: exact intersection of an orbit/action with a target locus. Lost information: effects of generators outside the alphabet. Needed sidecar: actual matrix/physical action for each proposed enlargement and its invariant cone. Cheapest test: compute preimages of v=(3,1,2) for each newly allowed matrix and test physical-state constraints. THM-5080 already has an exceptional non-tower preimage (45,1,5) excluded only by h<=3O; “all non-towers have no nonnegative preimage” is false for the57 alphabet.

A conjugation search is reasonable only when the inverse steps have a physical, authorized combinatorial realization. It is not justified to invert an integer transfer matrix formally and call the result a tournament operation.

### Bridge III: finite realization versus local/limit reachability

Source: t-adic closure contains shears absent from finite words. Target: THM-5091 residue surjectivity with exactly2M² states for the normalized full-count observer, and THM-5080 exact rational cones. Preserved predicate: finite modular/linear observations; destroyed predicate: exact integer orbit membership and finite operation cost. Needed sidecar: one coherent finite witness, quantitative height/cost bound, and compatibility across observation levels. Cheapest hostile: (10,4,5) satisfies Session236's exact nonempty cone and h=1 mod3 but h=10 is a permanent ternary gap. Moduli or closure cannot supply Chan–Pak completeness by themselves.

### Bridge IV: actual symmetry quotient versus formal content

Source: literal planar cyclic subgroups and free actions. Target: THM-5092's content. Preserved predicate only when an actual rooted action exists: divisibility and orbit interpretation. Destroyed information under gcd division: the action itself. Sidecar: permutations realizing the action and pointwise stabilizers. Cheapest hostile: rigid state(9,9,3,3). This is the strongest concrete theorem-level connection for the user's literal planar-group request.

## 5. Algebraic planar functions: useful but secondary

`01-canon/theorems/THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing.md` is PROVED with exact independent audit: complete chirp intensities plus labelled singleton energies reconstruct one signal up to a global phase. Separate planar graph signals retain independent phases; they do not determine the coarse target sum. `01-canon/theorems/THM-2363-planar-graph-detector-dominates-word-support-energy.md` supplies a one-way quadratic energy bound, not the converse or a physical target probe. MISTAKE-261 and ACTIVE-GUARDRAILS27d explicitly retain the cross-graph phase obstruction.

Potential connection to THM-5093 pair energy: use full pairwise/binary-position data rather than only an aggregate energy. No Fourier-mask realization of tournament operations or new Chan–Pak implication has been constructed. This is an analogy awaiting a map, not a proved bridge.

## Recommendation

Lead with the literal Maschke/cyclic symmetry lineage and its rigorous existing action-to-orbit transfer. Merge the transformation-group sequence as a second, sharply scoped inheritance: operation-relative obstruction, exact cancellation, generator enlargement, and closure failure. The next productive test should retain an actual action or all-root sidecar and challenge a fixed-alphabet obstruction, rather than accumulating another residue or cone census. Do not claim that any of these results settles Chan–Pak.
