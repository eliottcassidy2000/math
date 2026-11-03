# Five ways to keep the information that a clock forgets

Session173, 2026-11-03. **Synthesis of five independently audited packets.** The mathematical statements below inherit their proofs, hypotheses and finite controls from the linked packets. This note introduces no new global theorem, no spectrum-completeness claim, and no priority claim.

The shared question is constructive: when a numerical pattern repeats, what additional information lets us compose it into an actual object? The answer differs across the five settings. A phase, an integer displacement, an occupied vertex, a forbidden-shift mask, and a path-cover profile are not interchangeable.

## 1. Start with the actual provenance, then declare the model

The prior repository note `05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`, §3, records the marked valuation chain

    13 → 2 → 5 → 3 → 7 → 4,

coming from prime-to-valuation outputs of the displayed Monster/Baby Monster order ratio. Closing it requires the extra identification of the final exponent4 with the marked prime13. The source does not supply an ordinary six-state dynamical system or an integer-polynomial rule modulo9. In particular the full valuation list does not descend modulo9: primes11 and2 have the same residue but different successor residues.

The [clock packet](../clocks/PROOF.md) makes two explicit choices rather than silently repairing this loss. First, L(x)=527x+44 modulo72 has CRT factors −u+4 modulo8 and5r+8 modulo9. The identity5+8=13 appears literally at r=1, with13 displayed as digital root4. Second, the marked five-cycle becomes B=(4 2 5 3 7), fixing the other four mod9 states. Its extension away from the five marked states is chosen, not recovered from the valuation data.

These choices yield a substantial theorem: A(r)=5r+8 and B generate A9, and their72-state lift generates C2×A9 with exactly five reachable sectors of sizes9,9,18,18,18. The proof is elementary and explicit; complete finite group enumeration independently supports it. The affine map L alone has cycle type1²2¹¹6⁸. Its positive-integer version has no positive periodic orbit.

The distinction is essential. The useful new operation is a noncongruence-compatible switch added to an affine clock. This proves a group-generation result for the declared model, not a Monster theorem or a theorem about tournament Hamiltonian-path counts.

## 2. Relative phase and integer lift are different missing coordinates

Two synchronized finite cycles of lengths m,n have gcd(m,n) relative-phase orbits, each of length lcm(m,n). Thus common period does not identify a trajectory. This is already enough to prevent a false merger of two6-cycles.

The [square-sum packet](../squares/NOTE.md) adds a second distinction. A square-root word has an exact affine endpoint action x↦εx+d, together with signed occupied offsets. Mod8 and mod9 observe its reduction modulo72; they lose integer multiples of72. An actual certified path through[527] contains

    230 → 395 → 446,

whose square roots25,29 give translation216=3·72. Both clocks close, while the integer endpoints remain different. At height44 the combined residue label authenticates endpoint equality because reduction modulo72 is injective on the interval; it still does not certify integer squareness. For example28 is a square residue modulo8 and modulo9 but is not an integer square.

Retaining the exact roots, starting port, integer translation and occupied offsets supplies an explicit composable interface. The packet proves the full interval-placement criterion and the sharp statement that one exact root word has at most one Hamilton placement for odd interval size and at most two for even size. A36-vertex pair attains the latter bound; an alternating-block theorem classifies all double placements. These are placement results, not a new proof of the cited all-n square-sum existence range.

The same interface produces an optimal concrete connection of the motivating numbers:

    527–2–7–317–44–5–4–12–13–3.

Two added vertices are necessary, and{12,317} suffice;317 is also the least possible maximum addition when the minimum number of additions is used. Requiring8 adds a fourth old component; appending1–8 attains the corresponding optimum of three additions, again with maximum317. These are exact difference-of-squares constructions, not associations inferred from residues.

## 3. For runners, keep forbidden shifts and their overlaps

The [runner packet](../runner/result.md) uses the same numbers in an independently declared family. A common-divisor core stays fixed under shifts t₀+q/d. Each exceptional speed has a forbidden mask Bᵢ on those shifts. A maximum-weight spanning tree, with edge weights |Bᵢ∩Bⱼ|, gives an overlap-aware union bound. A prescribed cyclic order has a related bound with a full-intersection correction. The exact defects count disconnected selected subtrees or multiple selected cyclic runs, respectively.

This is an actionable bridge from clocks to a certificate: retain the masks or their Fourier vectors, then their Gram matrix supplies the required overlaps. Unrelated adjacency eigenvalues do not supply that matrix. The cycle(4,2,5,3,7,13) gives an exact successful9-sheet certificate. In the8-sheet example the prescribed cycle bound is inconclusive, while a better tree succeeds. The numbers44 and527 are valid composite moduli too; primality is unnecessary.

Here5 and13 on an8-sheet clock have the same step but shifted initial phases. The equality13=5+8 does not identify their forbidden events. The packet gives the precise extra phase condition under which one mask is a translate of the other.

The external status matters: the packet checked the statement and conventions of Allikvere's computer-assisted preprint [*Fourteen and fifteen lonely runners*, arXiv:2609.02604v2](https://arxiv.org/html/2609.02604v2). Its stated theorem includes both14 and15 total runners. That global computation was not independently rerun here. Our independently proved contribution is the composite-mask certificate and exact defect, not a new global lonely-runner proof; older repository descriptions of the target as open are historical context, not a reason to ignore the checked external claim.

## 4. A complete profile can compose correctly and still conceal states

The [tournament storage packet](../spectrum/NOTE.md) proves that the full path-cover profile is preserved by every fixed uniform-module substitution context. Equal-profile seeds therefore remain indistinguishable under a large, precisely defined family of compositions.

That does not imply the seeds are the same object. Cyclic substitution through a transitive feedback module constructs2^M pairwise nonisomorphic storage tournaments with identical complete profiles and arbitrary prescribed vertex-connectivity r, subject to the explicit side-size hypotheses. The unique minimum separating set recovers the intact ordered seed word. A fixed one-vertex, nonuniform attachment inside an addressed seed reads its bit through an exact Hamilton-count difference.

Two qualifications explain rather than weaken this result. Uniform substitution preserves the profile; the read head uses a more discriminating attachment. Connectivity survives fewer than r vertex deletions, but the stored bits need not survive deletion of a distinguishing seed vertex. Neither guarantee should be substituted for the other.

For the two-feedback family the explicit Hamilton counts have exact3-adic valuation M. They are therefore divisible by9 for M≥2, even though their mod8 residues have a separate simple pattern. The construction improves resilience and observability; its arithmetic itself shows that it does not fill the missing-count candidates7 or21.

## 5. The integer527 also has an exact composition-count interpretation

The [527 multiplicity packet](../counts527/RESULT.md) concerns a different observable: N(n,h), the number of ordinary tournament classes with actual Hamiltonian-path count h. It does not count states of L or path-cover-equivalent storage words.

Inherited audited inventories give A₁₇=2z⁶ and A₃₁=2z⁶+2z⁷. Since527=17·31, ordered strong-component composition yields

    F₅₂₇(z)=A₅₂₇(z)/(1−z)²+8(z¹²+z¹³)/(1−z)³.

The inherited total of184 strong527 classes fixes, for n≥11,

    N(n,527)=8n²+8n+1152−W,

where W=A′₅₂₇(1) remains explicitly unrecovered. The coefficients8,8 are known; the constant is not. This is a proved conditional consequence of the stated finite inventories, not an independent fresh census.

If44 denotes(H−1)/2, its actual count is89; if it denotes H directly, parity excludes it. For the former reading the complete inherited inventory gives N(n,89)=23n−149 for n≥8. These two polynomial sequences have exact eventual least periods9 and72 modulo72, respectively. The527 period is independent of the unknown W. This is a genuine arithmetic connection, but its clock advances tournament order n, whereas L advances a residue state. Sharing a modulus does not identify the transitions.

## 6. Constructive next targets, especially around{7,21}

The known exclusions7 and21 are inherited source results. None of the five packets reproves those exclusions or establishes that every other odd value occurs. The useful next questions specify the interface that would make a connection meaningful:

1. **Complete one census moment.** Recover W from the strong527 order histogram, retaining its exhaustive-completeness provenance. This closes the exact constant without repeating unrelated enumeration.
2. **Build an attachment calculus for actual counts.** Use endpoint and ordered path-cover data to search small nonuniform tournament extensions with an exact H-response law. A productive target is an explicit image theorem or congruence obstruction for a declared constructor, not another scalar pattern that ignores its domain.
3. **Escape the storage family's9-divisibility deliberately.** Change the seeds or feedback context, derive the new positive run-skeleton formula, and determine which odd values it realizes. The present construction's response machinery is available, but its current arithmetic cannot address7 or21.
4. **Demand an actual clock-to-tournament map.** Specify vertices, pairwise orientations, the permitted switch, and its action on Hamilton paths before using A9 or the marked cycle to explain a spectrum phenomenon. An explicit mismatch is useful information about the required extra state; an analogy alone is not a proof.
5. **Optimize real seams.** In square-sum constructions retain integer endpoint and occupancy information; in runner constructions retain translated mask overlaps. The optimal connector and mask examples show that these additional coordinates lead to constructive answers rather than merely more verification.

The common lesson is not that modular ideas fail. Each packet succeeds by stating what its observer preserves, retaining the specific information needed by its next operation, and proving the resulting composition law. Those laws are different, and that difference is what makes the combined research usable.
