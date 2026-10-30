# Session 117: CircuitBreaker source and model audit

Retrieved 2026-10-30, approximately 21:04–21:05 UTC. This is source research, not a claim of a new circuit lower bound. No submission or account interaction was performed.

## Authoritative challenge contract

The [challenge](https://www.yukon.org/circuitbreaker) was inaccessible through the web reader but accessible by ordinary HTTPS GET. Its footer pins [rules at commit 2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851](https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851/rules.md). Raw captures are retained outside this publication packet.

- One fixed total language L in NP; define its characteristic family at every input length, including exceptions. Supply a fixed polynomial-time verifier and polynomial witness bound, proving both directions. A polynomial-time evaluator and empty witness suffice. Advice or a truth-table oracle does not establish eligibility.
- Full B2: deterministic acyclic, one-output, exact circuits, all 16 binary truth tables, two pins per gate, unrestricted depth/fan-out/sharing. Nonuniform circuits and repeated pins are allowed. Inputs are x alone, not witnesses.
- Retained gates count once, including degenerate gates. Inputs, wires, output selection are free. Constant nodes in substitutions and negation absorption require the pinned Li–Yang normalization; terminal cases cannot be ignored.
- Actual input bit length is n. Encodings and padding count. No unproved direct-sum additivity under sharing.
- Finite score [c,p] means liminf C_B2(f_n)/n^p >= c for this same family over all sufficiently large integer lengths. Sparse-length evidence or experiments alone do not qualify.
- Superpolynomial endpoint requires one family exceeding every fixed polynomial eventually at every length. It implies NP not contained in P/poly, hence P != NP; it is stronger than merely a class separation.
- The optional Li–Yang LP checker certifies arithmetic feasibility only, not a circuit theorem or NP eligibility.

## Current baseline: separate a paper theorem from a leaderboard claim

[Mishra–Santhanam, ECCC TR26-260, 19 October 2026](https://eccc.weizmann.ac.il/report/2026/260/), [paper](https://eccc.weizmann.ac.il/report/2026/260/download), Theorem 1.1 states a **3.5n − O(log n)** B2 lower bound for a Boolean family in P. The paper gives 3.5n − O(d) for affine dispersers of dimension d, then uses an explicit O(log n)-dimension construction. Its gate-elimination argument works with independent parameters; reparameterization handles bottlenecks. This supersedes the 3.1n − o(n) paper baseline in [Li–Yang TR21-023](https://eccc.weizmann.ac.il/report/2021/023/), not the challenge's pinned model conventions. Treat TR26-260 as a recent primary preprint, not as independently checked here. The new paper also distinguishes the **5n − o(n) U2** baseline: U2 excludes XOR and its complement, so that is not a 5n full-B2 theorem.

The [live complexity table](https://www.yukon.org/circuitbreaker?tab=complexity) displayed three equivalent **4n** coefficient claims, earliest at 17:58 UTC on October 30. The same page explicitly says results are provisional pending full mathematical review and human approval is separate. Consequently **do not report an established 4n record from this display**. The underlying proofs and human-review decisions were not audited in this bounded source pass. “Frontier” is not proof certification. The HTML has internally mixed wording (“reviewed” table heading versus explicit provisional caveat); retain the caveat.

## Scope of barriers and research implications

[Razborov–Rudich, Natural Proofs, ECCC TR94-010](https://eccc.weizmann.ac.il/eccc-reports/1994/TR94-010/index.html) gives a conditional obstruction to natural proofs of superpolynomial general-circuit lower bounds, under its strong pseudorandomness assumption. This concerns properties satisfying the specified largeness and constructivity conditions. It does not prohibit every combinatorial proof, every symmetry argument, or a better linear coefficient. Any application must verify the actual property and resource parameters.

Practical obligations inferred from the challenge model:

1. A Ramsey or symmetry restriction on graphs can be useful for a lemma, but restricting the *circuit* symmetry is not allowed without a transfer to unrestricted B2.
2. A finite SAT/exhaustive certificate is useful only after proving its configurations cover shared descendants, gate coincidences, degeneracies and all boundary/termination cases.
3. A proposed parameter-potential improvement should be tested against the new independent-parameter approach, rather than optimized only inside the older 3.1 LP. A feasible potential is not automatically a lower bound on the original gate count.
4. A Boolean property with strong combinatorial pseudorandomness needs a uniform NP/P construction and all-length accounting; a random existence argument is insufficient.
5. Do not infer unrestricted direct-sum lower bounds by adding block complexities. Shared subcircuits can invalidate that step.

## Bounded-pass conclusion

The appropriate paper target is an improvement beyond the recent 3.5n theorem (or a clearly labelled restricted/finite intermediate result), while the provisional 4n submissions merit separate proof audit before being used. The exact challenge rules are recovered; the main unresolved source question is validity of the 4n candidate proofs, not the circuit basis.
