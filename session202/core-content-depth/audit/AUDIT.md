# Independent CONTENT audit and separate TRIANGLE_67 audit

Decision: APPROVE the mathematical claims of CONTENT.md at SHA256 947b483f77b80013429ad3db2eb27fc20c19653fa16eccbd40f0710cc09fc64e within their explicitly restricted grammar. Separately APPROVE the finite-exact all-depth coprime-to-three target census in TRIANGLE_67.md at SHA256 ed6caf6d6c0563732b32eac301a3d017c46ab14f4b6d77e53c691ea1b5e40707. No mathematical correction is required. This audit does not itself promote or publish either packet, approve the entropy companion, or establish a universal tournament occupancy threshold.

## Proof review

The frozen ROOTED_ORDERED_AMALGAMATION.md Section 2 path-fragment argument gives the four physical transfer coordinates: the reverse-through-root case contributes (R-d)(L_C-H_C), while the other three disjoint cases contribute (h-e)H_C, (H(C)-start_C)d and e start_C. Substituting x=R-d and y=d-e yields the displayed h row with precisely the stated coefficients. Root-designated labelled path covers give the x row, and the non-singleton-core endpoint condition gives y and d. All coefficients are nonnegative: root-singleton covers inject into both R_C and L_C, endpoint paths inject into root-deleted paths, and root-deleted Hamilton paths inject into full paths by deterministic legal insertion. The singleton identity is not accidentally subjected to the non-singleton endpoint law. The two actual rooted edges give exactly O and Z.

For arbitrary integer c>=2, c divides H_C and beta_C is sufficient. At depth m the grades are d:m, x/y:m-1 (m>=1), h:m-2 (m>=2). The next core raises all grades as asserted; the first h constraint is exponent zero and needs no premise. O and Z preserve the grades, including the initial boundary. No finite-bank hypothesis or bound on spacer length enters this proof: an infinite bank is valid provided the SAME c divides both constants of EVERY allowed core. Each graph uses a finite word. This is not a claim about an infinitely long graph product.

Root deletion is an ordered join, so d is exactly the product of the core H_C values; singleton factors contribute one. Positivity of tournament Hamilton counts plus c|H_C gives H_C>=c, and deterministic insertion gives h>=d>=c^m. Thus for m>=3 the divisor c^(m-2) is greater than one and strictly smaller than h. This proves the prime-depth bound for all primes, including primes dividing c. The separate gcd(n,c)=1 consequence is immediate. Condition (4) is only necessary; for composite c it must not be rephrased using a single prime valuation.

For deletion consisting of an ordered triangle and k singleton components, ordered-join Hamilton counts multiply and labelled two-bin covers factor independently across components, giving H_C=3 and beta_C=12*2^k. The infinite bank of such rooted cores therefore qualifies. Although this bank has more primitive cores than the four K_j, each of its cores factors exactly into O/Z rooted-edge factors and one K_j: order the deletion SCCs, adjoin the same root with its original incidences to each SCC, and take their ordered root amalgamation. Singleton SCCs give O or Z, while the triangle SCC gives K_j up to cyclic relabeling. All cross arcs are forward by the SCC order, so the resulting graph is exactly the original core, not merely a count-equivalent graph. Associativity then expands any word of bank cores and spacers into the fixed K_j+O/Z grammar. Therefore the target67 hostile DOES apply to this infinite bank. If the bank contains all four K_j and permits arbitrary O/Z spacers, the generated languages are equal; for a subbank, its language is a subset.

For the fixed no-spacer triangle family, independent rotations give one incidence orbit for each j. Its four full coefficient profiles are correct. The sixteen two-core choices give {9,15,31,41,45,61,75,99}; depth one gives {3,5}, and deeper outputs are composite. Hence {3,5,31,41,61} is exact. Ordered joins of members cannot create a new prime because all Hamilton counts multiply and are positive. The one-core O1 graph has H=13 and correctly controls the no-spacer restriction.

## Separate TRIANGLE_67 completeness review

The common-content bound permits restricting to at most two cores only for the requested coprime-to-three target class. In that restricted grammar, h never decreases. Once h+x>N, each Z or j<=2 letter immediately exceeds N; O and j=3 do not lower h+x, so delaying the other letters cannot restore eligibility. O keeps h unchanged. A j=3 step makes h divisible by three, and another such step or O cannot undo that. The discarded suffix therefore cannot produce a new coprime-to-three target. This argument requires no bound on total length.

For expanded states, y<=d and d<=h follow from physical profiles, and x<=N-h, m<=2 give a finite exact-state space. Equal-state merging loses no continuation. The frozen normal replay confirms 91 expanded states and gaps {7,55,67}. The assertion that 67 is an omission of this grammar is proved; the prose assertion of realizability by other tournaments is not independently established by this audit. It is not needed for the hostile.

## Reproduction and independent controls

Both supplied manifests passed sha256sum -c. Normal and optimized replays of both supplied verifiers exactly match their respective frozen receipts. PYTHONDONTWRITEBYTECODE=1 was used to avoid touching frozen author inputs.

independent.py imports NONE of the author verifier/helper code. It independently builds all rooted labelled tournament adjacency matrices with root0, counts subset Hamilton paths using incoming-edge DP, computes labelled two-bin profiles by subset sums, and physically constructs ordered-root products. It checks all four output coordinates on 814 pairs: all input orders 1..3 against all non-singleton core orders 2..4. This includes both edge spacer cores and every four-vertex root incidence, and verifies every nonnegative physical coefficient. Independent triangle profiles match equation (5).

It also checks 44,000 steps of declared deterministic random algebraic controls for c=2..12, including composite c, arbitrary nonnegative response coefficients, and arbitrary spacers. These are filtration controls, not assertions that random coefficients are physical and not a substitute for the general proof.

A second, different finite inverse retains EVERY h<=67 state and saturates x at68 rather than applying h+x<=67. Every positive use of a saturated x makes h exceed67; zero use ignores it, so this quotient exactly preserves all outputs <=67. There are 229 reachable saturated states at core depth<=2, again giving precisely {7,55,67} as coprime-three odd gaps. This separately validates the author's 91-state cutoff result without copying that cutoff. Normal and optimized independent receipts agree.

Commands from workspace root:

    python3 session202-threshold-audit/content-independent/independent.py
    python3 -O session202-threshold-audit/content-independent/independent.py

The accompanying replay JSONs and input/artifact SHA256 files record the exact audit inputs and outputs. A reproducibility hardening recommendation: CONTENT-MANIFEST and TRIANGLE67-MANIFEST omit verify.py, although their verifiers import it. It was inspected and its hash is recorded in this audit's input manifest. This does not change the mathematical approval; any promoted reproducibility packet should pin that transitive executable dependency as well.

No files in the shared math checkout or frozen author packet were edited. No commit, push, or status promotion was performed.

## Scope clarification after parent review

The first audit version incorrectly warned that target67 could not be excluded for the infinite bank of cores whose deletions have one triangle SCC and ordered singleton SCCs. Parent identified the exact factorization above, and I independently verified it from the graph definition and associativity. The warning confused primitive-core variety with generated-language variety. It is corrected here; the frozen author claims, proofs, verifiers, receipts, and the two mathematical approvals are unchanged. This correction strengthens the stated obstruction scope: arbitrary compositions from that infinite bank, with arbitrary O/Z spacers, also cannot realize67. No claim is made for cores with a genuinely different deletion SCC type or a different cross-arc product.
