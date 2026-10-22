# Independent audit: one-apex continuation observer

Verdict: PASS. Frozen 2026-10-22. The mathematical arguments establish the stated all-word and all-modulus theorems. No substantive correction is required. This covers the final PROOF.md, including its CRT and strong-completion additions, and the inherited Session25 endpoint-counting recurrence. This is not a literature-priority assessment.

Audited source: `/workspace/shared/session26-tournament/PROOF.md`.
SHA-256: `6ce3a12477ffeb80f04dbaf4b40219b6a00bce8fc1e8ab73bf215bdb44f9faf0`.
The source was reread after executor recovery; this is also the hash recorded in the completed independent results. The previously suggested strong-completion wording refinement is present.

## Definitions and recurrence

The carrier is the free binary-word monoid, including the empty word. Vertices are an ordered transitive core and a marked apex; a bit 1 means apex-to-core orientation. H counts directed Hamiltonian paths without quotienting by automorphisms. The empty word has single-vertex path count 1. Session25's internal-apex endpoint split retains all orphan paths and gives the stated recurrence, including the previous-bit term on append 0. The state is relative to right concatenation, not arbitrary insertion, substitution, or metagraph wiring.

## Reachability and deterministic minimality

The reset proof is valid. For M=2^k m choose L at least max(1,k) with 2^L=1 modulo odd m. The reset block 1^L 0^L returns A to zero modulo M and p literally to zero, and advances H by D=(2^L-1)^2+1 modulo M. Its gcd with M is gcd(2,M). Positive block counts already sweep the whole permitted H fiber. A binary suffix with value congruent to a modulo M and last bit p exists by CRT for odd M, and exactly under the stated parity compatibility for even M. Its H increment is independent of the reset H. This proves every proposed state reachable, not merely every output residue attainable.

The probes empty, 0, and 10 separate all states when M>=3. Equal probes imply delta H=0, delta A+delta p=0, and 2 delta A=0. Therefore 2 delta p=0 modulo M. Since delta p is -1, 0, or 1, M>=3 forces equality of the bits, then equality of A. For M=1 and M=2 all outputs are constant. Raw-state counts and quotient-state counts are correctly distinguished. The minimal counts are 2M^2 for odd M>=3, M^2/2 for even M>=4, and one for M=1,2. This is a deterministic residue-output state count, not an unrestricted computational complexity theorem.

## Weighted dimension and integral fingerprint

Four homogeneous coordinates give a field-linear realization. The displayed minor has determinant -16, so rank is exactly four over Q and every odd-characteristic field. In characteristic two every entry is one, so rank is exactly one.

The full infinite Hankel determinantal divisors are 1, 2, 4, 16, not merely those of an isolated sampled minor. Odd entries force universal 2- and 3-minor divisibilities. State-row parity gives a factor four for every four-row determinant; suffix-column parity gives another factor four. Their product forces every four-minor divisible by 16. The displayed minor attains all required gcds, and larger minors vanish by the realization. Thus Smith factors 1, 2, 2, 4 are valid under the stated interpretation through gcds of finite minors.

An independent stronger lattice certificate was checked. With U=(empty,1,10,11), the state matrix R has rows (1,0,0,1),(1,1,1,1),(3,2,0,1),(1,3,1,1). With V=(empty,0,10,00), the response matrix C has rows (1,1,1,1),(0,1,2,3),(0,1,0,1),(0,0,2,0). Both are integral invariant bases for their respective reachable and observable lattices. The independent script verifies exact integer transition identities and R C equals the displayed minor. No floating-point rank inference is used.

## CRT and strong completions

For coprime m,n>=3, coordinate reduction gives precisely the fiber product over literal p. CRT lifts h and a uniquely, and the parity constraint at an even factor survives. Each p fiber is half the local state set, so N(mn)=N(m)N(n)/2. Modulus two loses p observably, and N(2m)=N(m) for odd m. This is a genuine exception to the two-nontrivial-factor formula.

Every raw state has a witness beginning with 1. Suffixes 0,10,110 end such witnesses with 0 and hence produce strongly connected tournaments. Their readouts H+A+p, H+2A+2, H+4A+4 yield the same separating equations. The state-count lower bound therefore survives these restricted queries. The final proof correctly avoids claiming every intermediate prefix is strong.

## Independent exact controls and artifact verification

`independent_verify.py` imports no producer code and uses no optimization-sensitive assertions. Reproduce with `python3 independent_verify.py` or `python3 -O independent_verify.py`. The completed `results.json` reports PASS and the audited source hash above. Its checks cover:

- All 1,023 words of lengths 0 through 9 against explicit adjacency and forward subset dynamic programming.
- Exact reachable sets and independent Moore partition refinement for M=1 through 100.
- Constructive reset-and-suffix reachability for M=1 through 30.
- Strong-completion probe separation for M=1 through 100.
- CRT fiber-set equality for all coprime pairs 3<=m<n<=12.
- Exact minor gcds and integer invariant-lattice identities.

The executor interruption left the original prose-write attempt unexecuted. After recovery its wrapper no longer existed and filesystem inspection confirmed REVIEW.md absent. This review was then written once, while the already-completed script/results were retained unchanged. No repository or Git files were changed by this reviewer.

Scope exclusions remain essential: no global Hamiltonian spectrum classification, prime-value attainment, arbitrary wiring completeness, LRC implication, or literature novelty is proved. Separate cut/extremal research and a possible modular row-span cardinality extension are outside this frozen audit.
