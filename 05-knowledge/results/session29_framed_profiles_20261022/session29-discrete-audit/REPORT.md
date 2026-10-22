# Independent Session29 discrete audit

Verdict: PASS on the stated Collatz and tournament-profile results, with one minor domain clarification required in the debt paragraph: assume lambda1, lambda2 > 0. No substantive counterexample found. This audit neither proves Collatz convergence nor complete tournament metagraph wiring.

## Method and scope

Read both producer reports and their scripts; read the root seed checker and repository THM-2183's uncrossing proof. Copied scripts into this audit directory before running, so no producer artifacts or repository files were modified. Exact audited-source SHA256 values appear in audited-source-sha256.txt.

Replayed Collatz verify.py normally and with -O. Both outputs were byte-identical to the producer's verification.json. Replayed root seed checker and compared its JSON byte-for-byte with the root output. Replayed profile second_loop.py and mirror_check.py normally and with -O; their imported audit.py also reran the complete n<=6 orbit checks. Result, second-result and mirror-result JSON files matched the producer files byte-for-byte.

Fresh independent.py imports no producer modules. It exhaustively finds each valuation-word residue using literal trajectories for all 511 positive words with summed clock 1 through 9, then independently enumerates modular powers rather than using the producer's binary lift. It verifies the oddpart bijection, cylinder formula, one exponent class, an explicit root-free positive lift, and exactly three distinct mod105 sidecars for each word. It checks annulus cardinality for integer N=1 through 10000. It independently uses adjacency matrices and direct permutations for the five-vertex seeds, all proper decks, Hamiltonian paths, automorphisms, anti-automorphisms, path orbits and reflection. Unlike a mere fixed-point count comparison, it verifies that every anti-isomorphism induces the same full reflection map on the tile set. It checks merged class/mass/blue coefficients through r=10 and all 16 allowed local binary assignments in the join-metric uncrossing inequality. Normal and optimized outputs are identical.

Commands (run inside this directory):

    python collatz_verify.py > collatz-replay.json
    python -O collatz_verify.py > collatz-optimized.json
    python root_seed.py > root-seed-replay.json
    python second_loop.py > second-replay.log
    python -O second_loop.py > second-optimized.log
    python mirror_check.py > mirror-replay.log
    python -O mirror_check.py > mirror-optimized.log
    python independent.py > independent.json
    python -O independent.py > independent-optimized.json
    cmp independent.json independent-optimized.json

## Proof audit: Collatz

The exact cylinder modulus is 2^(A+1), not 2^A. The extra terminal odd bit is necessary. The power-of-three subgroup and its order yield the stated admissible initial prefixes. The odd-index oddpart family covers all odd residues, giving every finite word, with conditional exponent density 2^(-A). The density is on exponents, not on positive integer sources, and cannot be used as source ownership or recycled credit.

The first-reset descent split and arbitrarily long growth prefixes are valid. The ROOT exclusion is valid because every finite prefix endpoint is injective affine in the initial source. At most finitely many members hit 1 within a fixed prefix, so removing them leaves infinitely many realizations. No assertion of an infinite positive orbit with arbitrarily prescribed infinite word follows.

The mod15/mod35 register is correctly joint: the word fixes k modulo 4, allowing exactly three of the six k modulo 12 states. Their source values modulo 105 are distinct. The division by 4 must occur after retaining modulus 420, as the report does.

The annulus [4N,256N) contains at most two native oddpart-family sources for N>=1, since three eligible consecutive terms have endpoint ratio strictly greater than 64. This is a correct obstruction to a naive cardinality/mass-preserving inference, not a theorem excluding a future weighted descendant/ownership construction. Uniform 153N coverage and the remaining R' estimate remain open.

Minor repair: state positive lambda explicitly in Section 6. Without positivity, the claimed inequality composition/weighted cone closure does not follow. With positive lambda, the algebra and closure are correct.

## Proof audit: tournament profiles

The seed statistics are independently confirmed: mask40 has H=15, Aut order3, 5 path-orbit tiles, 1 reflection-fixed tile; mask41 has H=13, Aut order1, 13 tiles, 3 fixed tiles. Both are self-converse and have identical proper decks and scores. The former rigid/rigid description does not survive in the audited report.

The SCC argument for ordered words, equal D<=4 profiles, product H and Aut, and reversal-as-converse are sound. The metric conclusion follows from THM-2183's valid equal-factor-order strict uncrossing proof. It concerns unrestricted unlabeled arc-reversal distance, not fixed-spine complement edges.

Weighted Burnside is used correctly: 394=15^2+13^2 for Hamiltonian weights and 194=5^2+13^2 for one-representative tile weights. Actual merged tile mass is instead 18^r because nonpalindromic nodes contain both classes. The symmetric incidence GF (1+4z)/(1-18z^2) follows from forced matching of paired blocks and one middle fixed tile. No nonempty family member is all-blue. These are incidences, not half as many internal edges; complement partners may leave the family. Empty-word r=0 is a formal GF convention.

Downward deck identity, score-binomial inversion, and bounded-union overlap identities are valid with unmarked motif data and n known. They do not imply low profiles refine scores or that arbitrary global markings are encoded by low-order decks.

## Literature check and limits

Independently opened Stockmeyer's primary PDF https://arxiv.org/pdf/2110.06158 and located Theorem 2 on PDF page index12: its nonisomorphic same-deck tournament construction is indeed for power-of-two orders at least4. The stated unmerged separating-order consequence follows. No merged-class conclusion is imported from that construction. The initial HTML route returned404; the PDF succeeded. This audit does not independently reconstruct Stockmeyer's full matrices, inspect the separate 2026 score-sequence paper, or claim historical novelty.

Finite tests audit implementations only. All-order conclusions above rely on the displayed proofs. No new Git action, dependency installation, broad repository scan, or hidden resource access was used.
