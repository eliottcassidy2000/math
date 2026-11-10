# Independent audit: exact arithmetic memory and synchronized residue images

Verdict: **PASS**. No mathematical correction required. The all-order statements follow from the reviewed proofs; the independent finite computations below corroborate them. Chan–Pak completeness and finite-bank density remain open.

Audited: `../arithmetic_states/PROOF.md` and `04-computation/experiments/session237_arithmetic_states/verify.py`. This audit did not import or rerun the producer implementation. The independent implementation is `independent_verify.py` beside this report.

## Proof review

### Reachability, including the division by two

Q is even because Q(2a)=Q(a) and Q(2a+1)-Q(a)=a+(a mod 2) is even. Write 2M=2^k u with u odd. A block of b ones contributes 2^b-2 to Q. Choose b>=max(k,2) divisible by the multiplicative order of 2 modulo u (no restriction when u=1). Then x=2^b-1 is divisible by u.

The claimed cross terms are correct. A low bit i and a high block beginning at B contribute x*2^(B-i-1); A>=k+bit_length(a) supplies the factor 2^k. Two high blocks with starts separated by D contribute x^2*2^(D-b), and D>=b+k. The stipulated gaps prevent additional adjacent-bit terms. Thus all cross terms vanish modulo 2M, not merely modulo M. That distinction permits dividing Q congruences by two to obtain q congruences modulo M.

The increment c/2=2^(b-1)-1 is odd and congruent to -1/2 modulo u, so is invertible modulo M. Varying the positive block count t fills all q residues while preserving m modulo 2M. Adding multiples of M to t proves arbitrarily large realizers. This remains valid when M=1, when u=1, and when M has both odd and even factors.

### Minimality is an output-automaton theorem, not an input-state count

The proof distinguishes any two reachable states by one of the suffixes ε, 0, 1. With equal current outputs, suffix 0 detects different r modulo M. The only remaining distinct pair has r values differing by M and equal z. Suffix 1 separates them by floor(M/2) or ceil(M/2), nonzero modulo M for M>=2. Therefore 2M² is genuinely minimal, not just the size of a convenient reachable carrier.

At M=1 there are two reachable arithmetic states, but their outputs and all future outputs agree. The minimal output automaton has one state. The statement correctly handles this exception. Optional leading zeroes fix the initial state and do not change the function; restricting to canonical nonempty-leading-one numeral representatives would be a different language/model question, not a contradiction.

### Joint seed image and normalization

Write m=2j+b and Q=2q. Then exactly

    (H_s(m)-1)/2 = (h-1+(R-e)b)/2 + Rj-dq.

Physical seed parity h odd and R-e even makes the first term integral; no inverse of 2 modulo an even M is used. The proof divides an integer first and only then reduces modulo M. Reachability of (m mod 2M,q mod M) makes j and q independently variable modulo M for either fixed b. Consequently the two affine module cosets are the full image, not just a necessary upper bound.

For the first example, Y_s=2+3Y_t-b is correct. For the second, det([[2,-1],[11,-5]])=1 proves full pair image within each fixed parity, over every Z/M including composite M. The synchronized-parameter assertion is expressly different from independent choices of parameters for different seeds. No integer coverage or finite-bank density conclusion follows.

## Independent exact checks

The auditor computed Q from bit-pair distances and adjacency, rather than its recurrence. The recurrence was checked only afterward. Pumped integer witnesses were then built explicitly for every pair (r,z) at all M=1,...,17, and their q values reevaluated from their literal binary digits. This tests sufficiency independently of automaton traversal.

A separate graph-closure and full output-partition-refinement implementation checked moduli

    1,2,3,4,5,6,8,9,10,12,15,16,18,21,25,27,30,32,35.

It recovered 2M² reachable states for each modulus, one minimal state for M=1, and 2M² minimal states for all others. One refinement using output and both successor outputs already separates all states for M>=2. This directly tests the suffix-length bound rather than assuming it.

For every M=1,...,17 the explicit pumped integers generated the same joint triple image as the affine formula. The root/rank-one projection equals exactly the two claimed affine lines. The root/unimodular projection has M² elements separately for each parity b.

The auditor also exhaustively searched labelled tournaments to find physical realizers of all three example seeds. End-counts, deletion counts, and the one-bridge response were computed from permutation enumeration, not the producer's endpoint DP. For 32 root tails and eight tails of each other seed, literal tournaments were constructed in physical low-bit-first order and Hamiltonian counts were independently enumerated by permutations. Their exact H counts and normalized residue formulas agree. The seed adjacency rows and all test counts are in `results.json`.

The implementation uses explicit exceptions, so its checks remain active under optimized Python. Matching normal/optimized results and hashes are recorded alongside this report after execution.

## Scope and provenance cautions

- The producer's pre-audit standalone checks are internal corroboration; this separate implementation supplies the independent audit.
- The historical binary tail law and pumping construction are credited rather than claimed anew. This audit checks the displayed mathematics and does not establish external publication priority.
- Physical seed verification is finite and explicit. The theorem's arbitrary physical seed list is conditional only on the known exact binary tail law and stated parity, with its new residue-image claim proved here.
- No correction to the main theorem, suffix theorem, examples, or normalization was necessary.

Reproduce:

    python3 05-knowledge/results/session237_chanpak_synthesis_20261109/arithmetic_states_audit/independent_verify.py
