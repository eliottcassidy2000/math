# Independent audit: infinite bounded-degree hole atoms

Verdict: **PASS**, with no mathematical defect found in NOTE.md or verify.py as inspected. This audit does not certify historical novelty or a major open-problem result.

## Mathematical checks

- Definitions match the supplied THM-2500 source snapshot: retained composite weak atoms; strict atoms additionally require nonempty ambient strict factor fibre, excluding prime squares.
- Weak cutoff 2K+1: a rank at least 2K+2 splits into two retained factors of ranks at least K+1.
- Strict cutoff 2K+2: rank at least 2K+3 has a split with unequal ranks and hence unequal factors.
- Weak equality: each proper split at total rank 2K+1 has exactly one factor of rank at most K; the target itself is automatically retained.
- Strict equality: if the multiset contains distinct prime entries, a balanced-rank block and its one-entry swap cannot both equal the square root. Therefore only pure prime powers can attain the bound; the remaining strict splits require precisely p through p^K to be holes.
- Rank-ideal formula correctly removes prime squares only when K=1 and adds pure powers at rank 2K+2. The union is disjoint by rank.
- Arbitrary prime-hole semiprimes count as Q(P-Q)+(Q^2+Q(2s))/2; triples count by the complete symmetric polynomial. The strict correction is exactly -Q(2s)+Q(4s). No multiplicity is confused with support membership.
- For real s down to 1, P(s)=L+O(1); all P(rs), r>=2, are bounded. The top complete symmetric term is P(s)^j/j!, and nonidentity cycle types have at most j-2 singleton cycles. This proves the stated leading logarithmic orders.
- Under the explicitly conditional Q(s)~c L^alpha, the only dominant competitors are QP and Q^3/6. Their crossover is alpha=1/2; the alpha=0 and alpha=1 endpoints are valid, and strict corrections are bounded.
- Nonempty prime-hole sets have divergent reciprocal atom sums by the fixed-p semiprime ray, including when Q(1) is finite. Absolute convergence for real part >1 follows directly since the sets are subsets of the positive integers.
- The THM-2500 snapshot really gives finite Dirichlet-polynomial coefficients times a prime tail plus finite remainders. Such a finite-hole expression has at most first-order logarithmic growth; the comparison is justified.

K should explicitly be called a positive integer. The phrase “linear-pole heuristic” is better replaced by “first-order logarithmic heuristic.” These clarity suggestions were sent to the author; neither affects the arguments under their intended meaning.

## Reproduction

Producer normal and optimized runs were independently rerun and byte-match each other and its stored normal.json. Counts: 159936 prime-hole target cases, 12495 rank-ideal target cases, 204288 bounded-rank target cases, 7196 exact series coefficients, 30 sharp endpoints, six hostiles.

The separate verify_independent.py imports no producer code and uses prime-exponent lattices and all proper exponent-vector pairs. It checks 20340 independent bounded-rank hole cases and 12 members of the family H={2^a:a>=1}, n=2^j*3, which explicitly witnesses why bounded DISTINCT-prime degree cannot replace Omega. Its normal and optimized outputs byte-match.

The independent finite lattice coverage is intentionally limited; the proof review establishes the unbounded claims. No installs, forbidden files, source edits, or external writes were performed by this auditor.
