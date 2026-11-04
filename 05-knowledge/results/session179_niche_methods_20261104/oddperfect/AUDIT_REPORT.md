# Independent audit: primitive-order domains

Verdict: **PASS**, with no mathematical correction required. Audited PROOF.md, check_primitive_domains.py, and RESULTS.json directly. No primary code was imported into the independent checker. No installs or prohibited ledger access.

## Proof checks

- Bang–Zsigmondy correctly applies to p^d-1 with d>1 and odd p. The positive-base exceptions are precisely d=2 with p+1 a power of two, and base 2/exponent 6. The first is excluded for ordinary sources by odd n, and for the Euler source by q=1 mod4; the second is excluded by odd support. An exceptional primitive divisor cannot be 2, because 2 divides p-1. All uses of ord_r(p)=d and r|sigma(p^(n-1)) are valid.
- The top primitive level gives n<=max(S)-1 and a<=max(S)-2. Uniqueness of multiplicative order gives tau(n)-1 distinct target primes, hence tau(n)<=s. No unsupported uniform-in-support bound is claimed.
- PC really must hold at every divisor, not merely at the top level. Both concrete factorization counterexamples check exactly.
- The exact row stripping is complete: no general factorization oracle is needed. Role-compatible candidates have positive exponents. Closure plus the 2-adic condition and all odd-prime column equalities is equivalent to sigma(N)=2N. The (s-1)^s count bounds row combinations, not bit complexity, correctly stated.
- The divisor-ideal / multicomplex interpretation is sound. In the phrase 'role-compatible lengths', retain the original n>1 condition; the document already defines it that way. No topological invariant is improperly promoted to an obstruction.
- The stronger primitive-domain layer does reject both Euler branches for S={3,5,7,11,13}, while the stated maximal-order graph is strong. This is a mechanism illustration, not a new prime-count theorem.

## Independent verification

`INDEPENDENT_AUDIT.py` uses multiplicative orders from divisor testing of r-1 with modular exponentiation, direct geometric sums, its own prime stripping, independent Boolean-Warshall reachability, and actual trial factorization of 25 small primitive levels. It reproduces:

- 189 support/Euler cases;
- 124 closed rows, all preserved by independently recomputed domains;
- every concrete source domain in RESULTS.json;
- strongness of both concrete maximal graphs;
- the n=9 and n=5 boundary factorizations;
- the d=2 exceptional cases for bases 3 and 7 and the absence of other tested exceptions.

Primary normal and optimized runs agree byte-for-byte with each other and RESULTS.json. Independent normal and optimized runs also agree byte-for-byte.

Primary SHA256: 6c83a01848c189a482c7d7f6f905f8b889a8ebe48690400d4a911c7d6338a17a.
RESULTS SHA256: 098a369c27722592bf39af4d62173d45c833118b3b9ac5864c09758010a7c58c.

## Source attribution

Opened the cited primary-host PDF directly: https://math.colgate.edu/~integers/k27/k27.pdf . Its Theorem 13 on printed page 325 (PDF page 7) states the primitive-divisor result with the two positive-base exceptions used here. The additional negative-base (2,-1,3) exception is irrelevant because b=1. The article title, authors, volume, and year match the proof. Historical priority is appropriately not claimed; this is an exception-safe integration of a classical theorem.

## Scope

This is a rigorous support-specific finite decision and sparse-domain result. It does not establish that no odd perfect number exists or provide a uniform finite support universe. Tests do not replace the cited classical theorem. All these boundaries are already respected in PROOF.md.
