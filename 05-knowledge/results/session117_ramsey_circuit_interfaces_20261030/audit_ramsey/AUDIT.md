# Independent audit: twin projection and alternating skew doubling

Verdict: **PASS**, subject to the stated restricted scope (no circuit lower bound, Ramsey-number record, or algorithmic priority claim).

Reviewed `ramsey/ramsey_twin_projection.md` against the exact inherited arc convention in `THM-483-zigzag-law-refutes-the-sandwich.md`. No producer implementation was imported by the independent checker.

## Uniform proof audit

1. A doubled chain with twin v must have all other primes strictly between v0 and v1, and all other plains outside. For a second potential twin x, an outside plain x0 forces an orientation against the corresponding inside prime x1. Hence there is at most one twin. All before-plain vertices dominate v and the inside primes in the base; all inside primes dominate all after-plain vertices; v dominates the primes. Reversing the prime block therefore yields a valid base transitive order.
2. The projected base order is unique. Its twin occupies position p; all primes occupy precisely a nonempty interval p..q. Every interval yields a valid chain, so the number of lifts is k(k+1)/2, including k=1. No lift exists for empty or nontransitive support.
3. The derivative identity counts those lifts once with size k+1; the empty subset remains solely in Z. The n-bound requires a nontransitive base: transitive bases attain n+1. The tighter minimum with 2t is valid since each copy part separately has size at most t.
4. Enumerating absent/plain/prime states has exactly 3^n states; testing each by triangle or degree checks is polynomial work, with polynomial auxiliary storage. This improves only the explicitly named naive 4^n scan. It does not establish a best-known bound or polynomial optimality certificates.
5. For A_l, excluding U positions straddled by W extremes proves t=l+1. The base is cyclic already at l=1, so t(D)<=2l+1 uniformly. The displayed full-support chain attains equality.
6. For l>=2 a twin chain has at most l+2<2l+1 vertices, so maxima have full support without twins. Each native cyclic triple requires exactly two primes. The j=1 and j=l equations force all W bits equal through an overlapping index; then all W are prime and all U plain. At l=1 this overlap is absent: three full-support assignments plus nine twin-containing lifts give twelve maxima.
7. The base polynomial counts zero W, one W, and pairs of extreme W indices. For fixed distinct extremes, the optional interior W positions and allowed outside U positions total l-1. All resulting subsets are transitive, so the coefficient formula is exact for every l>=1.
8. The Boolean clause encoding is exact: native cyclic triples permit only weight-two assignments; native transitive triples forbid only prime endpoints with a plain middle. At most five width-three clauses per triple suffice. This is an upper-bound construction, not a circuit lower bound. The alternating family for l>=2 is a single-pattern indicator.

## Historical quantifiers

THM-483's corollary wrote delta(A_l)=l for general l, but its A2 paragraph only supplied the lower bound and finite equalities through l=5. Its explicit Honesty section retained the uniform bracket 2l+1 <= t(D(A_l)) <= 2l+2. The new uniform proof genuinely closes that documented gap, while not refuting the historical unboundedness claim (which needed only the lower bound). It should be described as a repository proof-gap closure, not a literature novelty claim.

## Independent executable evidence

`check_independent.py` builds the double directly, tests transitivity by absence of cyclic triangles, scans all doubled subsets, and verifies each support's lift multiplicity. It also checks the full-support Boolean predicate without importing any producer code. All 1,099 labeled tournaments of orders 1..5 are covered. A_l base polynomials and doubled maximum multiplicities are exhaustively checked through l=4; full-support equations are exhausted through l=8, and every surviving assignment is checked in the actual doubled tournament.

`normal.json` and `optimized.json` are byte-identical successful receipts from ordinary Python and `python -O`. Explicit exception gates remain active under optimization. Uniform claims rely on the proofs above, not on these finite checks.
