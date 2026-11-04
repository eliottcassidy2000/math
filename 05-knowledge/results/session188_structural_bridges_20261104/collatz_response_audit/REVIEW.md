# Independent review of the four-letter reset compiler

**PASS for the mathematical claims and their stated scopes.** The earlier implementation concern was corrected and rechecked before this final review. This review includes the new rational-address theorem and the exponential-limsup corollary; the latter is a proof review, not a claim that finite computations establish its limsup.

## Proof review

- **Exact phase.** The carry sum is correct. Integrality is equivalent to one exponent class modulo 2*3^L because 4 generates the principal units modulo 3^(L+1). The last carry summand is a unit modulo 3, excluding phase 0 modulo 6. The legal phases are exactly 2 or 4 modulo 6.
- **Positivity and exact valuations.** Starting at m_b=(2^b-1)/3, the successive divisibility conditions obtained from the final congruence make every inverse step integral. An integral inverse of a positive odd input is again positive and odd. Thus every displayed letter is an actual valuation, not merely a formal affine exponent. For b=2 the route has reached 1 before its final reset; for b>=4 an earlier visit to 1 is impossible.
- **Alphabet and sharpness.** The two legal choices differ by two, and E_(a+2)=4E_a+1 prevents both predecessors from being multiples of 3. Consuming one ternary digit per step preserves legality. In the forbidden phase 14 modulo 18, the endpoint is 7 modulo 9, and alphabet {1,2,3} forces exponent 2 into a multiple of 3. This proves the sharp initial-alphabet threshold without a search assumption.
- **Extremal costs.** The formulas for all-1, final-3 and all-2 prefixes follow by their exact carry sums. The maximum first-hit exponent is 2*3^L+2, attained by all-2 prefixes. It is crucial that this is the one-reset problem; allowing several extra odd steps is a different optimization.
- **Periodic address.** D is a unit modulo 3 and B is a unit modulo 3. The exact repeated-carry identity gives approximation error of valuation rL+1, so the phase congruences are compatible. The binary mismatch formula is valid when reset b differs from a_0. The equal case is properly excluded.
- **Rational-address theorem.** Write j=p/q in lowest terms, q>0. Because j is a 3-adic integer, 3 does not divide q. The continuous exponent law gives C^q=4^p as an equality of rationals. Ordinary prime valuations imply C=+/-2^k, and the cyclic carry identity fixes k=v_2(C)=a_0>0. In particular p cannot be negative; also 2p=qk, so q divides 2. The congruence C=1 modulo 3 fixes C=(-2)^k. The center rho=(C-1)/3 is an odd integer. Cyclic carry identities make every rotated center an odd integer and realize the signed periodic orbit exactly. Its next value is +/-1, so periodicity forces the fixed orbit. Only all-1 and all-2 templates remain. The negative square-root branch is essential: 4^(1/2)=-2 in 1+3 Z_3, whereas the all-2 address is 1. No real-logarithm branch is used.
- **Digit nonperiodicity and limsup.** A rational 3-adic integer has eventually periodic base-three digits, and conversely. A nonfixed template therefore has infinitely many nonzero blocks of width L in its address. At each corresponding index r, the least phase is at least 2*3^((r-1)L), while every least first-hit reset is at most 2*3^(rL)+2. These bounds give limsup log_3(b_r)/(rL)=1. The all-1 and all-2 cases have the same limsup by their explicit formulas. This proves neither a uniform lower bound at every r nor positive density of large reset indices.

The scope statements are appropriate: the source is constructed after the word, not independently supplied. The results do not prove Collatz convergence, a natural-density bound, absence of cycles, digit equidistribution, transcendence, or descent across template changes.

## Implementation and complexity qualification

`reset_phase` and `prefix_for_phase` genuinely work with moduli having O(L) ternary digits; the digit lifting has L stages, not O(L) bit-operation cost. Modular exponentiation and integer arithmetic still have their usual nonconstant costs, including the encoded input-letter size.

The initially inspected `completion(expand=False)` and `periodic_phase` constructed exact carries. The author corrected both before freeze, and this final review rechecks the edits. `completion(expand=False)` now retains only a modular carry and symbolic formulas; `periodic_phase` computes both carry and denominator modulo 3^(rL+1). Exact B and the source are expanded only under `expand=True`. Huge encoded-letter controls exercise both repaired paths. The proof now explicitly includes the input bit length and cumulative-cost sidecars, so its compactness claim is accurate.

## Independent checks

`audit.py` uses a direct carry sum rather than the author's recurrence, enumerates each small exponent class independently, builds sources by inverse iteration, and replays every step with literal arithmetic. It checks:

- 1,554 words over {1,...,6} of lengths 1 through 4;
- 3,108 complete literal routes, including 45 non-first-hit b=2 controls;
- all 6,560 allowed phases through depth eight;
- 1,093 forbidden phase classes for the alphabet-three obstruction;
- 60 repeated-template center identities;
- 9,330 words for the rational-address arithmetic classification;
- huge encoded-letter phase controls without materializing their exact carries.

Normal and optimized receipts agree byte for byte. The inspected proof and verifier SHA256 hashes are retained. Finite checks support the implementations; the unbounded conclusions rely on the proofs reviewed above.
