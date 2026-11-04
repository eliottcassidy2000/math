# Independent audit: periodic-template parser and rank

**Verdict: PASS. No mathematical correction required.**

Reviewed `RESULT.md` and `periodic_template_rank.py`, independently rederived the formulas, and ran both producer and independent code under normal Python and `python -O`. Producer reruns are byte-identical to one another and `test_results.json`. Independent outputs also match byte for byte.

## Mathematical checks

- The cyclic affine carry identity is correct. Since every source is odd, D and B_i are odd, nonzero z has K≥1; thus division of K−1 in the parser is legitimate.
- The strict threshold is exactly S_j<K, not S_j≤K. At equality of current precision and the expected letter, the two odd terms cancel to yield a strictly higher valuation. Positive sources ensure 3n+1 is nonzero, so the mismatch valuation is finite.
- D cannot vanish because 2^A and 3^L are distinct positive prime powers. Both D>0 and D<0 are permitted. Signed cofactors cause no problem because the energy uses absolute value.
- The phase weights compensate the multiplier by precisely 3^(L a_i−A), including the cyclic wrap. Energy is a positive integer and bounds |Dn−B_i|. Thus each phase has finite energy sublevel sets; the finitely many phases also have finite joint sublevels. Matching steps strictly decrease the lexicographic pair.
- For z=0 at a positive odd source, the carry identity forces legal odd steps through all phases and returns after the word. The text correctly treats this as a cycle case, not an assumed route to 1.
- Repeating a word r times gives the stated odd H. More explicitly, on phase zero, E_repeat=H^(rA) E_original^r. This verifies the claimed equivalence of the phase-zero ordering.
- The neutral inverse tower is integral under the specified ternary condition, has genuine valuation-two steps, and exactly preserves the old energy. The density statement is in the original t parameter: 65536 is invertible modulo 3^r. It adds no source coverage, as correctly stated.
- No global rank or global template selector is established, and no inference from finite long shadows to Collatz convergence is made.

## Independent controls

`INDEPENDENT_AUDIT.py` imports no producer functions. It computes affine coefficients by rational function composition, constructs exact v2(Dn−B)=K source cylinders directly, and compares the cumulative-cost rule with actual U trajectories. It checks:

- 3,060 exact-precision cases across 340 words of lengths 1–4, letters 1–4;
- 317 positive-D words and 23 negative-D words;
- precision thresholds K∈{1,2,3,4,5,7,11,19,31}, phase transport, mismatch direction, and properness;
- powers r=2 and r=3, including the exact energy multiplier H^(rA);
- five exact long boundaries of 205, 205, 615, 820, and 820 matched letters;
- center 1 for repeated valuation-two words, a negative signed-cofactor case, and 11 neutral inverse towers with the old family energy inequality.

The producer's 1,364-word, 51,702-source, 145,920-check and five fixed-point counts reproduce exactly. Its longer tests are separate from the reported short-test maximum of 11, so that field is not a contradiction with 803 matched letters in the long controls.

## Provenance checks

Verified that snapshot `1b91742ce55a202b73c80d56561ba56d00b503a6` exists locally. Read the cited branch-toll rank's energy definition and larger-child family, the weak-reset cycle-center discussion, and run-block inverse-phase formula. The source descriptions agree with these texts. Also checked the bounded inverse-cover, guarded-affine, atomic-rank, backward-measurement, and polynomial-dual status statements, and THM-4493/4514/4523/4507 headers and stated scopes. They support the brief source descriptions without supplying the missing global theorem.

No forbidden-path reads or package installations were performed. All audit claims are restricted to the specified templates and theorem hypotheses.
