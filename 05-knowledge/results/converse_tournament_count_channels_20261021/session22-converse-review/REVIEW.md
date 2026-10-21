# Independent converse-channel audit

Verdict: PASS. No blocking mathematical or implementation issue found. No repository or remote mutations performed.

## Scope and findings

- Reviewed `session22-converse-channels/README.md`, `channels.py`, `check.py`, frozen `incoming_local_scale.py`, and all four provenance sources.
- The raw carrier for F_n(a,b) is correctly simultaneous isomorphism of a ordered tournaments on one vertex set with 2^b vertex colors. Its Burnside fixed-point product is valid for literal integers a>=1, b>=0. The document correctly does not assign a new object-count interpretation to arbitrary signed core substitutions or no-fixed restrictions.
- Degree-local transfer retains the correct odd-factorial scaling. The response remains unweighted, and the code multiplies the actual cross contribution by a exactly once. The b contribution counts actual cycles, not substitution excess. The nonnegative-shift and integral-transfer arguments are valid; inherited carry width and remaining-degree masks are unaffected.
- The anti-automorphism monodromy proof is correct. The centralizer factor 2^K cancels the additional 2^K orientation factor, proving C_(2m)=F_m(2,0) and C_(2m+1)=F_m(2,1), including orders zero and one.
- Unique ordered SCC decomposition proves A=1/(1-S), while converse reverses and converses the sequence, proving C=A(z^2)(1+U). The code's coefficient recurrences implement these ordinary generating-function identities correctly. Empty and singleton conventions are consistent. Both converse-merged averages are correct.
- Strict literal integer and Boolean typing holds in public channel_dp; counts delegates its bound validation before using it.
- This is newly supplied reproducibility, not recovery of historical implementation bytes, a new classical count formula, or a measured acceleration claim.

## Canonical explanatory correction

The actual THM-283 clause 3 incorrectly says a cycle divisible by four creates an odd pair orbit, and allowed cycles have only even pair orbits. Its partition filter, formula, and sequence remain correct.

A minimal replacement for clause 3 is:

“3. **Parts ≡ 2 mod 4 (not ≡ 0 mod 4):** For an unordered-pair orbit of length L, consistency requires g^L to preserve the two endpoints individually when L is even and to swap them when L is odd. In an even vertex cycle of length r, its antipodal pair orbit has length r/2 and swaps endpoints, so r/2 must be odd. Thus r must be 2 mod 4; non-antipodal internal pair orbits then have even length and preserve endpoints.”

This agrees with the independently pinned HYP-2078 monodromy expression `swap XOR (L mod 2)`. Direct counterexamples: the transposition (01) has an allowed unordered-pair orbit of length one and fixes both labelled tournaments under anti-action; the 4-cycle (0123) has a forbidden antipodal orbit of length two and fixes zero tournaments under anti-action. Endpoint swapping, not unordered-pair length alone, decides consistency.

## Executed controls and receipts

- `python session22-converse-review/audit.py` and `python -O session22-converse-review/audit.py`: 114 checks each; JSON receipts identical.
- This independent checker enumerates every labelled tournament bitmask and every permutation for n=0..5, computes the literal relabelled mask, and directly tests equality or equality to its bitwise complement. It also checks all five raw multi-tournament/color channels via powers of literal fixed-mask counts and literal permutation cycle counts. It does not use the producer's partition or edge-orbit consistency logic for those comparisons.
- Checks all four source pins against current bytes, both explicit permutation counterexamples, both generating-function convolutions through n=14, and rejection of expanded invalid integer/mode inputs.
- Fresh producer checker executions: 578 checks each in normal and -O; receipts identical. These cover all four core/no-fixed modes, five parameter channels, all bounds 0..20, frozen-kernel coefficient and telemetry equality, full edge-action Burnside through n=7, literal isomorphism class and reachability census through n=6, the published C prefix through n=19, and all six vectors through n=40.

Artifacts: `audit.py`, `audit-normal.json`, `audit-optimized.json`, `producer-normal.json`, `producer-optimized.json`. Audit receipts pin the exact reviewed source bytes. No large census or performance benchmark was performed.
