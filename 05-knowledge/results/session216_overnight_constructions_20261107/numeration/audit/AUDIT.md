# Independent numeration audit

Verdict: PASS for the mathematical claims and retained finite receipts, with the implementation-domain qualification below. No all-target tail theorem follows.

Inputs are frozen in frozen/ and individually identified by INPUT_SHA256SUMS. The core state file is byte-identical to session216-three-core-audit/frozen/CORE_STATES.txt. The separately audited construction formula is the dependency; this audit does not independently re-prove that whole packet.

## Analytic review

For positive P, recurrence induction proves 0 <= f(P) <= P-1, with equality on the right exactly for positive binary repunits. The even branch is strictly below P-1 except the separate P=1 base; the odd branch requires an odd half-argument and equality there. Substitution yields x <= KP and x >= (K-h)P+h-S, hence exactly the stated necessary interval.

The constant bound C >= h is valid even after the endpoint correction: for R=0 that correction vanishes and AL-T >= 0; for R>0, f(R) <= R-1 and E <= h ensure h(2R-f(R))-E first >= 0. Also L >= 2^b+1 and G <= q-h imply K >= q2^b+h. For a nonempty candidate interval x>0. Its unrounded length is [hx-K(h-S)]/[K(K-h)] <= hx/[K(K-h)] <= hN/(q^2 4^b). Thus the strict uniqueness threshold is correct. The search-window complexity is an arithmetic-operation estimate for a fixed core bank and fixed half-width, not a bit-complexity theorem or tail result.

Alpha stripping is exhaustive: odd h forces each last suffix bit, and the inverse update is unique. A genuine remaining seed alpha0=qP+A is positive for every strong recipe, so stopping at a nonpositive seed discards no admissible solution. The upper and lower alpha bounds follow respectively from tau<=alpha and tau>=0. Python integers avoid overflow in the all-length inverse. The rejected alpha0=0 cases cannot meet strongness.

The curvature comparison is valid for the explicitly described unchanged free block and fixed context, including reversal/global complement and transfer of the entire free block to the suffix. The quadratic coefficient is -h times the same distance factor in either location. Therefore h'=dh is necessary. Neither 2 nor 3 maps a member of {5,9} into that bank. This excludes only the stated universal fixed-context rewrites. Restricted languages, multiple duplicated blocks, target-dependent contexts and nonlinear encodings remain outside the result.

## Reproduced and independent checks

- Recompiled balanced.cpp and byte-reproduced all three [23,10000] band receipts: 1123, 215 and 59 misses.
- Byte-reproduced all five width-1 large-target outputs, the width-2 trillion-target repair, and both short-side repair outputs.
- Replayed inverse_alpha.py and reproduced every output, including 74 tested candidates for 47 and 129 for 63.
- Replayed verify.py: 320 literal subset-DP counts, 20,000 interval controls, 5,004 uniqueness controls, and seven large graphs. Its generated graph file remained byte-identical to the frozen input.
- check_independent.py uses returning prepend/append coordinate updates instead of the inverse's f-based formula. It directly enumerates 8,119 candidates under the independent alpha upper bound and checks exact membership for every integer 1 through 1000 against alpha stripping.
- Independently checked every pairwise edge of all seven retained large graphs against the core/cut construction, each Hamilton cycle, and each count via the returning compiler.
- Checked the elementary f bound and equality characterization for P=1 through 99,999, plus every core/suffix constant and denominator inequality for b<=7.

The large graph counts are certified by the audited construction formula and independently evaluated compiler, not by whole-graph subset DP. The bounded-short-side absence proof is outside this audit; only its two positive repairs were checked here.

## Implementation domain qualification

balanced.cpp uses signed long long arithmetic for target/formula evaluation and only uses __int128 for its threshold calculation. It is not an arbitrary-precision implementation, and the b<27 ceiling is not its sole practical domain restriction. The retained inputs (largest about 1.1e12) are safely within the integer range and all receipts pass. Any claim of executable correctness for unrestricted numeric inputs should additionally require that all intermediate values fit signed 64-bit arithmetic, or change the evaluator to wider/arbitrary-precision arithmetic. This does not affect the abstract interval theorem or Python all-length membership procedure.

No source files were edited and nothing was committed or published by this audit.
