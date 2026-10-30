# Independent hostile audit: phase–defect upper-envelope quotient

Auditor: separate session113 worker; producer and NumPy checker were read but neither is imported by this scalar-integer verifier. No Git writes made.

## Mathematical bridge audit

1. **State closes without resetting a source.** From `m'=3m/2^a` and `b'=(3b+1)/2^a`, `z=b+1-m` has exactly the advertised update. Both2 and3 are nonsquares modulo5, so the square flag changes by `a+1 mod2`. Prefix clock A is never normalized away.
2. **Content and guard.** At level e, content requires `m=1,z=0`. Nonsquare m cannot be1 even modulo5. Ones preserve this flag and multiply z by the unit3/2, preserving its truncated height. Discarding finer m is legitimate only because the consumer is an upper bound, not exact return recognition.
3. **Epoch upper envelope.** For a non-one followed by L ones, the L+1 positions have at most `floor(L/(2*5^(e-1)))+1` hits at each e≤h. The order formula follows from `v5((3/2)^2−1)=1`. Summing and removing floors gives `h+lambda_h L`. This remains valid for L=0 and a final partial epoch. Inactive epochs cost0. Using `max(h,k)` only increases the bound.
4. **Initial/empty/all-one paths.** Starting from identity, the initial ones have exact level counts `floor(L/(2*5^(e-1)))`, hence cost≤lambda_E L. Their full clock remains in A. Subsequent quotient paths start at(0,0) after that initial run, with the separate entry inequality paying the first non-one. No unjustified `V(identity)=0` is needed. Empty and all-one words satisfy the bound without invoking a vector.
5. **Cost and sign.** One-edge cost is `D(lambda_h−lambda_E)` in active phase and `−D lambda_E` otherwise; the stated integer formula agrees. Non-one cost is at least `D h−N a` when active. Thus `D S_E−N A≤C`. The direction `V(t)≥V(s)+c` and the entry condition give `C≤V(final)≤0`, not the reversed inequality.
6. **Infinite alphabet.** For every omitted a>64, truncated h,k≤E implies `c≤DE−Na≤DE−65N`. The two separately checked tail inequalities cover all preceding states and the first-entry case. No periodicity assumption or bound on an untruncated valuation is used.
7. **Sharpness.** For t=2*5^(E−1) ones, exact valuation counts sum to `(5^E−1)/4=lambda_E t`. The finite positive source `2^(t+1)−1` has precisely those first t valuation letters. This does not assert an infinite positive all-one orbit.
8. **Resource claim.** Full affine state count is `phi(5^E)*5^E=4*5^(2E−1)`, versus `2*5^E` quotient states. At E6 these are195,312,500 and31,250. The quotient deliberately forgets exact identity incidence, so it does not contradict an exact-observer lower bound.

No defect in these bridges was found. Acceptance is only for precisions whose finite integer certificate passes. No induction in E, all-prime statement, ROOT bound or Collatz convergence follows here.

## Independent implementation

`check.py` uses Python arbitrary-precision scalar integers, struct-decoded signed64-bit certificate values, direct `2**a` numerators and Python modular inverses. Its one-edge reward slope is recomputed as a geometric sum rather than using the producer's closed formula. It checks every edge for a=1,...,64, each entry, vector maximum, and both tails. It also enumerates all262,143 nonempty positive words with total valuation clock≤18 by exact integral(P,Q,B) recursion. At every prefix it checks the actual gcd content bound, envelope bound, and accumulated path bound. These word checks are adversarial diagnostics, not substitutes for the all-word certificate.

Reproduce, for E in6,7,8:

    python check.py E /tmp/collatz_VE.bin --out auditE_independent.json

Producer reproduction is documented in `../collatz/THEORY.md`. The vector hash in each independent receipt must agree with the producer receipt. Vectors are not assumed optimal or canonical; their checked inequalities suffice.

## Completed results

| E | states | independently checked edges | min V | max V |
|---|---:|---:|---:|---:|
|6|31,250|2,000,000|−1,206,060|−156,240|
|7|156,250|10,000,000|−6,392,860|−781,240|
|8|781,250|50,000,000|−32,942,760|−2,289,280|

All three passed, with minimum edge slack0; both tail margins are strictly positive. All three hashes agree with the respective producer receipts. Each precision also passed the262,143 exact-word diagnostics. Thus the fixed-precision theorem is independently accepted at E6,7,8.

### Extra finite-precision corollary (E6 and E7 only)

At E6,E7 the audited maximum is exactly `−2N`. Consequently every word containing at least one letter greater than1 satisfies the stronger

    S_E(w) <= lambda_E (A−2).

The proof is the same path envelope, retaining its negative maximum rather than weakening to0. It is sharp: take `1^t,2` with `t=2*5^(E−1)`. The initial ones attain `lambda_E t`; appending2 makes m nonsquare so adds no content charge. Its total clock is `A=t+2`. Source realizability of finite words is not needed for this upper bound. The ordinary equality example remains the clean source-realizable sharpness control.

**Do not extrapolate this2-clock gap to E8:** its maximum is `−2,289,280`, larger than `−2N=−3,906,240`. The present certificate proves instead the gap `2,289,280/(8*5^8)` for words containing a non-one. This does not disprove a stronger actual-content gap; it limits this particular envelope certificate.

## Second audit: exact obstruction and low-phase repair

The later `collatz/THEORY.md` §6 and `uniform_content/THEORY.md` were independently audited. **Approved with their stated scope.**

- U3=(3,1,3,1^6,3,1,3), followed by231 ones and a second copy: clock271,length255, actual untruncated content154, original budget excess769/5000, E9 envelope excess60112/390625.
- U5=R1^233R1^8R1^233R, R=(5,1^8,5): clock546,length514, actual content317, original budget excess479371/312500, E10 envelope excess1043012/1953125.
- The positive excesses refute the universal budget hypothesis and any potential satisfying the existing weights, all transitions, entry inequalities and nonpositive maximum. They do NOT refute the actual-content bound, nor the already-verified E6–8 certificates.
- Resonant doubling is exactly `z_(X1^LX)=z_X(1+m_X(3/2)^L)` by the affine composition law. The three displayed valuations3,4,3 and all listed intermediate heights are confirmed by exact integer iteration.
- Both words have independently constructed positive integer sources. The auditor uses reverse CRT modulo P=3^length, choosing a positive odd endpoint then solving for the initial source, unlike the producer's forward CRT modulo2^(A+1). Every individual valuation letter is checked afterward.

For the low-phase repair, squares modulo5^e form precisely the cyclic group generated by3/2, of orderT_e=2*5^(e−1). Thus the first offsets d_e exist uniquely in[0,T_e−1], with d_e congruent to d_r moduloT_r and hence d_e≥d_r for e>r. The stated high-layer estimate has the correct upper-bound direction. Low counts are exact, and no tail is charged before d_r. At d_r the pending-bit machine charges h−r once; each subsequent one charges lambda_h−lambda_r, exactly reproducing the formula. A new non-one resets the flag and evaluates the new epoch's first position immediately. Inactive epochs contribute0. The state count at most2*phi(5^r)*5^E is valid for this online epoch consumer; it does not certify a feasible global potential.

The independent checker constructs d_e by successive5-adic lifting rather than the producer's direct scan. It verifies exact epoch counts, both refined bounds, and equality of online and closed-form charges. Refined totals match:

- U3,r1:198661/1250.
- U5,r1:26448046/78125.
- U5,r2:5094453/15625.

All are above actual content and below5A/8. A further24,800-case sweep checks every unit multiplier modulo125, h=0,...,3, r=1,2 and L=0,...,30, including inactive phases. All upper bounds pass.

Reproduction: `python check_new_boundaries.py`. Exact receipts, first offsets and independent positive sources: `new_boundaries_independent.json`. This closes the requested audit; no universal feasibility claim is inferred from the repaired examples or sweep.

## Third audit: repaired one-digit E9 certificate

**Approved:** the repaired consumer now independently proves the sharp E9 bound

    S_9(w) <= (488281/781250) A

for every finite positive valuation word. This is a fixed-precision result, not all-E feasibility or convergence. It does not retract either counterexample to the older, different envelope.

The valid state domain is all four ordinary u modes with z moduloq, and pending-u4 states only with z divisible by5:4q+q/5 states. After a non-one, u1 charges h immediately; u4 with positive h defers the first charge until the next one; nonsquares and height0 charge0. Pending-u4 followed by one charges h and clears the flag. Other active one-steps charge the exact level1 return plus the high-level geometric density. This reproduces the previously audited F_1 formula, including L=0 epochs and words ending before a pending first hit. Initial all-one runs end at u1 oru4,z0; giving either entry phase zero residual deficit is a safe overestimate since those runs were already bounded bylambda_E times their original clock. Each non-one resets epoch logic without changing the original source or clock.

The independently derived array checker uses modular inverses from direct Python powers, obtains h by divisibility tests, and computes high-density weights as a sum of reciprocal periods, rather than the producer's closed formula. It checks the full valid domain and excluded sentinel slots, destinations and closure, every a=1,...,64 edge, both entry phases for every a=2,...,64, nonpositive maximum, and both infinite-letter tails. Closure for omitted non-one letters is also structural: a pending destination is created only if its new z is divisible by5. Gross non-one rewards never exceed E, so both tails apply to every omitted letter.

All gates are explicit `if not ...: raise`, active under python -O. For the supported E≤9 and chunked signed64 arrays, q≤1,953,125 and transition products are below4q²; reward and clock terms are below10^9. Potential magnitudes are explicitly gated below2^50, leaving all additions/subtractions strictly inside signed64 range. Thus the array inequalities use exact integers, not floating point or overflow arithmetic.

| precision | valid states | checked edges | minV | maxV |
|---|---:|---:|---:|---:|
|6|65,625|4,200,000|−1,343,640|−156,240|
|9|8,203,125|525,000,000|−199,232,020|−19,531,240|

At E9 the126 entry inequalities pass. All valid values are finite; the1,562,500 excluded pending-height0 slots have precisely the sentinel. Minimum edge slack is0. Entry-tail margin is294,908,280, transition-tail margin314,439,520. Vector SHA256:

    83e016db40ab86a1ae8d7e6aef605a4710c0d1279f23f1741e11f617a187670f

A separate standard-library scalar Fraction checker independently tested all41,664 valid edges at E1,E2,E3. A corrupted certificate was rejected under python -O. The earlier exact-word and refined-epoch semantic controls continue to support the mathematical bridge, independently of these array tests.

Since maxV9=−2N9, the stronger sharp inequality `S_9<=lambda_9(A−2)` holds for every word containing a non-one. Equality follows from t=2*5^8 ones then a2. This corollary uses the new E9 certificate; it must not be attributed to the old E8 certificate whose maximum did not establish that gap.

Reproduce:

    python -O audit_repaired.py 9 /tmp/repairedV9.bin --output repaired9_independent.json
    python -O check_repaired_scalar.py

The receipt and scalar controls are `repaired9_independent.json` and `repaired_scalar_controls.json`. No E10 experiment was undertaken in this audit.
