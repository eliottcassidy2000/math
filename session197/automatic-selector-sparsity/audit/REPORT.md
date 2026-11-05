# Independent review: automatic-selector sparsity

**PASS. No correction.** This review binds the frozen automatic-selector-sparsity PROOF.md SHA2568b89586e6ea48e19fd8fe0d2a955aee52f912e98767e02013037ce8ce2472b7a, verify.py d92e452ea39eaf7a6ea70c2bdc64f2e9f2931dceb27aad6d0e99e95add6d000c, and normal/optimized receipts b5177e9b270df95281188abbbd7d3394ab3056857a0d6287eaf5dfd405a790e3. The source was copied under source/ and left unchanged.

## All-order argument

The theorem applies to a synchronous regular subrelation of a fixed physical final-gadget function, not merely to its complete graph. Since that relation is functional, every input in its projected domain brings its exact function value. Regular addition handles the L=n+1 shift. The common-zero left quotient before canonical intersection is correctly included, so projection cannot silently lose pairs when an eliminated auxiliary track was longest.

A complete live strongly connected component supplies the full two-parameter rectangle p1^(A+Ti)0^(B+Tj)s by deterministic eventual cycles under1 and then0, with one fixed accepting suffix. A leading-one reaching prefix exists for a canonical positive-integer domain. The suffix may be empty. I independently expanded the objective on this rectangle: B0=h(P+1)>0, C0=h(1+q_s)+R2^t>0 and A0=(P+1)C0; the stated D0 is correct.

The auxiliary Y=A0 2^(a+b) is genuinely synchronous regular to recognize jointly with the restricted input: its binary pattern is fixed-prefix/zero-tail and its length differs from the input by a fixed constant. This does not require a numerical register computing a+b. Fixed-integer addition and existential projection then recover the unique positive Z=B0 2^a+C0 2^b. Applying common-padding normalization at this second projection is again essential and correctly stated.

For bounded odd part, positivity of both terms excludes cancellation at unbounded unequal exponents. The displayed valuation bound contains every allowable exponent difference, including equal valuations. Choosing the bound from the diagonal makes the finite difference set nonempty. Fixed prefix/suffix quotients preserve regularity, while the explicit widely separated one-runs and corresponding zero-run suffixes give infinitely many distinguishable prefixes. Thus the required regular diagonal-band language is impossible.

If no complete live SCC exists, every live state has a path out of the live set. A shortest such path visits at most the number of live states before exit; after exiting a reachable coaccessible set it cannot regain acceptance. At least one length-s0 block is therefore lost uniformly at every live state. This proves the strict entropy bound, and summing lengths yields a strict subunit power of the numeric input bound. The lambda=1 case is correctly absorbed into any fixed positive power below1. Shifting n by1 does not change that conclusion.

The target-count bound uses F(L)>=RL, not an assumed injectivity of F. Finite explicitly recoverable tags preserve regularity under tag restriction, and a finite sum of the resulting strict power bounds remains a strict power bound. The central-binomial valuation argument proves the required elementary lower bound pi(X)>=cX/log X; no prime number theorem is needed. Consequently the covered fraction of prime or seven-times-prime targets tends to zero.

The scope exclusions are necessary and correctly maintained: this does not rule out a cofinite range of a fixed prefix, a finite bank with a nonautomatic witness selector, an untagged union lacking regular prefix recovery, a growing bank, or unrestricted low-cycle coverage. In particular it does not contradict the reviewed arithmetic square-root inverses.

## Independent controls

check.py imports no author code. Its rectangle evaluator uses an unfolded sum over zero-digit positions rather than the source's halving recurrence. It verifies274176 full rectangle values, including4352 empty-suffix cases, and274176 forced-padding projection controls. It separately checks2895 odd-part bands, including coefficients with large unequal2-adic valuations, and312660 explicit distinguishing-prefix comparisons.

The automaton controls enumerate all5898 binary transition/acceptance fixtures with one through three states, with a separate explicit leading-one canonical wrapper. They cover4318 fixtures with complete live components,371 leaky fixtures and1209 empty fixtures. There are69088 accepted rectangle points and527 shortest-escape/uniform-block checks. These tests corroborate the finite-state lemmas; the all-order conclusion comes from the proof above, not enumeration.

The independent checker also verifies the elementary binomial/prime-count bound through n300. Both independent normal and optimized runs match exactly. Copied source normal/-O runs reproduce the frozen receipt exactly.

Reproduce:

    python check.py > independent-normal.json
    python -O check.py > independent-optimized.json
    python source/verify.py > source-normal.json
    python -O source/verify.py > source-optimized.json

None of these scripts writes other files. They require only the Python standard library.

Hashes:

- check.py:176eea274bc53af892070c7755e3a387d0d01f938a7b1185d7e7e255fbed6cf9.
- independent-normal.json and independent-optimized.json:39bb4bf5e70fa0129ee512c167192a9eae06f1f87b72e9c28cb666fddc0aa34f.
- source-normal.json and source-optimized.json:b5177e9b270df95281188abbbd7d3394ab3056857a0d6287eaf5dfd405a790e3.

No mathematical correction or unresolved blocker remains.
