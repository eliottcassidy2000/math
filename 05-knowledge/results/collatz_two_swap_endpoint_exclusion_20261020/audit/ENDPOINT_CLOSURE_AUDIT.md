# Independent audit of the endpoint large-range closure

20 October 2026. Status: **VERIFIED proof application and FINITE-EXACT independent certificate.** This is separate from the exact-content packet. It concerns only the two specified legal endpoint chains with p>=3, p/2<r<p, gcd(p,r)=1, q=p+r.

## Verdict

The absolute-gap direction, logarithm intervals, continued-fraction completeness, Legendre threshold, and Matveev cutoff all pass. The producer's large-range result excludes integrality for p>=100, with either sign of D. The separate all-clock residual below 100 also passes independent direct-word verification. It is a separate dependency, not part of the large-range producer's advertised 63 checks.

The final text of `ENDPOINT_TWO_EDIT_CLASSIFICATION.md` and `PRIMARY_SOURCES.md` was reviewed after both certificates were complete. Its entire eighteen-row theorem table matches the independent receipt exactly. The stated conclusion is justified: no integral output for either endpoint family at p>=3 in the displayed coprime band, at any cut or either clock sign; the p=2 boundary retains precisely the two rotations of the declared negative cycle. No mathematical or scope correction is required, and the candidate theorem may be promoted to **PROVED** with the cited external dependency retained.

## 1. Correct quantities and inequality directions

Assume an endpoint word is integral. The already audited content theorem then gives, for p>=25,

    |U|=|2^q 3^(-p)-1|=|D|/3^p=g/3^p<exp(-p/9).

The nonzero multiplicative expression U and the additive logarithm L=q log2-p log3 must not be conflated. Since U>-1 and |U|<1/2,

    |L|=|log(1+U)| <= |U|/(1-|U|) < 2 exp(-p/9).

This works for both positive and negative D. A positive-only bound log(1+U)<U is insufficient for the negative-clock branch.

Writing alpha=log3/log2 gives

    |q/p-alpha| < 2 exp(-p/9)/(p log2).

For p>=100, 4p exp(-p/9)<log2. It is enough to verify this at 100 because p exp(-p/9) decreases for p>9. An independent exact threshold check uses

    400 exp(-100/9) < 400/2^11 < log2.

Thus the approximation error is strictly less than 1/(2p^2). As gcd(q,p)=1, Legendre's criterion forces q/p to be a principal continued-fraction convergent of alpha. Intermediate convergents cannot supply an untested extra case.

## 2. External theorem and its exact application

The real case of Bugeaud–Mignotte–Siksek, [Theorem 9.4, printed page 16](https://arxiv.org/pdf/math/0403046), was inspected independently on 20 October 2026. It is a bound for a nonzero product of powers minus one. Apply it over Q to bases 2,3 and exponents q,-p, taking degree 1, exponent bound q<2p, and height parameters 7/10 and 11/10. Since 2^(9/2)<23, the prefactor is below

    (7/5)*30^5*23*(7/10)*(11/10)=602494200<10^9.

Unique factorization makes U nonzero. Consequently

    log|U| > -10^9(1+log(2p)).

This source-dependent lower bound is applied to U, not to L. Comparison with log|U|<-p/9 requires

    p/9 < 10^9(1+log(2p)).

At p=10^12 the reverse strict inequality holds, since log(2*10^12)<30 and 10^12/9>31*10^9. The difference

    p/9-10^9(1+log(2p))

is increasing from that point because its derivative is 1/9-10^9/p>0. Hence every possible integral endpoint must have p<10^12. The boundary p=10^12 is excluded, not accidentally left in either finite range.

## 3. Exact continued-fraction certification

The producer bounds log2 and log3 using rational truncations of the positive series

    log((1+t)/(1-t))=2 sum_(j>=0) t^(2j+1)/(2j+1).

After N terms the omitted tail is positive and bounded above by

    2 t^(2N+1)/((2N+1)(1-t^2)).

The ratio interval [L3/U2,U3/L2] has the correct direction. Equality of the floor at both interval endpoints certifies each partial quotient, and reciprocal inversion reverses the interval endpoints as required.

The independent script uses a different decomposition,

    log3=(3 log2+log(9/8))/2,

with 112 terms for t=1/3 and t=1/17, rather than the producer's 220-term direct t=1/2 series. It separately computes convergents from both rational endpoints and checks agreement through the first denominator at least 10^12.

The resulting prefix is

    1,1,1,2,2,3,1,5,2,23,2,2,1,1,55,1,4,3,1,1,15,1,9,2,5,7.

There are exactly 18 convergents with 100<=p<10^12. The first is q/p=485/306. The last in range is 1193652440098/753110839881; the next, terminal convergent is 8573543875303/5409303924479. Therefore every denominator in the bounded candidate range is covered, including the empty interval 100 through 305.

For each candidate the independent interval for L is

    [q L2-p U3, q U2-p L3].

It avoids zero, and its absolute lower bound exceeds

    2 * 2^(-min(floor(p/9),256)).

This rational quantity is an upper bound, not a lower bound, on 2 exp(-p/9): e>2, floor(p/9)<=p/9, and the cap at 256 can only enlarge it. Thus each exact comparison excludes the candidate for both clock signs.

## 4. Independent replay and residual cross-check

`audit_endpoint_large_range.py` gives 52 explicit checks and byte-identical normal/-O receipts. It verifies the producer script hash and receipt equality, independently reconstructs the logarithm bounds and all 18 convergents, and verifies every absolute-gap contradiction.

`audit_small_endpoint_residual.py` independently reconstructs all 1,502 small shapes and 3,004 paths, including the two p=2 boundary cases. It uses occupancy-set moves and a reverse carry scan, imports no producer implementation, and matches every producer case through the full trace digest

    c76be7354efe56ce1706f26ccbefdc58f152700680f75e32f1eb9a2364468be4.

All 3,002 endpoint instances with 3<=p<100 have g<|D|. The smallest reduced denominator |D|/g is 5, at p=3,r=2 for both families. The p=2,q=3 exceptions have endpoint fixed points -7 and -5; direct parity iteration verifies both belong to the same three-step negative cycle.

Both audit scripts have byte-identical normal and optimized receipts. Reproduce from this audit directory:

    python3 audit_endpoint_large_range.py > large-range.receipt.replay.json
    python3 -O audit_endpoint_large_range.py > large-range.receipt.optimized.replay.json
    python3 audit_small_endpoint_residual.py > small.receipt.replay.json
    python3 -O audit_small_endpoint_residual.py > small.receipt.optimized.replay.json
    cmp large-range.receipt.json large-range.receipt.replay.json
    cmp large-range.receipt.json large-range.receipt.optimized.replay.json
    cmp small.receipt.json small.receipt.replay.json
    cmp small.receipt.json small.receipt.optimized.replay.json

The large-range audit defaults to the sibling `../large-range` producer directory; `--producer-dir` can override it without changing mathematical inputs.

No claim for other two-edit displacement patterns follows from this closure. The proof depends on the special linear common-root collapse and the endpoint content theorem.
