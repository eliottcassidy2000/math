# Why 1001 is the unique three-distinct-prime predecessor

**Theorem.** Among squarefree integers having exactly three prime factors, **1001=7·11·13 is the unique integer whose aliquot successor is a prime cube with square divisor sum**.

Here s(n)=sigma(n)−n. The relevant steps are

    1001 -> 343=7³ -> 57,
    sigma(1001)=1344,   sigma(343)=400=20².

The second step is square-sum; the first is not. The theorem describes a predecessor of the exceptional square-edge source, not a square-edge source itself.

## Proof

The classical square-repunit theorem of Ljunggren, in the precise scope checked in the [companion prime-power classification](SQUARE_SIGMA_EDGES.md), gives

    sigma(ell³) square, ell prime  <=>  ell=7.

Thus it suffices to classify distinct primes p<q<r satisfying s(pqr)=343. For this squarefree product,

    s(pqr)=1+p+q+r+pq+pr+qr >= 1+3p+3p².

Since 1+3·11+3·11²=397>343, the least prime p belongs to{2,3,5,7}. If p=2, both q and r are odd, and sigma(2qr)=3(q+1)(r+1) and2qr are even, making s(2qr) even. This excludes the odd target343.

For each remaining p, regroup the exact equation as

    (p+q+1)(p+r+1)=343+p(p+1).

Both factors are positive; their order is strict because q<r. All factor pairs are:

| p | right-hand side | increasing positive factor pairs |
|---:|---:|---|
| 3 | 355 | (1,355),(5,71) |
| 5 | 373 | (1,373) |
| 7 | 399 | (1,399),(3,133),(7,57),(19,21) |

For p=3 the smaller factor must be at least3+5+1=9, excluding both pairs. For p=5 it must be at least5+7+1=13, excluding its only pair. For p=7 it must be at least7+11+1=19, leaving exactly(19,21). These give q=11,r=13. Direct substitution gives s(1001)=343. This proves uniqueness.

**External dependency.** Only the first reduction imports Ljunggren's theorem; the classification of squarefree three-prime sources of343 is elementary and complete. The unrestricted square-repunit theorem is explicitly stated in [Bugeaud–Mignotte–Roy, Pacific J. Math.193(2),2000,p257](https://msp.org/pjm/2000/193-2/pjm-v193-n2-p.pdf), citing Ljunggren, Norsk Mat. Tidsskr.25(1943),17–20. No new proof of that classical theorem is claimed.

## Independent complete inverse check and scope controls

Any source n with s(n)=343 is composite. If a is its least prime factor, the proper divisor n/a is at least sqrt(n) and differs from1. Therefore

    343 >= 1+n/a >= 1+sqrt(n), so n <=342²=116964.

A full proper-divisor sieve through this proved bound gives exactly

    578, 1001, 1685, 3641, 9077, 9641, 16697,
    17141, 19241, 19637, 20777, 24341, 24617,
    25397, 25877, 27641, 28757, 28841, 29177.

Of these19 sources, precisely1001 has three distinct prime factors. The source578=2·17² is not squarefree; the other17 sources are products of two distinct primes. Thus the stronger claim “1001 is the only predecessor of343” is false. Every listed source is independently checked by the multiplicative divisor-sum formula, and its complete trial-division factorization is recorded.

As a secondary consistency check, the formal three-distinct-prime divisor expression gives

    s(p(p+4)(p+6))−p³ = −(p−7)(p²+4p+5)

whenever p,p+4,p+6 are distinct primes. Its only positive equality case is p=7, since p²+4p+5=(p+2)²+1>0. The prime guards are indispensable: the expression is not a formula for s of an arbitrary composite polynomial input. This restricted constellation identity is weaker than the generic three-prime classification above.

## Reproduction

    python3 check_anchor.py > certificates.json
    python3 -O check_anchor.py > certificates.optimized.json
    cmp certificates.json certificates.optimized.json

Standard library only; explicit runtime checks survive `-O`. The script exhausts the proved inverse bound, re-encodes all19 sources via independently factored sigma, enumerates every factor pair in the proof, and checks the constellation identity by polynomial coefficient multiplication. Finite checking supplements the displayed proof; it does not replace the cited Ljunggren dependency.

Independent audit reproduced the full19-source fiber with a separate divisor sieve and checked every factor-pair case and the polynomial identity. The uniqueness statement is restricted to the stated squarefree three-prime source class.
