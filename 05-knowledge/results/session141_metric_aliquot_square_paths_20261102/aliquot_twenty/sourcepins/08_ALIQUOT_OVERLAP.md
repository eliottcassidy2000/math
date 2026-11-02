# Inner/outer divisor-pair compression: why 20 belongs to two grammars

2026-10-28. **PROVED elementary classification and one-tower criterion;
CITED classical context; FINITE-EXACT tests; OPEN global consequences.**
Internal adversarial review, not an independent-agent audit. No priority claim.

The user's clarification supplies the missing multiplication guard. We interpret
a **cubic representation** as positive integers a,b satisfying

    ab=n,                 s(n)=1+a+b,

where s(n)=sigma(n)-n. The corrected arithmetic at 20 is
`20*1=20` and `20+1=21=2+4+5+10`. Thus the product remains the original
source while the pair sum reproduces the sum of its nonunit proper divisors.
The outer pair (1,n) is allowed explicitly. If both entries must be nonunit
proper divisors, 20 is excluded. This distinction is the mechanism.

The listed 16 does not have this representation: s(16)-1=14, whereas its
factor-pair sums are 17,10,8. The nearby 6 does. We do not silently replace
16 by 6; the theorem states exactly which interpretation has been settled.

## Inheritance and board

Closest proved mechanism: [THM-4805, three-proper-divisor regime](../../../01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md)
already classifies divisor counts, and [Session78](../session78_fourth_power_aliquot_20261027/ALIQUOT.md)
already proves exact source reconstruction from (s(n),s(2n)). Neither is new here.
June [HYP-2220 / Vieta carriers](../../hypotheses/HYP-2220-vieta-perfect-aliquot-carriers.md)
and [HYP-2221 / aliquot fixed points](../../hypotheses/HYP-2221-perfect-aliquot-fixed-point-carrier.md)
are OPEN syntheses; only their displayed elementary identities are inherited.
The current [Session82 anchor](../session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md)
and [shared-port obstruction](../session82_decimal_aliquot_square_sum_20261028/FIBER_PORTS.md)
are stronger than a fresh small census and are retained intact.

Board: **factor-pair involution / divisor count / abundance / Vieta roots /
odd-part carry / named ports**. Canonical hostile: 16 has a square-root fixed
divisor; counting it twice invents a divisor occurrence. Corrected near miss:
an arbitrary expression t=(t-2)+1+1 says nothing about its source; the product
constraint makes it meaningful. Least-used sidecar: whether the witnessing
pair is the outer pair or an actual inner pair.

## 1. Complete classification of the cubic representation

**Theorem.** For every positive integer n, the following are equivalent:

1. There are positive integers a,b with ab=n and s(n)=1+a+b.
2. Either tau(n)=4, or sigma(n)=2n+2.

The alternatives are disjoint. In the first, the unique unordered witness is
the nonunit divisor pair; n is p^3 or pq for distinct primes. In the second,
the unique unordered witness is (1,n). No witness has a=b.

**Proof.** The cases n=1 and prime n fail directly. If a=1 or b=1, the
equation is precisely s(n)=n+2. Otherwise a,b are proper nonunit divisors.
For a!=b, their sum equals the sum of *all* nonunit proper divisors only
when there are no others: tau(n)=4. Positivity is decisive.

If a=b=sqrt(n)=r, a prime-square source has only the one nonunit proper
divisor r, so its sum is r, not 2r. If n is a square with at least five
divisors, it has another distinct nonunit pair d,n/d. Their sum is at least
2r, and the fixed divisor r contributes once more. The total is at least
3r>2r. Thus the diagonal cannot witness the rule. Conversely tau=4 gives
exactly the required pair, and sigma=2n+2 gives (1,n).

For tau=4 the inner pair a,b>1 has a+b<n+1, since
`n+1-a-b=(a-1)(b-1)>0`. Thus it cannot also be the outer case. This also
proves witness uniqueness. ∎

**Exact Vieta test.** Put I=s(n)-1. The representation exists iff I>0 and
`I^2-4n` is a nonnegative integer square. Its square root has the same parity
as I, so `(I±sqrt(I^2-4n))/2` are positive integral factor ports. For the
exceptional outer case the discriminant is exactly `(n-1)^2`. This is an
exact reuse of the June Vieta idea, with a different polynomial and a proved
image criterion; it is not a transfer of the June LRC conjectures.

## 2. Cubic/quintic overlap is uniquely 20

The literal divisor rule has tau(n)-1 summands, including 1. Three summands
mean tau=4; five mean tau=6. The latter sources are exactly p^5 and p^2 q,
where p,q are distinct primes. This is already THM-4805's divisor-count fact.

**Corollary.** The only source with five proper divisors and a cubic
representation is 20.

**Incoming attribution:** during this session, Session83 independently
published the stronger complete five-term abundance classification
`{12,18,20}` with unique perfect source28 in
[DIVISOR_TERM_CLASSES.md](../session83_aliquot_divisor_grammar_20261028/DIVISOR_TERM_CLASSES.md).
Our displayed short proof of its20 corollary is retained as a cross-check,
not a new result over that incoming atlas. The all-n factor-pair equivalence
in §1 permits the outer pair and addresses a different representation question.

**Proof.** By the theorem it must have sigma(n)=2n+2. A prime power is
deficient: sigma(p^a)/p^a < p/(p-1)<=2. Thus n=p^2 q. Substitution gives

    q(p^2-p-1)=p^2+p-1,     q=1+2p/(p^2-p-1).

At p=2 this gives q=5; at p=3 it gives q=11/5; and at every p>=5 it gives
1<q<2. Only n=20 survives. ∎

| n | s(n) | actual nonunit proper divisors | cubic witness |
|---:|---:|---|---|
|6|6|2,3|(2,3)|
|8|7|2,4|(2,4)|
|12|16|2,3,4,6|none|
|16|15|2,4,8|none|
|18|21|2,3,6,9|none|
|20|22|2,4,5,10|(1,20)|
|28|28|2,4,7,14|none|
|104|106|2,4,8,13,26,52|(1,104)|
|650|652|2,5,10,13,25,26,50,65,130,325|(1,650)|

Prime-power status, actual summand count, and compressed-pair status are
different coordinates. For example 8 is both a prime power and literally
cubic. Prime powers may be assigned first in a disjoint reporting table,
but doing so should not erase their other flags. The intended three classes
are therefore best stored as an overlapping classification with labelled
representations, not as three mutually exclusive intrinsic species.

### The next completed layer: eleven proper divisors

**Theorem.** The sources with eleven proper divisors and an outer cubic
representation are exactly650 and1952. Consequently the complete list of
outer collisions with at most eleven proper divisors is

    20,104,464,650,1952.

This extends the incoming [at-most-nine-term atlas](../session83_aliquot_divisor_grammar_20261028/LOW_TERM_ATLAS.md)
specifically for excess2; it is not a full eleven-term abundance atlas.

**Proof.** Tau=12 gives precisely the shapes p^11, p^5 q, p^3 q^2, p^2 qr,
with distinct prime letters. Prime powers are deficient. For the two-prime
shapes, use the complete excess-two classification in §3 (also independently
proved by [Session83's fixed-excess factorization](../session83_aliquot_divisor_grammar_20261028/FIXED_EXCESS_AND_FIBERS.md)):
the source must be2^a p, with p=2^(a+1)-3 prime. Tau=2(a+1)=12 forces a=5,
giving1952=32*61.

It remains to solve the p^2 qr shape. If all three primes are odd, the
abundancy ratio is at most

    (1+1/3+1/9)(1+1/5)(1+1/7)=624/315<2.

To justify the placement of the square, upgrading exponent1 to2 multiplies
its factor by `1+1/[p(p+1)]`, decreasing with p; hence the maximum uses the
smallest prime for the square. Monotonicity then selects3,5,7.

If the squared prime is2, write n=4qr with odd q,r. Then sigma(n) is divisible
by4, whereas2n+2 is2 mod4, impossible. In the remaining even case write
n=2p^2 r with odd distinct primes p,r. Direct substitution in sigma(n)=2n+2
gives

    r=(3p^2+3p+1)/(p^2-3p-3)=3+(12p+10)/(p^2-3p-3).

At p=3 this is negative. At p=5 it is13, giving650. At p=7 it is169/25.
For every p>=11 it is strictly between3 and5, so cannot be an odd prime.
This proves the eleven-term assertion. Ten terms means a single prime to
the tenth power, so deficiency excludes it. The incoming at-most-nine atlas
supplies exactly20,104,464 in the preceding layers. ∎

Here the old p^2 qr theme acquires an explicit new predicate: excess2, with
650 as its unique solution in that exponent shape. No identification with
the repository's unrelated F=U+S quantities is inferred from the notation.

## 3. Classical abundance-two context and the exact doubling transversal

The exceptional equation is `sigma(n)-2n=2`. For even n>2 it means that
deleting the proper divisor 2 makes the remaining proper divisors sum to n:
a near-perfect number with redundant divisor 2, in the terminology of
[Pollack–Shevelev, *On perfect and near-perfect numbers*, §1–2](https://www.pollack-math.net/pnp-4.pdf).
This established branch is not a newly discovered number class. For a
hypothetical odd solution, 2 is not a divisor, so that terminology does not
apply with the same definition. Li–Ma–Zhang use a narrower convention in
[their 2013 paper](https://publi.math.unideb.hu/paper/1759/download/); their
Theorem 1.1 requires any odd solution to have at least six distinct prime
factors. We import no all-odd nonexistence theorem.

Write n=2^a m, m odd. Session78's exact carry law gives, by subtraction,

    sigma(2^a m)-2^(a+1)m = 2^(a+1)s(m)-sigma(m).

**One-tower criterion.** For each odd m>1 there is at most one a>=0 for
which the outer cubic representation exists. It exists exactly when

    (sigma(m)+2)/s(m) = 2^(a+1),

or equivalently `(2^(a+1)-1)s(m)=m+2`. The right-hand side of the defect
formula strictly increases with a, because s(m)>0. For m=1 its value is
always -1, so there is no hit. This is a transversal of doubling towers,
not a statement about forward aliquot or Collatz orbits.

For m=p odd prime, s(m)=1 gives p=2^(a+1)-3. Whenever that value is prime,
`n=2^a(2^(a+1)-3)` is an outer collision. This is the k=1 case of
Pollack–Shevelev's Proposition 1, not a new prime-producing theorem.
The stronger converse restricted to exactly two distinct prime factors is
Li–Ma–Zhang Proposition 2.3. The source 650=2*5^2*13 is outside that family:
its m=325 has s(m)=109, sigma(m)=434, and 436/109=4 selects a=1 exactly.

The complete interval n<=10^6 has outer collisions

    20,104,464,650,1952,130304,522752.

This reproduces the corresponding prefix of the classical table. It neither
proves infinitely many such numbers nor excludes larger odd solutions.

## 4. Square-sum edges and the 1001 anchor retain different predicates

For any genuine aliquot edge, its endpoint sum is sigma(n). In particular
1001->343 is not square-sum, while 343->57 is, with total400. Session82
already proves that 1001 is the unique *squarefree three-prime* source feeding
a prime cube with square divisor sum; the deep input is Ljunggren's theorem,
explicitly checked in [Bugeaud–Mignotte–Roy, p257](https://msp.org/pjm/2000/193-2/pjm-v193-n2-p.pdf).
We retain that dependency rather than rebranding a finite search as uniqueness.
Our complete sieve independently recovers all19 sources of343 and isolates
1001 as the only squarefree three-prime one.

The identity `1001+1001000=1001^2` is a different square-sum graph edge,
an instance of `n+n(n-1)=n^2`; its second endpoint is not s(1001).
The exact orbit of1000001 is the already-certified

    1000001,10003,1437,483,285,195,141,51,21,11,1,0.

There is also a clean exclusion: an **even** outer collision cannot be an
aliquot square-sum edge, since sigma(n)=2(n+1)=2 mod4. For an odd outer
collision, square sigma would require n=2u^2-1. Incoming Session83 proves
the much stronger necessary condition of at least356 distinct prime factors
and1156 digits in [SQUARE_EXCESS_OBSTRUCTION.md](../session83_aliquot_divisor_grammar_20261028/SQUARE_EXCESS_OBSTRUCTION.md).
Its proof uses the restriction p=±1 mod8 and an exact prime-product crossing;
we inherit that result without re-running its separate rational certificate.
Existence is not settled here.
Thus the exceptional compression at20 and the square edge at343 describe
different structures, despite both having small certificates.

## 5. What transports to the ongoing ROOT and tournament work

The divisor lattice is the exponent box `product_i {0,...,e_i}` for
`n=product_i p_i^e_i`. Complementing a divisor sends b_i to e_i-b_i.
Its two-cycles are factor pairs; a fixed point exists precisely when n is
a square. The outer orbit is labelled separately. Evaluation sends a
vertex b to the integer product p_i^b_i and then sums the selected vertices.
At20, two inner pair sums `(2+10)+(4+5)` equal the one outer pair sum1+20.
This equality does not identify their vertices or their external incidence.

**Connection contract.** Source: named divisor occurrences and the fixed
integer n. Target: a shorter expression with the same source product and
scalar aliquot sum. Preserved: n and s(n). Lost: proper-divisor membership,
multiplicity, inner/outer role, and graph attachment ports. Repair: keep the
factorization, selected pair, and the certified expansion equality. Cheapest
hostile: the claimed inner ports1,20 are not both proper divisors of20.
For square-sum attachments, Session82's complete fiber F13={27,35} is stronger:
9 attaches to27 but not35, so retaining the value13 alone does not transport
the attachment. The new representation must not erase the named source.

The same principle informs [the chart ensemble](README.md#r1-three-charts-remove-artificial-intermediate-debt):
temporary descriptions may change when their endpoint observations agree;
source weights and external ports must not be reassigned. Divisor compression
is a scalar identity, whereas the chart identity is quadratic. No map between
their dynamics or proof of universal ROOT cancellation is claimed.

**Next precise questions.** Can a prescribed port signature be retained while
compressing an abundance-two divisor packet? Can the odd-part criterion
bound the number of outer collisions using the carried pair (s(m),sigma(m))?
Can the source-preserving ROOT block controller admit a comparable monotone
carried defect? The first has explicit finite tests; the second leads to
existing near-perfect-number questions; the third remains a research target.

## Verification

`check_aliquot_overlap.py` compares a complete divisor-incidence sieve through
10^6 with the Vieta test at every integer, and directly enumerates divisor
pairs plus independently factored sigma for all n<=20000 and every larger
exception. Square fixed points and prime powers are included. It checks
doubling towers, restores the complete343 fiber using its proved342^2 bound,
replays both short orbits, and tests false source/port identifications.

    python3 check_aliquot_overlap.py > aliquot.json
    python3 -O check_aliquot_overlap.py > aliquot.optimized.json
    cmp aliquot.json aliquot.optimized.json

All validation uses explicit exceptions. Finite checks supplement the displayed
proofs and do not replace their quantifiers or the cited classical theorem.
