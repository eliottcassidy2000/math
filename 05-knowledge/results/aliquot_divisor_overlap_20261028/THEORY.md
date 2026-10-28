# Divisor shape, additive arity and the realization guard

## Inheritance and interpretation

Anchor: formalize the user's20 overlap. Niche: near-perfect excess2.
Wildcard: transport between aliquot and square-sum edges. The live board is
**prime-exponent shape, proper-divisor arity, factor-pair witness, excess,
square total, continuation**.

Closest proved mechanism: [THM-4805 / aliquot-sources-of-2-5-7-21](../../../01-canon/theorems/THM-4805-aliquot-sources-of-2-5-7-21-twin-prime-collisions-the-powers-and-doublings-chain-structure-and-abc-hits-on-the-row-of-one.md),
Lemma1, already classifies the small divisor-count regimes. Canonical hostile:
16->15 has four proper divisors; a count of summands cannot be inferred from
the value alone. Corrected near miss: the [aliquot124 correction lineage](../aliquot_power_124_20261026/CORRECTIONS.md)
keeps distinct-prime and conjectural prime-supply guards visible. Least-used
sidecar: the actual prime-power exponent vector behind the displayed sum.
The [Session82 source-preserving synthesis](../session82_decimal_aliquot_square_sum_20261028/README.md)
already proves the relevant1001 theorem and warns against merging its graphs.

Write s(n)=sigma(n)-n, with s(1)=0. Prime powers mean p^e with e>=1.
Their proper-divisor sum is 1+p+...+p^(e-1). For n=product p_i^e_i,
the number of proper divisors is t(n)=tau(n)-1=product(e_i+1)-1.
The user's prime-cube form has three summands; the prime-fifth-power form
has five. These are arities of a divisor expression, not polynomial degrees
of an invariant map on all integers.

The literal lists for12,18,28,24 all include1. Their sums are16,21,28,36.
The user's clarification is encoded by the pair of requirements

    ab=n,       s(n)=1+a+b,       a,b positive integers.       (C)

This interpretation preserves both multiplication and addition. Arbitrary
regrouping into three summands would have no comparable classification.
When a,b>1 and distinct, these are actual proper-divisor terms. With a=1
the display repeats1 and includes n; it is a numerical compression, with a
different realization guard. Twenty legitimately has both representations.
Sixteen does not pass(C); it remains in the prime-power family.

## P1. Complete classification of the factor-pair compression (PROVED)

For all n>=1, condition(C) is equivalent to

    tau(n)=4  OR  s(n)=n+2.                               (1)

The alternatives are disjoint. Equivalently, the nonnegative discriminant

    Delta(n)=(s(n)-1)^2-4n                                (2)

must be a square. Its roots give a,b=(s(n)-1 +/- sqrt(Delta))/2.

**Proof.** If one factor is1, the other is n and(C) is exactly s(n)=n+2.
If1<a<b, all of1,a,b are distinct proper divisors. Equality with their sum
forces them to be the entire proper-divisor list, so tau(n)=4. Conversely
tau4 means n=p^3 or n=pq with distinct primes, giving pairs(p,p^2) or(p,q).

The repeated-factor case a=b cannot occur. If a is prime, s(a^2)=1+a,
which is smaller than1+2a. If a is composite, choose its least prime p.
Then a>=2p, so a^2/p>=2a; the distinct proper divisors1,a,a^2/p already
have sum at least1+3a, larger than1+2a. This exhausts all factor pairs.

Prime cubes are deficient: s(p^3)<p^3. For a distinct-semiprime pq,
`s(pq)-pq=2-(p-1)(q-1)<=0`, with equality only at(p,q)=(2,3).
Thus tau4 cannot also have s(n)=n+2. Finally a square Delta has the correct
parity because Delta=(s-1)^2 mod4; n1 and primes have negative Delta, while
all other cases have positive root sum and product. This proves(2). QED.

The correct reading of the user's20 equation is therefore

    a*b=20,  a+b=21,  1+a+b=22,  (a,b)=(1,20).

No claim that20 or the repeated1 is a new proper divisor is needed.

## P2. Twenty is the unique five-term overlap (PROVED)

An integer with exactly five proper divisors has tau(n)=6, hence the shape
p^5 or p^2 q with distinct primes. Prime powers are deficient, and cannot
have s(n)=n+2. For n=p^2 q,

    s(n)=p^2+p+1+(p+1)q,
    s(n)=n+2 iff q(p^2-p-1)=p^2+p-1,
    q=1+2p/(p^2-p-1).                                   (3)

Since q>=2, (3) implies p^2-3p-1<=0, so prime p is2 or3. At p2 it gives
q5 and n20. At p3 it gives11/5, not an integer. P1 handles the other possible
compression branch. Therefore **tau(n)=6 and(C) iff n=20**. QED.

This is stronger than a finite observation that20 appears early. Its
uniqueness is relative to FIVE actual proper divisors. There are further
unit-pair overlaps at higher arities.

One can retain the user's three bins by giving prime powers priority:

1. all prime powers;
2. remaining(C) sources: distinct-semiprimes together with s(n)=n+2;
3. all remaining composites.

Number1 is a separate base state at row0. Every member of bin3 has at least
five proper divisors, but bin2 also includes20 and other high-arity overlaps.
For the intrinsic classification, retain both the factorization shape and
the compression witness instead of forcing the two notions of arity to agree.

## P3. The affine divisor rows and classical near-perfect family (PROVED / CITED)

For distinct primes p,q, exponent a>=1 and r=s(p^a), geometric summation gives

    s(p^a q)=1+r(p+q),       r=(p^a-1)/(p-1).             (4)

Indeed sigma(p^a)=pr+1 and sigma(q)=q+1; subtract p^a q. The coprimality
guard is essential: treating p=q=2 as a distinct-prime input would give
the wrong divisor sum for8. The exponent vector chooses the correct chart.

Specialize p2,a2:

    n=4q,       s(n)=3q+7,       s(n)-n=7-q.

At q=3,5,7, this is exactly12->16,20->22,28->28, with excess4,2,0.
Writing x=q+2 gives s(4q)=3x+1, a literal affine identity on this prime-guarded
chart. It is not a conjugacy of iterations: the output usually leaves the
chart and is not re-encoded as the next x. No Collatz conclusion follows.

More generally ALL sources of the form p^a q with distinct primes and
s(n)=n+2 satisfy

    p=2, a>=2, q=2^(a+1)-3 prime.                        (5)

Proof: p^a=(p-1)r+1 turns the excess2 equation into
`q=((p*r)-1)/((p-2)r+1)`. For p>=5 this is less than2. For p3 it is less
than3, so q would be2, requiring r3; but s(3^a) is1,4,13,... and never3.
For p2 the formula gives q=2r-1=2^(a+1)-3; a1 gives1 and is excluded. QED.

The family in(5) is classical: it is the redundant-divisor2 specialization
of [Shevelev, *On perfect and near-perfect numbers*, Theorem3](https://arxiv.org/pdf/1011.6160).
An even n with s(n)=n+2 is near-perfect because omitting the proper divisor2
leaves sum n. The name refers to one omitted divisor; “2-near perfect” is
not used here for this condition. The paper already gives650 as an example
outside(5). We independently check650=2*5^2*13 and s(650)=652.
The prime condition in(5) is retained; no infinitude of that prime family
or completeness for all exponent shapes is claimed.

As further all-size arity checks,104 is the unique seven-proper-divisor
overlap, and464 the unique nine-proper-divisor overlap. For tau8 the possible
shapes are p^7,p^3q,pqr. The first is deficient; (5) gives104 in the second.
An odd pqr is deficient since sigma(n)/n<=4/3*6/5*8/7=64/35<2. If pqr is
even, sigma(n) is0 mod4 but2n+2 is2 mod4, excluding excess2. For tau10,
the shapes are p^9 or p^4q, so(5) gives only464. These are elementary
classifications; the later eleven-term overlaps650 and1952 already have
different exponent shapes.

The finite scan through1000001 finds exactly

    20,104,464,650,1952,130304,522752

with excess2. This is a complete bounded list, not the complete infinite
sequence. Every one is even in this range; no all-order parity claim is made.
For comparison, P1 shows the only PERFECT number satisfying(C) is6. Thus
the common normalized perfect-number triple1+1=2 still loses the distinction
between the three-term source6 and the five-term source28.

## P4. The1001 connections retain their separate operations (INHERITED + VERIFIED)

[Session82's anchor characterization](../session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md)
already proves the user's uniqueness claim, in this exact scope:

    n squarefree with exactly three primes,
    s(n)=ell^3 for prime ell, and sigma(ell^3) square
    iff n=1001=7*11*13, ell=7.

The prime-cube reduction imports Ljunggren's square-repunit theorem. Its
statement is confirmed on printed p257 of the primary research article by
[Bugeaud–Mignotte–Roy](https://msp.org/pjm/2000/193-2/pjm-v193-n2-p.pdf).
We inspected that statement, not a new proof of Ljunggren. The associated
complete prime-power square-edge list is3->1,81->40,343->57, as explained
in [Session82's square-sigma report](../session82_decimal_aliquot_square_sum_20261028/SQUARE_SIGMA_EDGES.md).

After forcing ell7, the remaining uniqueness proof is elementary. For primes
p<q<r with s(pqr)=343, the lower bound1+3p+3p^2 rules out p>=11, and parity
rules out p2. The factor-pair identity is

    (p+q+1)(p+r+1)=343+p(p+1).

For p3,5,7 the right sides355,373,399 have only one admissible factor pair
after imposing q>p and r>q:19*21 at p7, giving q11,r13. We reproduce all
factor pairs and the complete19-source inverse fiber of343. Exactly1001
in that fiber has three distinct primes; it is not the only predecessor.
Completeness has the proved bound n<=(343-1)^2: any composite source has
the proper divisor n/p>=sqrt(n), as well as1.

The independently rerun short orbits are

    1001 ->343 ->57 ->23 ->1 ->0,
    1000001 ->10003 ->1437 ->483 ->285 ->195 ->141 ->51 ->21 ->11 ->1 ->0.

In the first,343+57=400 is a genuine aliquot square-sum edge. The preceding
edge1001->343 has sum1344 and is not square-sum. Separately,
`1001+1001000=1001^2` is the universal square-sum construction
`m+m(m-1)=m^2`. It requires no aliquot relation. We do not re-run or extend
Session82's long orbit from1001000.

The arithmetic signs from the preceding [duality packet](../liouville_mobius_eleven_20261027/README.md)
make the exponent change visible:

| source | prime-exponent shape | proper-divisor count | mu | lambda |
|---|---|---:|---:|---:|
| 1001=7*11*13 | (1,1,1) | 7 | -1 | -1 |
| 343=7^3 | (3) | 3 | 0 | -1 |
| 57=3*19 | (1,1) | 3 | +1 | +1 |

Liouville misses the first support collision because total exponent parity
stays odd. Möbius detects the repeated prime. Neither current sign predicts
the entire next aliquot state. Even the common shape(1,1) at6 and10 has
different successor shapes:6 stays semiprime and10 goes to the cube8.

## Composition lesson and exact scope

The proper-divisor polynomial for p^2q is
`D(x,y)=1+x+x^2+y+xy`. At x2,y5 it equals `2+x^2y` numerically. These are
different formal polynomials: one records every proper divisor once, the
other repeats the unit and includes the source. P2 proves the unique
prime-guarded coincidence of this shape. Evaluation can identify different
incidence data, just as the earlier path and Vitali summaries do.

The relevant carrier is consequently **prime labels and exponents, divisor
incidence, and the chosen compression witness with its guard**. The inverse
aliquot row loses some of those coordinates: row22={20,38}, row15={16,33},
row16={12,26}. A numerical equality in a row is not a replacement for the
factorization or divisor list. This is a precise use of the repo's
source-preserving composition principle, not an asserted isomorphism to a
tournament, Turing tape, or Collatz dynamical system.

No global classification of excess2 numbers, prime supply, aliquot basins,
or LRC(14) is established. The new closed statements concern(C), its exact
five/seven/nine-term overlaps, and the stated finite universe.

## Incoming integration: literal arity, compression, and two parity spaces

The [Session83 low-term atlas](../session83_aliquot_divisor_grammar_20261028/LOW_TERM_ATLAS.md)
independently classifies abundance through nine proper divisors. Its literal
three-divisor class and our C predicate must remain distinct: the former
requires three actual proper divisors, while C explicitly permits the unit
pair(1,n). Thus20 lies outside the former but inside C. No contradiction or
change to either definition follows. Its five-term abundance list12,18,20
and perfect source28 independently explain our P2 boundary.

The [fixed-excess report](../session83_aliquot_divisor_grammar_20261028/FIXED_EXCESS_AND_FIBERS.md)
proves a stronger scope for(5): it classifies ALL excess2 sources with at
most two distinct prime factors. Its elementary reduction is worth retaining.
Two odd prime factors give abundancy strictly below(3/2)(5/4)<2, and prime
powers are deficient. Otherwise n=2^a p^b, with p odd prime and a,b>=1. Put
R=1+p+...+p^(b-1). Then

    sigma(n)-2n+1=(2^(a+1)-p)R.

Excess2 forces the product to be3. If b>=2 then R>=p+1>=4, impossible.
Thus b=1 and p=2^(a+1)-3, recovering(5) with a>=2. This classification is
inherited from that report, which credits Li–Ma–Zhang(2013), Proposition2.3;
our earlier Shevelev citation already supplies the family itself. The source
650 has three distinct prime factors and remains the necessary scope guard.

The [144 bridge](../session83_aliquot_divisor_grammar_20261028/CARRIER_AND_144_BRIDGE.md)
adds another exact connection with the earlier{2,3,11} work. Define Q(x) to
be the vector(v_p(x) mod2)_p, of finite support. It records a square class,
not merely the bit saying whether x is square. For coprime a,b,

    Q(sigma(ab))=Q(sigma(a))+Q(sigma(b)) over F2.

The source signs mu(n),lambda(n) concern the exponents of n; this second
parity space concerns the exponents of sigma(n). They are different data.
A cheap hostile example for merging them is a=2 versus a=5. Both source
sign pairs are(-1,-1), and both sigma values3,6 are nonsquare. The common
legal multiplier11 has sigma12, with square class3. Consequently

    sigma(2*11)=3*12=36 is square,
    sigma(5*11)=6*12=72 is not square.

The guards gcd(2,11)=gcd(5,11)=1 hold. Thus even the joint summary of both
source signs and the square/nonsquare bit of sigma fails this continuation.
Keeping Q(sigma) suffices for this predicate under coprime multiplication.
This is not a claim of a smallest carrier for every arithmetic operation.

At the next step sigma(3)=4 has zero square class, so multiplying22 by3
gives66->78 with square total144. The alternative triple70=2*5*7 also has
sigma144 and edge70->74. The precise connection is cancellation in a
prime-indexed parity vector, with coprimality as realization guard. Source
factorization/divisor incidence, total square class, and the chosen observed
predicate remain distinct coordinates even when numerical totals agree.
