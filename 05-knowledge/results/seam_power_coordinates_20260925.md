# Scale, fourth powers, and exponent exchange: which information can be discarded?

**Date:** 2026-09-25. **Status:** PROVED elementary identities and scoped
reductions; INHERITED PROVED exponent-order and divisor classifications;
FINITE-EXACT controls; CITED Mills and odd-cycle results in the linked
lanes. No equivalence of Collatz, square-sum Hamiltonicity, and planarity
is asserted. No novelty claim or new canon identifier.

The useful common question is operational: after compressing a number,
power expression, or graph, can the next allowed operation and its target
predicate still be recovered? This continuation supplies positive decoders
and exact boundaries for several suggested compressions.

## Inheritance, portfolio, and concept board

Closest proved mechanism:
[THM-4107, gcd-normalized exponent tournament holonomy and LRC blindness](../../01-canon/theorems/THM-4107-gcd-normalized-exponent-tournament-holonomy-and-lrc-blindness.md).
It already proves the scalar nature of raw exponent comparison, the
rational equality curve, and stabilization of higher primitive power
comparisons. The divisor classification is inherited from
[arithmetic braids, DB1--DB3](arithmetic_braids_20260917_divisors.md).
The [intrinsic halving construction](decoder_halving_20260925.md) and its
[target-dependent repair](decoder_pair_repair_20260925.md) supply the
current decoder boundary. Canonical hostile: the raw tie 2^4=4^2 is
destroyed by common dilation. Corrected near miss: both earlier unconditional
degree-two forcing and a label-free divisor-profile quotient lose needed
data. Least-used sidecars: signed gap, prime-support overlap, and the
support/center labels of edge objects.

Anchor: explain the user's scale/power/exponent relations. Niche: recover
a sound finite controller for divisor balance. Wildcard: make triangular
ten into an actual graph decoder. The six live concepts are:

| Concept | Predicate and operation | Decisive test |
|---|---|---|
| Sum/product gap coordinates | recover A,B and transport fourth powers | A=2,B=4 |
| Raw and primitive exponent orders | compare after scale or power change | tie versus a gcd-induced triangle |
| Divisor balance | multiplication of named primes | 6 and 15 followed by multiplication by four |
| Odd projection | remove even terms or even cycle lengths | fifth-order terms in a regular five-tournament |
| Prime power shells | choose and decode a prime child | cube gaps 3,30,6 and quartic exclusion |
| Graph growth / triangular ten | degree-two ear, edge-owner pairing | Q15 succeeds while G15 remains planar |

The META-PATTERNS card used is "Type every analogy and every implication";
the repair-quotient card also applies to the different costs of arbitrary
and prescribed tournament quotients. No new meta-pattern is promoted here.

## 1. The exact answer for 4A, 4B, A^4, and B^4

Assume positive integers A<B, and set d=B-A. The user's identities are

    A+B=2A+d=2B-d,
    AB=A^2+Ad=B^2-Bd.                                  (1)

At exponentiation the corresponding exact relation is

    R(A,B):=A^B/B^A = A^d/(1+d/A)^A,
    log R(A,B)=d log A-A log(1+d/A).                    (2)

The extra factor in (2) is the correction missing from a linear analogy.
It determines the direction of the comparison.

There is also a direct transport between fourth powers of bases and
quadrupling exponents:

    (A^4)^B=A^(4B)=(A^B)^4,
    (B^4)^A=B^(4A)=(B^A)^4.                            (3)

Thus applying the same fourth power to the two compared values preserves
their equality and inequality. It is different from multiplying the bases
by four. More generally, for every positive real scale t,

    (tA)^B/(tB)^A = t^d R(A,B),
    R(tA,tB) = [t^d R(A,B)]^t.                         (4)

In particular R(4A,4B)=[4^d R(A,B)]^4. The scale correction can reverse
the comparison: 2^3<3^2, while 8^12>12^8. It destroys the integer tie:
R(2,4)=1 but R(8,16)=65536. These are exact integer comparisons.

### The equality curve and the particular doubles seam

Put r=B/A>1. Equation (2) gives

    A^B=B^A iff A^(r-1)=r,
    A=r^(1/(r-1)), B=r^(r/(r-1)).                       (5)

This is the positive-real parameterization in the user's
[linked equation page](https://en.wikipedia.org/wiki/Equation_xy_%3D_yx),
derived directly here and inherited with its rational classification
from THM-4107. Only (A,B)=(2,4) is a distinct positive-integer solution:
log(x)/x decreases strictly for x>=3; A=1 always loses; for A=2,
B=3 loses, B=4 ties, and 2^B>B^2 for B>=5. The last assertion starts
at 32>25 and propagates since 2B^2>(B+1)^2 for B>=3.

On a doubling edge the entire comparison simplifies to

    R(A,2A)=(A/2)^A.                                   (6)

So A=1 loses, A=2 ties, and every A>=3 wins. Along the column
1,2,4,8,..., the tie occurs exactly between 2 and4. Other odd-root
columns q,2q,4q,... with q>=3 have no such tie. This is a genuine
operation-specific connection to the first doubles seam, without a claim
that it causes either graph threshold.

## 2. A valid high-power simplification, scoped to comparison

Simultaneously replacing both arguments by their r-th powers gives

    log R(A^r,B^r)
      =r A^r B^r [log(A)/A^r-log(B)/B^r].               (7)

For every integer r>=2, the function log(x)/x^r decreases strictly on
[2,infinity), because its derivative is

    (1-r log x)/x^(r+1)<0.

Hence for positive integers A<B and every r>=2,

    (A^r)^(B^r) > (B^r)^(A^r) iff A>=2.                (8)

If A=1 the reverse inequality holds; there are no ties. **For this
specific orientation observable, all powers r>=2 give the same answer.**
The stabilization already occurs at the square level, not first at four.
It does not identify the magnitudes of the expressions or their prime
factorizations, remainder data, or Collatz routes.

### The primitive comparison and its surviving cycles

As in THM-4107, normalize each pair by g=gcd(A,B), obtaining coprime
p=A/g and q=B/g. Orient A->B when p^q>q^p. Unlike the raw scalar
preorder, this pair-dependent normalization can create directed cycles.
For A<B the smaller vertex loses exactly when A divides B or 3A=2B.

Common scaling preserves this primitive relation, because the common
factor cancels. Raising both inputs to the same integer power r>=2
instead removes the 2:3 exception and leaves precisely

    A^r -> B^r iff A does not divide B,       A<B.       (9)

This follows from gcd(A^r,B^r)=g^r and the same decreasing potential
for p,q>=2. If p=1, the larger vertex wins. It is exactly the higher
power stabilization proved in THM-4107, section4.

Odd cycles need not disappear: at r=4 there is the explicit primitive
cycle

    16 -> 81 -> 256 -> 16.

The last edge uses the reduced pair 1,16, whereas the first two use
coprime pairs. Raw exponent comparison has no directed cycle at all.
Thus the cycle comes from edge-specific gcd normalization, not merely
from exponentiation or oddness.

## 3. Oddness removes even terms, not all higher terms

Write A=c-h and B=c+h, with 0<h<c. Then (1) becomes

    A+B=2c,      AB=c^2-h^2.

Fourth powers split into an even part and an odd part in the signed gap:

    B^4+A^4=2c^4+12c^2 h^2+2h^4,
    B^4-A^4=8c^3 h+8c h^3.                             (10)

For this antisymmetric quartic polynomial, the h^4 term really cancels.
This is a valid version of the suggested simplification. Its hypotheses
are the polynomial degree and the chosen antisymmetric observable.

The exponent-swap logarithm is also odd in h, but it is not a polynomial:

    log R(c-h,c+h)
      =c log((c-h)/(c+h))+h log(c^2-h^2)
      =2h(log c-1)
       -sum_(j>=1) [(4j+1)/(j(2j+1))] h^(2j+1)/c^(2j). (11)

The expansion converges absolutely for |h|<c. It follows by expanding
log(1-h/c) and log(1+h/c); every even power cancels, but every higher
odd coefficient is nonzero. The first are -5/3, -9/10, -13/21.

There is a decisive boundary example: A=2,B=4 gives c=3,h=1 and
log R=0 exactly. Every finite truncation of (11), retaining its linear
term, gives a strictly positive answer, because every omitted term is
strictly negative. Thus even an arbitrarily long finite odd truncation
fails to recover this exact tie. A rigorous approximation needs a tail
bound, not just the fact that the function is odd.

### Odd tournament cycles likewise retain lengths five and above

[THM-002, odd-cycle collection formula](../../01-canon/theorems/THM-002-ocf.md)
expresses a tournament's Hamiltonian-path count by collections of disjoint
directed odd cycles. It includes all odd lengths, not just triangles.
The regular five-vertex tournament H_5 has exactly five directed triangles,
two directed five-cycles, and fifteen Hamiltonian paths. Since two
nontrivial odd cycles cannot be vertex-disjoint on five vertices,

    H(H_5)=1+2*5+2*2=15.

Deleting the length-five contribution would give11. The exact count is
independently checked by enumerating all120 vertex orders. This is also
the five-vertex quotient of the preceding canonical halving construction.
Odd degree in (11), odd cycle length here, and odd prime exponent in a
factorization are different involution-based selections; no equality
between these predicates has been supplied.

## 4. Fourth powers give exact factorization filters

The sum/product viewpoint supplies the Sophie Germain identity directly:

    A^4+4B^4=(A^2+2B^2)^2-(2AB)^2
      =(A^2-2AB+2B^2)(A^2+2AB+2B^2).                  (12)

For positive integers, the first factor is (A-B)^2+B^2. It is one
only at A=B=1, where the sum is5. Every other positive pair gives a
composite number. In particular p^4+4 is composite for every prime p.
This is a valid quartic exclusion for a **specified prime candidate**.
It does not exclude other remainders in the fourth-power shell.

The [Mills lane](seam_mills_20260925.md) explains the correct retained
coordinate: q=p^c+d, with the actual gap digit d. For example
347=7^3+4 is prime. A fourth-power factorization or a prime-fourth-power
divisor filter cannot license discarding all gap digits divisible by four.

## 5. What the other recovered decoders contribute

### Divisor balance: a finite controller with prime-support registers

The [divisor lane](seam_prime_balance_20260925.md) recovers exactly

    F=S+U iff N=p, p^3, p^2 q r,                       N>=2,

with pairwise distinct primes in the last family. Here F,S,U count
proper nontrivial divisors, squarefree divisors, and prime divisors;
the values of F are0,2,10. They are not tournament orders. The example
60=2^2*3*5 has10 proper nontrivial divisors, seven squarefree and three
prime. Its total prime degree is four and cannot be omitted.

A number divisible by any prime fourth power can be safely rejected for
this equation. The stronger useful result is an exact controller with
eight live exponent profiles and one absorbing reject state under prime
multiplication. It retains at most three named primes while live. Names
are necessary:6 and15 share profile(1,1), but multiplication by four
gives unbalanced24 and balanced60. Rejected16 can halve to balanced8;
therefore this multiplication controller is not a Collatz quotient.

### Thresholds: the same ear operation with different attachments

The [threshold lane](seam_threshold_20260925.md) identifies the common
operation behind the nearby numbers. Q_14 is a tree with three leaves;
the ear1--15--10 repairs its third-leaf obstruction. In the sparse
operation graph, G_15 remains planar; the ear8--16--14 completes its
K_(3,3) subdivision. Thus Q_15 already succeeds while G_15 remains planar.
The shared degree-two insertion is an exact structural analogy. Its
attachments decide whether it repairs Hamiltonicity or breaks planarity.

### Triangular ten: a support-preserving map to regular five

T_4=10 counts the ten edges of K_5. Treat those edges as vertices of
L(K_5). A regular five-tournament pairs each vertex's two outgoing arcs,
giving five connected pairs among the ten edge objects. With their
support labels and common centers retained, this pairing decodes the
original tournament exactly and contracts to K_5. There are144 incident
edge pairings, precisely24 with distinct centers, corresponding to all
24 regular five-tournaments.

This gives a useful concrete relation to the desired H_5 quotient. The
earlier ten-vertex object came from3*3+1; it is not automatically this
edge-object model. Indeed K_5 cannot partition its edges into three
triangles and one extra edge, by vertex-degree parity. The different
roles of ten require the declared support map.

After these connections the board has two honest high-layer reductions:
the power-comparison orientation stabilizes, and the divisor classifier
has an absorbing multiplicative reject region. It also has three reasons
to retain more data: log-gap tails, Mills gap digits, and graph attachment
or owner labels. None of those reductions deletes arbitrary higher
states from a Collatz trajectory.

## 6. Reproduction and audit scope

    python 04-computation/experiments/seam_power_coordinates_20260925.py
    python -O 04-computation/experiments/seam_power_coordinates_20260925.py

[Recorded output](seam_power_coordinates_20260925.out). All checks use
exact integers or rational arithmetic; no floating-point sign decisions.
Universes: all2016 pairs1<=A<B<=64 for the coordinate/scale identities,
integer comparison and square-level controls; all66 pairs through12 for
literal fourth-power comparison; all4096 positive pairs through64 for
the factorization; formal log coefficients through degree21; every
Hamiltonian vertex order and odd directed cycle of the fixed H_5.
No random sampling or hidden filters. General proofs above establish
the all-integer claims independently of the finite universes.

Independent audit checks the gap expansion, its tail at the exact tie,
the primitive/raw distinction, and the concrete H_5 count. The other
three lanes contain their own proofs, controls, and primary source
boundaries. Universal Collatz root generation remains OPEN.
