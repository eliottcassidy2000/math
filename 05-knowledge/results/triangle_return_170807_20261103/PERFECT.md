# A complete inverse for one missing prime-power factor

2026-11-03. **PROVED elementary compiler; CITED historical examples;
FINITE-EXACT controls.** Odd-perfect existence is OPEN. The restricted
exclusion below is not an improvement to known global factor-count bounds.

## PC1. Fixed-core completion removes the prime and exponent search

Fix a positive integer C. Seek a prime q not dividing C and an exponent
a>=1 such that C q^a is perfect. Set A=sigma(C), D=2C-A.
If D<=0 there is no completion, since sigma(q^a)/q^a>1. Otherwise any
completion necessarily has

    q=ceil(A/D),  B=A-D(q-1)>0,  q^a=A/B.                (PC1)

Conversely these formulas give a perfect completion exactly when q>=2
is prime, gcd(C,q)=1, and A/B is an integral positive power of q.
Thus a fixed core has at most one prime-power completion, without a bound
on q or a. For an odd-perfect search require additionally that C is an odd
square and q=a=1 modulo4, as in Euler's classical form.

Proof. Multiplicativity gives D/A=1/q+...+1/q^a. Hence

    1/q <= D/A < 1/(q-1),  so A/D <= q < A/D+1.

This interval has exactly the displayed integer candidate. Summing the
geometric series gives q^a[A-D(q-1)]=A. All reverse implications follow
from the same equation with the prime and coprimality guards present.
Uniqueness is exact; no approximate logarithm is needed for the power test.

The source still includes the factorization of C. The pair(A,D) computes
the candidate but does not itself certify square shape or coprimality.
The unbounded choice of C is the remaining odd-perfect problem, not a
finite computational corollary of PC1.

## PC2. A source-preserving square-factor recurrence

Adjoin p^2, with p prime and p not dividing C. The exact update is

    A'=(p^2+p+1)A,
    D'=p^2 D-(p+1)A.                                   (PC2)

Its integer matrix has determinant p^2(p^2+p+1). Retain the prime support,
exponents and this denominator loss along with(A,D). In particular this
clock is not an invertible field action. It is an actual multiplication
operation on guarded cores, extending the repository's low-term divisor
grammar and [perfect-core affine law](../perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md).

For the terminal exponent a=1 one can bound the next square prime without
an arbitrary cutoff. If p does not divide A, reduction of

    (qD-A)p^2=A(q+1)(p+1)

shows p^2 divides q+1. For odd q,p this gives q>=2p^2-1. Substitution yields

    2D p^2-2A p-(3A+D)<=0.                              (PC3)

The exceptional primes p dividing A form an explicit finite list. PC3
therefore compiles the entire C p^2 q family into finitely many candidates.
This is a valuation argument; an unlabeled cycle in a factor graph would
lose the squared divisibility and its factor2 parity constraint.

## PC4. The11/13 seam and an all-exponent example

Take C0=3^2*7^2*11^2=53361. Then

    A=98553, D=8169, A/D=4693/389, 12<A/D<13.

There is no odd perfect number C0 p^2 q^a for distinct primes p,q outside
the core and Euler exponent a=1 modulo4. Here is a complete finite reduction
covering unbounded p and a, independent of any assumed global prime bound.

PC1 on C0 p^2 forces q>=13. If q=13, the ratio at exponent5 already exceeds
2C0/A, while the desired ratio is smaller than2C0/A. Thus only a=1 is
possible. Its quadratic in p is

    (26C0-14A)p^2=14A(p+1).

The discriminant is not a square (integer certificate in the receipt).
For q>=17 the prime-power abundancy is strictly below17/16. Necessarily

    (32C0-17A)p^2 <17A(p+1).

The coefficient on the left is32151>0; the inequality fails at54 and is
strictly increasing thereafter. Hence p<=53. PC1 rejects each of the finite
prime candidates, with its precise failed guard retained in the receipt.

At p=13=12+1 the forced terminal value is22021=19^2*61. It balances the
formal equation but is composite. This is **Descartes's classical spoof**,
not a newly discovered candidate. The factor transmissions are literal:

    sigma(3^2)=13,       sigma(7^2)=3*19,
    sigma(11^2)=7*19,    sigma(13^2)=3*61,
    22021+1=2*7*11^2*13.

Thus the old12-1/12+1 inspiration becomes the actual11/13 factor ports,
with multiplicities retained. Replacing the formal22022 by its true divisor
sum23622 gives the exact defect

    sigma(C0*13^2*22021)-2*C0*13^2*22021
       =1600*sigma(C0*13^2)>0.

The positive completion map is meaningful precisely because it rejects
this formal balance. The literature already studies spoof factorizations
and demonstrates the limitations of multiplicativity-only approaches:
[BYU Computational Number Theory Group, Odd, spoof perfect factorizations](https://arxiv.org/abs/2006.10697).
Euler's form and Touchard's12/36 congruences are recalled in
[Starni, Some Extensions to Touchard's Theorem](https://arxiv.org/abs/1709.05286).
Neither citation is used as a claimed new theorem of this packet.

## A usable synthesis and next problem

The tournament and arithmetic operations share a precise design rule:
positive realizations yield restrictive inverse tests. For births, the
positive directed cut bounds all future prefix counts; for perfect cores,
the positive geometric tail traps q in an interval of width one. Retain
the source relations/factorizations before applying either inverse.

The arithmetic target worth pursuing is a uniform bound on the depth of
prime-square extensions that can keep passing PC1-compatible completion
guards, or a structural obstruction to all such infinite branch families.
PC2 alone does not supply that bound. The tournament analogue is to find a
source grammar stable under both ternary returns, rather than repeat a
finite successful output bank. This is an algorithmic analogy between two
guarded maps; no bijection of perfect numbers and tournaments is asserted.

## Incoming integration: eliminate a row after fixing support

During this session, [Session178 prime-flow](../session178_prime_flow_honeycomb_20261104/primeflow/PROOF.md)
proved an effective exponent box for a fixed prime support and Euler prime,
using LTE capacities, LCM conjunctions and exact cut balances. This is
complementary to PC1, not a result of this packet. In its finite decision
procedure, choose the non-Euler exponents first and form their actual core
C. PC1 then either rejects the choice or determines the Euler prime and
its exponent uniquely. There is no need to enumerate the final exponent
row. All earlier support, valuation, LCM and outside-prime guards still
apply. This is an exact variable elimination; no benchmark advantage over
the incoming early-rejection sieve has been established.

The incoming [quotient-debt note](../session178_prime_flow_honeycomb_20261104/unification/NOTE.md)
also distinguishes arithmetic balance from an actual graph circulation.
PC1 asserts sigma(N)=2N, not that prime-factor arrows form a closed walk.
The Descartes reaction list above must retain its prime realization guard.
