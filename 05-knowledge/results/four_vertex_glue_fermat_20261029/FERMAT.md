# Fermat completeness: exact squaring recurrences and the remaining obligation

2026-10-29. **PROVED elementary formulas / CITED classical dynamics /
FINITE-EXACT controls. Exactly five Fermat primes remains OPEN.**

## Recovery before derivation

The source is [S58](../s58_fermat_pepin_collatz_20261026/README.md), including
its [independent audit](../s58_fermat_pepin_collatz_20261026/AUDIT.md) and
MISTAKE-673. Its classical squaring criterion, Collatz/Pepin reformulation,
and obstruction to a finite bank of prime divisors are inherited. We do
not use squarefreeness of Fermat numbers, which is not known generally.
Nor is an affine mirror identity a certificate that two points occupy
the same orbit. The useful new integration here is an explicit CRT
counting interface retaining tails, periods, and relative phase, plus a
small odd-order certificate on the first composite Fermat number.

Primary sources consulted:

- [Somer and Krizek (2004), On a connection of number theory with graph
  theory](https://dml.cz/handle/10338.dmlcz/127904).
- [Somer and Krizek (2011), The structure of digraphs associated with the
  congruence x^k = y modulo n](https://dml.cz/handle/10338.dmlcz/141538).
- [Boklan and Conway (2016), Expect at most one billionth of a new Fermat
  Prime!](https://arxiv.org/html/1605.01371v2). Their small expected count
  is an explicitly model-dependent heuristic, not a finiteness proof.

Only 3,5,17,257,65537 are known Fermat primes. None of the computations
here establishes that there are no others.

## Exact formulas with their quantifiers

Let N>1 be odd, with its actual prime-power factorization N=product q_i^a_i.
On the full unit group, iterate S(x)=x^2. For t>=0 and r>=1 let C_N(t,r)
count units whose tail has length at most t and whose eventual period
divides r. Then

    C_N(t,r) = product_i gcd(phi(q_i^a_i), 2^t(2^r-1)).     (F1)

Proof: S^(t+r)(x)=S^t(x) is x^[2^t(2^r-1)]=1. The units modulo an odd
prime power form a cyclic group, where the number of roots of x^h=1 is
gcd(h,group order). CRT multiplies the counts. No squarefree assumption
is made. This is an elementary derivation of a classical power-digraph
formula, not a claimed new primality theorem.

Write phi(q_i^a_i)=2^alpha_i d_i with d_i odd. Then

    C_N(t+1,r) = 2^#{i:alpha_i>t} C_N(t,r).                (F2)

The number with tail <=t and exact period r is

    E_N(t,r) = sum_(d|r) mu(r/d) C_N(t,d).                 (F3)

Subtract E_N(t-1,r) for exact tail t>0; the number of actual r-cycles is
E_N(0,r)/r. These formulas give complete tail/period statistics from the
prime-power groups. They do not retain the identity of a designated
source; that still needs its CRT coordinates.

Two specific local cycles of lengths a,b combine into gcd(a,b) cycles of
length lcm(a,b). Their relative phase is indispensable. For example the
local cycles (2,4) modulo7 and (3,9) modulo13 produce two cycles modulo91:

    16 -> 74 -> 16,       9 -> 81 -> 9.                   (F4)

Both have the same two local unbased cycle types. Forgetting their phase
would identify distinct marked return orbits. This is a precise analogue
of keeping the C5 action on the gluing quotient, not an identification of
the respective groups or dynamical laws.

## What a proof of five primes would still have to construct

Put F_n=2^(2^n)+1. A nonidentity unit x with x^d=1 modulo F_n for an odd
d>1 certifies compositeness: a prime F_n has unit group of pure two-power
order. Conversely, if F_n is composite, such an odd-order unit exists.
Here is a direct proof of the converse, also giving the cited squaring
criterion without assuming squarefreeness.

If the whole unit group were a two-group, no repeated odd prime could
divide F_n, and every prime q dividing it would have q-1 a power of2.
Such a prime is a Fermat prime F_j. But Fermat numbers are pairwise
coprime, so F_j cannot divide F_n unless j=n, forcing F_n itself prime.
Thus a composite F_n has a unit group with an odd prime divisor of its
order; Cauchy's theorem supplies an odd-order unit. Under squaring it
lies on a nontrivial cycle.

The checker supplies the factor-independent verification of the witness

    N=F_5=4294967297,
    x=3303966555,       x != 1,       x^3 = 1 mod N,
    x -> 3041328343 -> x under squaring.                  (F5)

The construction used Euler's factorization 641*6700417 and CRT. The
certificate itself needs only two modular powers. This is a concrete
control for a uniform proof target: construct such a witness, or another
compositeness certificate, for every n>=5 without presupposing the desired
compositeness. Doing that uniformly is unresolved here; reformulating
the obligation as cycle existence does not solve it.

Two tempting shortcuts fail for exact reasons. First, the distinguished
source2 always obeys

    S^n(2)=-1,       S^(n+1)(2)=1 modulo F_n,             (F6)

for prime and composite F_n alike. Its short tree sees no odd-order
coordinate. Second, a fixed finite divisor bank cannot cover all n:
pairwise coprimality lets each prime divide at most one Fermat number.
The recurrence F_(n+1)=F_n^2-2F_n+2 does not supply a ring quotient between
successive moduli; they are coprime. A valid recursive certificate must
introduce new arithmetic information, not merely transport an old residue.

Verification covers every unit for every odd modulus3..151, t=0..5 and
r=1..12 (5400 tail/period parameter probes, with separate Mobius checks),
the relative-phase pair (F4), and the literal witness (F5). The positive
prime controls and composite controls for (F6) run through n=10. No
probability estimate or finite prime search is substituted for the
universal quantifier n>=5.
