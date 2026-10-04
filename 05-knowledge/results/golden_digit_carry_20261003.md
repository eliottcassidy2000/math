# Literal golden digits, the number 105, and reversible carry memory

2026-10-03 (America/Denver).

**Status:** PROVED elementary pair-reader, polynomial-kernel, carry-composition,
and coordinate-change identities. FINITE-EXACT for the explicit computational
universes below. Classical base-phi numeration is CITED, with no novelty claim.
This note supplies arithmetic storage and verification, not a primality
criterion, a Collatz convergence argument, or a bounded-state normalizer.

## Inheritance and the three digit types

The closest established mechanism is `phi^2=phi+1`, together with the
two-coordinate Fibonacci reader in
[the guarded Zeckendorf note](zeckendorf_guard_automaton_20261003.md).
The canonical hostile is the missing radix: `10.01_phi=2`, whereas the
concatenated word `1001_phi=2+2phi` and the Zeckendorf word `1001` represents 6.
The corrected near miss is treating the same binary string in different
numeration systems as the same integer. The least-used sidecars are the
radix position and the signed polynomial discarded by normalization.

The owner's earlier golden prompt was recorded in
[the golden holonomy note](collatz_golden_holonomy_20261001.md), and its
current proved route is
[THM-4528, Collatz parity in base phi](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md).
That route reads an integer's **orbit parity itinerary** as a base-phi
fraction. It does not read the **source integer's positional digits**.
MISTAKE-556 corrects its dynamical wording to semiconjugacy and records the
exceptional endpoint. No parity-to-source-digit identification is used here.

The live objects are literal Laurent digits, integer Fibonacci digits, orbit
parity digits, integer residue guards, and carry witnesses. Their readers
share the quadratic recurrence, while their input meanings remain distinct.
The confirmed red/black/blue construction uses marked extra unit atoms; its
markers are representation coordinates, not an identification with the
auxiliary four-state charge or with literal powers of phi.

**CITED background.** Bergman's
[A Number System with an Irrational Base (1957)](https://math.berkeley.edu/~gbergman/papers/base_tau.pdf)
gives the local `100=011` rule, the examples `2=10.01` and `3=100.01`, and
arithmetic with overlapping digits. Shallit's
[Proving Properties of phi-Representations with the Walnut Theorem-Prover](https://arxiv.org/abs/2305.02672)
provides a modern automata treatment. These are primary sources for the
classical setting; the finite greedy implementation below does not establish
a new global termination theorem.

## Exact pair reader and radix memory

Work in `Z[phi]=Z[t]/(f)`, where `f=t^2-t-1`. Represent `A+Bphi` by `(A,B)`.
Multiplication and multiplication by phi are

    (A,B)(C,D) = (AC+BD, AD+BC+BD),
    phi(A,B) = (B,A+B).

Since `phi^-1=phi-1`, the inverse shift is `(A,B)->(B-A,A)`. Appending a digit
`d` on the right sends `(A,B)->(B+d,A+B)`. Thus Horner reading is exact, uses
two integer registers, and requires no floating point.

If a positional word has `R` digits after its radix, concatenating all digits
first reads `phi^R` times its value. Apply the inverse shift `R` times to
recover the literal value. A literal integer has exact pair `(n,0)`.
The vanishing of the second coordinate only modulo one modulus does not
certify literal integrality.

All these operations reduce modulo every positive integer `m`; phi is a unit
even when 2 or 3 is not. For a fixed `m`, the pair reader is finite-state.
The radix correction may be stored modulo the order of phi in
`(Z/mZ)[phi]`; the exact integer reader still requires its actual radix.
The order and prime-splitting analysis is in
[the companion clock note](golden_prime_clocks_20261003.md).

This is a faithful connection from a literal word with radix to its
quadratic-ring value and then its residue. It preserves arithmetic value;
it loses the chosen representation and event history. The radix and the
carry certificate below repair specific losses. An ordered factor list is
needed if the original factorization is part of the object.

## The literal calculation for 105

Exact greedy reading gives

    3 = 100.01_phi       = phi^2+phi^-2,
    5 = 1000.1001_phi    = phi^3+phi^-1+phi^-4,
    7 = 10000.0001_phi   = phi^4+phi^-4,
    105 = 1001010101.0101001001_phi.

The normal word's occupied exponents are

    9, 6, 4, 2, 0, -2, -4, -7, -10.

The raw product `C(t)=sigma(3)sigma(5)sigma(7)` is

    t^9+2t^5+t^2+2t+ t^-2+2t^-3+t^-6+t^-7+t^-10.

There are repeated digits at exponents 5, 1, and -3. Discarding their
multiplicities before normalization changes the value. Let `B(t)` be the
normal Laurent polynomial displayed above. With `R=10`, exact long division
gives

    t^10(C-B) = (t^2-t-1) Q,
    Q = -t^14+t^13-t^12-t^10+t^9-t^8-t^6+t^5-t^4
      = -t^4(t^2-t+1)(t^8+t^4+1)
      = -t^4 Phi_3(t) Phi_6(t)^2 Phi_12(t).

Here `Phi_j` denotes the cyclotomic polynomial. The last factorization is an
exact finite feature of this carry polynomial, not a primality test. It
exhibits roots of unity in the record of normalization; the value relation
itself is evaluated at the real root phi of `f`, which is not a root of unity.
Indeed `Phi_3(phi)=2phi^2`, `Phi_6(phi)=2`, and `Phi_12(phi)=2phi^2`, so
`Q(phi)=-16phi^8=-208-336phi`. Its quadratic norm is 256, using
`Norm(A+Bphi)=A^2+AB-B^2`. This carry specialization does not recover the
original factors 3, 5, and 7; the factor parse is a separate witness.

The owner's identities are the Lucas identity

    phi^(2k)+(-1)^k = L_k phi^k,
    L_k = phi^k+(-1)^k phi^-k.

It follows by multiplying the second equality by `phi^k`; the second follows
from the two conjugate roots of `t^2-t-1`. The instances `L_4=7` and `L_5=11`
are followed by `L_6=18`. In particular, the sparse canonical word
`18=1000000.000001_phi=phi^6+phi^-6` is composite. Binary digits or sparse
Lucas formulas alone therefore do not distinguish primes.

## Unique carry certificates and composition

**PROVED.** Let `C,B` be finite integer-coefficient Laurent polynomials with
`C(phi)=B(phi)`. Choose a fixed nonnegative `R` clearing the negative powers
of both. There is a unique `Q in Z[t]` such that

    t^R(C-B)=fQ.

Proof: division by the monic polynomial `f` leaves an integer linear
remainder `a+bt`. Evaluating at phi gives `a+bphi=0`; irrationality forces
`a=b=0`. Uniqueness follows because `Z[t]` is an integral domain. Conversely,
the displayed identity implies equal values at phi. Thus the certificate
provides an iff check, with nonzero remainders rejecting false carries.

The tuple `(B,Q,R)` recovers `C=B+t^-R fQ` exactly. It does **not** recover a
particular chronological rewrite sequence or the factorization of `C`.
Store the rewrite/operation tree when those distinctions matter. The
certificate is an aggregate signed carry, with no bounded-coefficient or
online-normalization claim.

Two local operations have explicit certificates at every integer offset k:

    t^k+t^(k+1)-t^(k+2) = -t^k f,
    2t^k-t^(k+1)-t^(k-2) = t^(k-2)(1-t) f.

The first is `11 -> 100`; the second is the overlap repair
`2phi^k=phi^(k+1)+phi^(k-2)`. A lone binary merge rule is insufficient to
describe convolution because its digits may exceed one.

Let `sigma(n)` be a selected normal Laurent representative of a nonnegative
integer. Define the radix-independent Laurent carry functions

    A(a,b) = [sigma(a)+sigma(b)-sigma(a+b)]/f,
    C(a,b) = [sigma(a)sigma(b)-sigma(ab)]/f.

They satisfy

    A(a,b)+A(a+b,c) = A(b,c)+A(a,b+c),
    C(ab,c)+C(a,b)sigma(c) = C(a,bc)+sigma(a)C(b,c).

Both equalities follow by expanding their numerators and using ordinary
associativity, then canceling the nonzero `f`. They allow normalization
certificates to compose through different arithmetic parse trees. They do
not assert that the parse trees themselves are recoverable after quotienting.

## The exact bridge to Fibonacci digits, and its arithmetic defect

For a binary word on Fibonacci weights `1,2,3,5,...`, read the same digits
without a radix in base phi, obtaining `(A,B)`. Its Fibonacci value `X` and
one-place Fibonacci shift `Y` are

    (X,Y) = (A+2B, 2A+3B),
    (A,B) = (2Y-3X, 2X-Y).

The matrix is `M^3`, with `M=[[0,1],[1,1]]`, and determinant -1. It is an
invertible integral change of coordinates, including after reduction modulo
any integer. It conjugates the positional append rule to

    (X,Y) -> (Y+d, X+Y+2d),

exactly the inherited two-register Fibonacci reader. This is the faithful
bridge to an arithmetic guard: retain both registers before projecting to
the integer residue. It does not identify the literal value `A+Bphi` with X.

For the canonical Fibonacci digits of n, set
`b(n)=floor((n+1)/phi^2)`. The inherited Beatty-register formula gives

    L(n)=n-2b(n)+b(n)phi.

The missing coordinate is necessary for transporting arithmetic. With
`delta=b(a)+b(b)-b(a+b)`, direct substitution gives

    L(a)+L(b)-L(a+b) = delta(phi-2) = -delta phi^-2.

For the projection `ell(A+Bphi)=A+2B`, multiplication instead obeys

    ell((A+Bphi)(C+Dphi)) = ell(A+Bphi)ell(C+Dphi)-BD.

The smallest useful hostile is `L(2)=phi` and `L(2)^2=phi^2=L(3)`, although
`2*2=4`. Therefore this lift is a recurrence-preserving coordinate map, not
an embedding of the ordinary integer ring. The literal positional expansion
of n, whose pair is `(n,0)`, does preserve integer multiplication.

## Collatz interface and exact audit

For a literal integer word `sigma(n)`, the numerator `3n+1` is obtained before
normalization as

    (t^2+t^-2)sigma(n)+1.

This follows from `3=phi^2+phi^-2`, and retains all overlap multiplicities.
The carry certificate verifies the subsequent normalization. Division by
`2^k` still requires the actual integer congruence and exact valuation guard;
being able to shift powers of phi does not supply it.

For a finite Collatz parity prefix with ones at positions j and length L,
the same pair API evaluates `sum phi^(-j-1)` and the multiplier `phi^-L`.
The source/endpoint identity then uses the endpoint itinerary as a separate
tail. This is an interface between two typed readers, not an equality of
the source integer's digits with its orbit digits.

Reproduction:

    python3 04-computation/experiments/golden_digit_carry_20261003.py

The saved [exact output](golden_digit_carry_20261003.out) records:

- 1,001 exact greedy normal forms for integers 0 through 1,000;
- 13,013 independent Laurent-sum/Horner modular comparisons over the thirteen
  stated moduli, with no primality filter;
- the 105 convolution, recovery, and cyclotomic factorization above;
- 82 local positive/negative-offset rewrite controls and 900 product carries
  for factors 1 through 30, plus an unequal-value rejection;
- 1,000 additive and 1,000 multiplicative carry-composition identities for
  triples with entries 1 through 10;
- 2,002 reader-append coordinate controls, and 10,201 checks each of the
  additive and multiplicative Zeckendorf-lift defects;
- Lucas identities at k=1 through 20 and the sparse composite 18;
- 333 exact Collatz-numerator normalization certificates.

All arithmetic is integral. The greedy search has a declared finite
fractional-digit bound and reports failure if exhausted; its successful
finite universe is not promoted to a general algorithmic complexity bound.
The general pair and carry statements are justified by the proofs above.
