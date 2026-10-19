# Exact divisor lattices in the cofactor map's inverse basins

19 October 2026. **PROVED, with independent audit PASS.** Extension of the global cofactor classification. This concerns F(n)=n-p+2c-1, where n=cp, c=2^a is the exact dyadic factor and p is odd. An aliquot interpretation additionally requires p to be an odd prime at every relevant step.

## 1. The inverse factorization

The forward formula is equivalent to

    F(cp)-1=(c-1)(p+2).

Thus a positive even integer m has one preimage with core c=2^a exactly when

    2^a-1 divides m-1,
    (m-1)/(2^a-1)>=3.

That preimage is

    n=2^a[(m-1)/(2^a-1)-2].

The quotient is odd, so its difference with 2 is positive odd and the proposed c is its exact dyadic factor. Different exponents give distinct preimages. This is an exact inverse compiler using Mersenne divisors, not an arbitrary branch selection.

## 2. Every immediate preimage of a fixed point

Fix A>=1 and the formal fixed point

    m_A=2^A(2^(A+1)-1).

**Theorem.** The immediate F-preimages of m_A are in bijection with the positive divisors a of A. Explicitly,

    n_a=2^a[(m_A-1)/(2^a-1)-2],     a|A.

There are exactly tau(A) immediate preimages, including the fixed point itself at a=A.

**Proof.** Sufficiency follows from

    m_A-1=(2^A-1)(2^(A+1)+1).

If a|A, the first factor is divisible by 2^a-1, and the quotient is at least 2^(A+1)+1>=5. Thus the inverse compiler applies.

For necessity, reduce A modulo a and write r with 0<=r<a. Modulo M=2^a-1,

    m_A-1 ≡ 2^(2r+1)-2^r-1.

Suppose r>0. If 2r+1<a, this integer is positive and strictly below M. Otherwise reduce 2^(2r+1) once using 2^a≡1; the residue becomes

    2^(2r+1-a)-2^r-1.

Here 0<=2r+1-a<=r, so the residue is strictly negative and strictly larger than -M. In either case it is nonzero modulo M. Hence divisibility forces r=0, that is a|A. The case a=1 already has r=0. ∎

This necessity uses only elementary modular arithmetic; no primitive-prime-divisor theorem is needed.

## 3. Native perfect-number basins

When 2^(A+1)-1 is prime, m_A is an even perfect number. Its immediate predecessors that remain in the dyadic-times-odd-prime class are obtained by testing the finitely many divisor-indexed odd values

    p_a=(m_A-1)/(2^a-1)-2,     a|A.

The a=1 alternative is never a new native predecessor for A>=2, since

    p_1=m_A-3=(2^A+1)(2^(A+1)-3)

is composite. At A=1 it is the fixed point 6 itself.

For example m_4=496 has divisor indices 1,2,4. Their odd parts are 493,163,31. Thus the native immediate preimages are exactly 652=4*163 and 496=16*31, while the core-2 preimage 986=2*493 is not native.

The theorem classifies immediate predecessors, not the entire inverse tree. Repeated application of the general compiler describes that tree, but subsequent levels need not have divisor indices of the same A. Nor does it establish infinitely many Mersenne primes or infinitely many nontrivial native predecessors.

## 4. Every immediate preimage of a power of two

**Theorem.** For b>=1, the positive even solutions to F(n)=2^b are in bijection with the proper positive divisors a of b:

    n=2^a[(2^b-1)/(2^a-1)-2],     a|b, a<b.

Thus 2^b has exactly tau(b)-1 immediate preimages.

The identity 2^a-1 | 2^b-1 iff a|b follows by Euclidean division of the exponents. At a=b the quotient is one and the odd part would be -1, so it is excluded. At a proper divisor b/a>=2 the quotient is at least 1+2^a>=3, giving a positive odd part. This proves the theorem, including b=1, which has no preimages.

For b=2t, the divisor a=t gives n=2^t(2^t-1), the first resonance in the descent proof. For b=3t, the divisor a=t gives n=2^t(2^(2t)+2^t-1), the intermediate value reached from the second resonance n_0=2^(2t)(2^t-1). Thus the resonance families behind the sharp transient-height law sit inside exact Mersenne-divisor inverse trees.

## 5. Structural lesson

The forward proof uses valuation changes and decreasing checkpoints. The backward proof uses divisibility of m-1 by Mersenne numbers. These are two exact descriptions of the same map. Their overlap explains why powers of two and formal perfect-number points organize the exceptional trajectories. The appearance of 2^(A+1)+1 in fixed-point preimages is an algebraic factorization; it is not a claim that this factor is prime or that only the five known Fermat primes exist.
