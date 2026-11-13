# Incoming valuation pinning and the joint branch family

**PROVED IN NOTE**, extending an incoming theorem; independent review is [AUDIT-joint-refinement](../AUDIT-joint-refinement.md). The source is [Session292 arithmetic, section3](../../session292_triple_bridges_20261113/arithmetic/PROOF.md), commit8e0debd6b2, fetched at the first checkpoint. That session independently obtained the same minimal single-count wheel and then strengthened its hidden-factor theorem. The common wheel is credited as concurrent recovery; no priority claim is made.

Let A=99*2^k+23, B=3A+2, C=3A+4, k>=1. All three are odd, and

```
gcd(A,B)=gcd(A,2)=1,
gcd(A,C)=gcd(A,4)=1,
gcd(B,C)=gcd(B,2)=1.
```

Thus their prime divisors are pairwise disjoint at every exponent, independent of any sieve.

## Arbitrary hidden factor complexity in every accepted joint phase

Fix any accepted residue class of the joint first25-prime wheel and any three nonnegative integers h_A,h_B,h_C. There is an infinite arithmetic subprogression on which A,B,C have at least h_A,h_B,h_C distinct prime factors, respectively. All these primes exceed97 and the three lists are pairwise disjoint.

Proof. The incoming valuation-pinning theorem applies to each positive exponential a*2^k+c with c odd on any infinite progression. Here is its elementary mechanism. Suppose finitely many primes S are already forced. At one exponent k0, record each exact valuation e_p. Refine the exponent progression by all conditions k=k0 mod ord_(p^(e_p+1))(2). They are compatible because k0 satisfies all of them. This fixes each S-valuation exactly. Since the exponential grows without bound on the refined progression, eventually it exceeds the product of these fixed prime powers and has a new prime factor q outside S. Pin q by its multiplicative order and repeat. An already accepted wheel class excludes all primes<=97 throughout every refinement.

Apply this argument h_A times to A, then h_B times to B, then h_C times to C, always retaining the previously imposed congruences. Every intermediate progression is infinite and contains its chosen witness. Earlier divisibility conditions persist. The gcd identities above force disjointness of the three lists. This proves the simultaneous theorem. It does not produce prime values or a tournament-return construction.

The explicit all-composite progression in [J3](JOINT_BRANCH_SIEVE.md) is a small certificate for one forced factor in each coordinate; the new statement quantifies over arbitrary complexity in every accepted joint phase.

## The entire formal triple misses A002496

The earlier single-count exclusion extends to both formal branches. For k>=2, A and B are3 mod4, so neither is a square plus1. For k>=4,

`C-1=297*2^k+72 = 8 mod16`,

which is not a square residue modulo16. The remaining positive exponents are exact finite boundaries:

* A_1-1=220;
* B_1-1=664;
* C_1-1=666, C_2-1=1260, C_3-1=2448.

None is a square. Consequently no entry of (A_k,B_k,C_k) belongs to A002496 at any k>=1, even when an entry is prime. This is an exact family classification, not a claim about the infinitude or enumeration range of A002496.

For A135244 the [inherited sharp bound](JOINT_BRANCH_SIEVE.md#j4-an-inherited-sharp-aliquot-interface) transports the subprogressions with h_A,h_B,h_C>=2 into strict inequalities below the corresponding prime-square extrema. The map preserves that extremal predicate; the finite prime wheel alone does not.

## Shared-variable comparison

The valuation refinement preserves a common exponent, just as the dihedral CRT construction preserves a common sign and the square-kernel birth preserves a named component partition. These are three different compatibility problems. Their common successful move is to refine the retained state before composing constraints. No equality of their state spaces, spectra, or prime labels is asserted.
