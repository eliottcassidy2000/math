---
id: THM-5105
title: Two-burst shifted-product inverse and exact dyadic coverage
status: PROVED all orders within the specified physical language; independently audited, with separate finite exact tests. No global coverage or historical-priority claim.
source: 05-knowledge/results/session241_prime_relations_20261110/coverage/PROOF.md
session: Session241, 2026-11-10
depends_on:
  - THM-4888-finite-oz-banks-mod3-law-involution-rigidity-residue-confinement
  - THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family
---

# THM-5105: exactly invertible, exactly thin two-burst arithmetic

For a,b>=1, form a transitive spine with root-edge word **O O^a Z O^b Z**.
The root beats O vertices and loses to Z vertices. The initial O is essential;
it is the rooted-edge seed, not part of either displayed burst. This tournament
T(a,b) is strong and has order a+b+4.

With r=a+2>=3 and s=b+1>=2,

    H(a,b)=1+2^(a+1)+2^(a+b+2)−2^b,
    2H(a,b)−3=(2^r−1)(2^s+1).

The map(a,b)→H is injective and has an exact factorization-free inverse:
inspect the binary runs of2H−3, decode(r,s), and verify the domain and product.
The shifted product has one all-ones run for r=s, two equal r-runs for r<s,
and leading run1 with trailing run s and total length r+s+1 for r>s.
These descriptions are disjoint in the stated domain.

For every integer k>=4, the number of distinct H in this family at most2^k is

    (k−4)(k−3)/2 + floor(k/2)−2.

Pairs with r+s<=k qualify, those with r+s>=k+2 do not, and the boundary
r+s=k+1 qualifies exactly when r<s. Hence the full image count up to X is
(1/2)(log_2 X)^2+O(log X), of density zero even among odd numbers. This is
an all-size obstruction for this language, not an obstruction to Chan–Pak.

The product coordinate retains precise arithmetic sidecars. For odd prime p
with t=ord_p(2), divisibility of2H−3 by p is equivalent to t|r or to t even
and s=t/2 modulo t. Odd-prime LTE gives both factor valuations without assuming
non-Wieferich behavior. Moreover gcd(2^r−1,2^s+1)=2^g+1 if r/g is even,
and1 otherwise, where g=gcd(r,s). These are not primality tests for H.

The inherited O/Z calculus proves the count; the independent audit also
derives it by placing the root between directed spine runs. For a fixed physical seed(D,E,U,d), appending nonempty O^a Z O^b Z
gives2(H−D)−d=(2^(a+1)U−d)(2^(b+1)+1), yielding a divisor-based exact
inverse test without asserting injectivity for arbitrary seeds.

[Full proof and source paths](../../05-knowledge/results/session241_prime_relations_20261110/coverage/PROOF.md)
include physical conventions, the binary inverse, precise counting boundary,
prime-power sidecars, and positive/hostile controls.
[Independent proof and exact audit](../../05-knowledge/results/session241_prime_relations_20261110/coverage_audit/REVIEW.md).
Chan–Pak completeness remains OPEN.
