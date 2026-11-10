---
id: THM-5118
title: "Synchronized source-block quadratic lifts, palindrome defects and their complete physical target inverse"
status: "PROVED, independently audited. Exact family construction and inverse; Mersenne coverage inherited."
session: Session245, 2026-11-10
source: 05-knowledge/results/session245_physical_spectrum_20261110/descent/PROOF.md
audit: 05-knowledge/results/session245_physical_spectrum_20261110/audit/AUDIT.md
---

# Synchronize two physical scales without assuming a free anchor

Use Session199's actual rooted O/Z/J grammar from t=(1,0,2,1). For a fixed
middle word W, define formal coefficient vectors

    b=W(1,0,0,1), a=W(0,0,2,0).

Neither vector is asserted physical. For s>=1 and k>=0 the actual graph
word O^k W O^(k+s) J is strong with exactly one Hamilton cycle and has

    H = 2^s a_U (2^k)^2 + (a_D+2^s b_U)2^k + b_D.

Its order is 2k+|W|+s+3. This follows by applying W to the actual state
(1,0,2*2^k,1), then resetting E with the nonempty last O block. Thus exact
vanishing of the linear coefficient is a finite certificate for an infinite
physical ray. It does not freely prescribe the other two coefficients.

For a binary middle of length l, with Z-bits epsilon_i, put

    A=sum epsilon_i 2^i, R=sum epsilon_i 2^(l-1-i),
    C=sum_(j<i) epsilon_j epsilon_i 2^(i-j-1)
       +sum_(i>=1) epsilon_(i-1) epsilon_i.

For s=1 the count is

    H=2^(2k+l+2)+2(A-R)2^k+1-C.

The linear term vanishes iff W is a palindrome. In that case the exact
strong C1 family is 2^(2k+l+2)+1-C(W), on 2k+l+4 vertices. The empty word
is included. C is always even; the adjacent-pair correction is essential.

For any half-word V, let R be its reversed binary integer, C=C(V), and e
its final bit (zero when empty). Exact mirror identities are

    C(V reverse(V))=2C+R^2+e,
    C(V O reverse(V))=2C+2R^2,
    C(V Z reverse(V))=2C+2R^2+2R+2e.

The cross-half terms factor; central adjacency contributes the remaining e
terms. For length l>=1, 0<=C<=2^l-2. Hence a represented h>=9 has a unique
possible ambient power satisfying

    h-1 <= 2^m <= 4(h-3)/3.

Enumerate palindrome lengths l<=m-2 of parity m and test C=2^m+1-h;
acceptance supplies k=(m-l-2)/2 and an actual graph. Failure excludes h
from this entire family at all orders. There are O(sqrt(h)) half-word
candidates, each with polynomial-in-log(h) exact arithmetic cost. The sole
output below 9 is 5, handled separately. This is numerical-target, not
input-bit-polynomial, complexity.

The Mersenne examples ZZ/ZOZ are exactly the inherited THM-4748 Fibonacci
fans, not new all-Mersenne coverage. Genuine middle-root-change examples
include ZZOJOOZ, ZZJOOOZ, ZZJJOOZ, with all-k counts respectively
340*4^k-35, 436*4^k-37, 420*4^k-43. The complete four-block coverage lemma
and arbitrary-target physical descent remain OPEN.

[Full proof, compiler and exact receipts](../../05-knowledge/results/session245_physical_spectrum_20261110/descent/PROOF.md).
