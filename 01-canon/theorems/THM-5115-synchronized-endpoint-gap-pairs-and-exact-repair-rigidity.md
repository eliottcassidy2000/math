---
id: THM-5115
title: "Pure-count and synchronized physical gap pairs with controlled root sidecars and exact continuation restrictions"
status: "PROVED, independently audited. Universal gap existence inherited; synchronization and scoped repair statements only."
session: Session244-even-gap-anchoring-20261110
source: 05-knowledge/results/session244_even_gap_anchoring_20261110/composition/PROOF.md
audit: 05-knowledge/results/session244_even_gap_anchoring_20261110/audit/REPORT.md
---

# Prescribed gaps with a synchronized physical root state

At every rooted source(D,E,U,d), the inherited physical O/Z/J maps give

    JJJJZ -> (9D-3E+11U,3D-E+3U,9D-3E+11U,3D-E+3U),
    JOZJZ -> (9D-E+11U,3D-E+3U,9D-E+11U,3D-E+3U).

For a strong source these are strong same-order graphs. Their gap is2E,
and they have identical endpoint and deletion counts, with response equal
to total count. Unlike Session243's endpoint compiler, the lower coefficient
row has no fixed factor3; literal prime pairs193/199,239/241,313/317 verify
that this is a genuine arithmetic distinction, not a prime-infinitude claim.

Applying THM-5113 with E=N gives sourceg=g(N),epsilon=N mod2 and counts

    L=31g+37N-11epsilon,  R=31g+39N-11epsilon <=70N,

with shared endpoint/deletion9g+11N-3epsilon and order at mostceil(log2N)+8.
All-even-gap existence remains inherited from Sessions199/221.

On the saturated physical face(h,e,h,e), the two five-letter maps are

    Fminus(h,e)=(20h-3e,6h-e),
    Fplus(h,e)=(20h-e,6h-e).

For odd h,e, Fminus preserves gcd(h,e), and the gcd of the two branch total
counts is gcd(5h,e). Its inverse is h=(h'-3e')/2, e=3h'-10e'. These arithmetic
conditions never replace physical predecessor realization.

For a synchronized difference q(1,0,1,0), q>0, a common O/Z/J continuation
returns to that same synchronized difference direction if and only if it
is Z^k. Its gap becomes2^k q. The proof tracks the nonnegative physical
difference cone: a J makes the deletion difference positive permanently;
without J, any O causes a strict mismatch between response and total.

The natural asymmetric branch swap is also obstructed: apply Fplus to the
lower new compiler output and Fminus to the upper. Their difference is

    (6(7E-D-U),12E,6(7E-D-U),12E).

Every coordinate is divisible by6, and every subsequent common integer
transfer preserves that divisibility. This route can never repair to gap2.
These statements concern the specified words and synchronization, not all
possible graph operations or arbitrary asymmetric repair.

The complete zero-defect classification, H>=3e and e=1 dyadic slice used
here are INHERITED from Session231 descent/ROUND2.md; the degree bound is
inherited from THM-4873/4875. The packet retains independent proofs and
scoped corollaries, not a new claim to those theorems.

## Stronger pure-count pair

For every rooted source, two seven-letter words have states

    OJJJJOZ: (17D+40U,3D+6U,21D+50U,3D+6U),
    OJOZJOZ: (17D+42U,3D+6U,21D+50U,3D+6U).

Only the total count changes; all three other rooted coordinates agree.
Choosing the inherited canonical response source U=N gives gap2N, upper
count17g(N)+42N<=59N, and strong same-order graphs on at most
ceil(log2 N)+8 vertices. N=1 uses the singleton source; the opening O
handles that boundary. This is a stronger synchronization condition than
the five-letter pair above, not a replacement for its smaller word length.

Every common O/Z continuation preserves the exact gap and all three equal
sidecars. Thus the anchor may change along an arbitrary physical binary
tail without changing the prescribed difference. A first common J sends
the pure-count difference to a positive multiple of the singleton state
(1,1,1,1); later physical transfers cannot make its deletion coordinate
zero again. Hence a common word preserves this pure-count synchronization
exactly when it contains no J. This gives no arbitrary-anchor selector.

[Pure-count proof and literal checks](../../05-knowledge/results/session244_even_gap_anchoring_20261110/composition/ROUND3.md).

[Compiler and checks](../../05-knowledge/results/session244_even_gap_anchoring_20261110/composition/PROOF.md) ·
[Efficient substitution](../../05-knowledge/results/session244_even_gap_anchoring_20261110/composition/ROUND2.md) ·
[Exact repair restrictions](../../05-knowledge/results/session244_even_gap_anchoring_20261110/logic/ROUND2.md).
