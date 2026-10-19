# Independent review: original-height repair for the K/V–F/H bank

19 October 2026. **Verdict: PASS for the stated bank-specific theorem.** No counterexample or substantive proof gap was found. The conclusion is a descent-or-smaller-periodic-obstruction dichotomy, not a global Collatz proof and not an endpoint zero-head theorem.

## Scope audited

The source x is a positive integer following a completely guarded Terras parity word. Its macroblock factorization is a nonempty word over {K,V=K^4F}, followed by a nonempty word over {F,H}, and its full affine multiplier is strictly less than one. Arbitrary ordering is allowed within each segment; arbitrary interleaving between segments is not covered. N is the full number of Terras steps.

Under precisely these hypotheses, either some iterate at a time 1 through N+92 is below x, or the lower word is a cycle of least period N and minimum x, with a distinct nontrivial positive cycle containing (x+1)/4<x.

## Proof audit

1. **Displacement envelope.** The two composition inductions give E(x)=Ax+B with B≤C(A−1), and R(x)=rx+d with d≤1−r. If s=rA<1, the normalized carry e satisfies e≤Cs+1−(C+1)r<C+1<534. Consequently a nonnegative integer endpoint displacement is one of 0,…,533. This does not require a lower bound on 1−s.

2. **Forced prefix.** An expansion word either is K, K², K³, or its full parity word starts with K⁴. Each short case descends at its first reset, for every lawful source x≥815464443. Both reset maps are nonincreasing on positive integers. Thus every non-descending full endpoint forces the K⁴ cylinder. Independent binary lifting recovered exactly

   r4=588943162011782273059357097228697595 modulo 2^120.

3. **Positive displacement.** Independently computed all 533 offset certificates using the closed sum

   b = sum over odd positions i of 2^i 3^(number of later odd positions),

   rather than the writer's coefficient recurrence. Each word length j is at most 92, its modulus divides 2^120, and the exact threshold (a d+b)/(2^j−a) is below 9000. The anchor's word therefore remains lawful for every inherited source lift, and x≥r4>9000 yields actual descent below x. This is a finite symbolic cover of offsets with an all-lift conclusion; replay at sample lifts is not the proof.

   Maximum depth: 92, uniquely at d=4.
   Maximum threshold: 116927/13, attained at d=480.
   The 533 checks are independent of bank-word length.

4. **Zero displacement.** For every nonempty K/V prefix its slope exceeds one. For every proper F/H prefix the difference from the block endpoint has positive slope and a nonnegative value at input 1; hence it is strictly positive for every input greater than 1. In a return to x≥r4, reset endpoints strictly decrease to x. It follows that T^t(x)>x for all 1≤t<N, proving both the minimum assertion and least period N. This checks intermediate Terras steps, not merely macroblock endpoints.

5. **Companion lawfulness.** Independent binary-lift reconstruction verified each upper guard is the corresponding lower guard plus one. Closed coefficient calculations verified equal length, equal multiplier numerator, and upper carry b+m−a for every bank block. The adjacent relation therefore propagates through actual block concatenation, giving T^N(x+1)=x+1. Its first two bits are 00, so (x+1)/4 is a positive periodic integer below x. It exceeds 2 and lies in a cycle distinct from the lower cycle by the latter's minimum. No claim of a descending lower-orbit iterate follows from the companion alone.

6. **Least-counterexample corollary.** In the descent case, minimality forces convergence through the smaller lower iterate. In the periodic case, the smaller companion periodic point contradicts minimality. This only excludes a least nonconverging integer from starting such a contracting bank word; it supplies no bank-entry theorem and does not transfer the original-height guarantee to a later entry after growth.

## Boundary and wording

The proof should use “nonempty proper lower prefix,” or explicitly “1≤t<N,” when asserting strict height above x. The empty prefix equals x. This is a terminology correction, not a failure of the argument.

The standard F pair 10/01 at sources 1 and 2 shows why affine adjacency alone cannot establish distinct cycles: both sources lie in the same trivial cycle. Here the large inherited source, independently lawful 00 start, and strict lower-cycle minimum are essential.

## Reproduction

The independent checker check_escape.py and its receipt.json accompany this review. Its normal and optimized executions agree. The writer's verifier was also replayed normally and with optimization: both agree byte for byte with verification-receipt.json and report PASS, 5,539 explicit checks, 54 compressed classes, and repair depth 92.

Only the new theorem, its short bank interfaces, and its bounded offset certificates were audited. No old large bank census was rerun. The literature paragraph was not part of this mathematical audit.
