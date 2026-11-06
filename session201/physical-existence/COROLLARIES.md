# Asymmetric budgets and a sharper J crossing

Session 201. **PROVED, independently audited in the stated scope.** This companion uses the physical cone, adjoint and exact letter matrices in PROOF.md. It does not change the frozen main proof. No universal coverage or globally optimal pair-height constant is claimed.

## 1. Exact asymmetric physical inventories

**Corollary A.** Let the target be n=F(s) for a full C1 source, or n=G(s) for a nontrivial-precursor balanced source. Write h0=2 in the full case and h0=3 in the balanced case. Let L,R be positive real numbers such that

    L>=2, R>=h0, LR>=3(n-1).

Then some cut of the source word gives a physical pair K(v,w)=n satisfying

    U_v<=L and U_w<=R.

For synchronized joint membership (m,h), substitute n=m and h0=3 and transport (b,z) together. The full auxiliary pair and total J cost are preserved. Since cover heights are integers, the corresponding inventory cutoffs are floor(L),floor(R). These floors need not themselves have product at least3(n-1): the assumptions are imposed on the real budgets before rounding. An integer-budget version uses integer L,R with the same displayed inequalities.

**Proof.** Along the actual word let i be the last cut with left height at most L. It exists because the seed has height2. If i is terminal, its right boundary has height h0<=R and the cut works. Otherwise the next left height is greater than L. If the right height at i were also greater than R, their product would be greater than LR>=3(n-1). Those two heights are exactly U(c v) and U(sigma(c) w) in the adjacent-letter inequality, which bounds their product by3(n-d_v d_w)<=3(n-1). This contradiction proves the right height at i is at most R. The word is only repartitioned. QED.

The separate balanced value3 remains outside the triangle-precursor claim and is handled directly, as in the main inverse. This corollary allows unequal inventory limits when one side is cheaper to enumerate. It does not assert the existence of a pair for a target not already in the spectrum.

## 2. A letter-specific improvement at J

Let

    C=P^T J^T Q P
      =[[2,1,4,3],[1,1,2,2],[4,2,7,6],[3,2,6,5]],
    A=(2,1,4,3).

Then

    16C-7A^T A
      =[[4,2,8,6],[2,9,4,11],[8,4,0,12],[6,11,12,17]]

is entrywise nonnegative. Therefore, on the full positive cone and in particular for actual physical sources,

    U(Jv) U(Jw) <= (16/7) K(Jv,w).              (1)

**Corollary B.** If the left/right height crossing of a specified source word occurs at J, the better of the adjacent cuts has maximum height at most

    floor(sqrt(16n/7)).

This follows from min(A,B)^2<=AB and (1), with the same cut orientation as the main proof. Integer implementation can use isqrt((16*n)//7). The assertion is conditional on the crossed letter being J; there is no claim that every target admits such a crossing.

The factor16/7 is sharp for this J-only inequality and for cuts of specified words. For k>=1, let d=2^k+1 and

    v=w=t.Z^k.J=(2d-1,d-1,2d,d).

Then

    n=K(Jv,v)=D_v^2-E_v^2+U_v^2=7d^2-2d,
    U(Jv)=4d-1,
    (4d-1)^2/(7d^2-2d) -> 16/7.

The literal source word Z^k J J J O^k has this F count. Its two middle cuts, immediately before and after the central J, both have maximum height4d-1. All earlier cuts have at least this right height and all later cuts at least this left height, so the best maximum across every cut is exactly4d-1. The k>=1 restriction is essential: the empty Z block leaves t, not the balanced formula for t.Z^k.

Again these sharpness statements do not rule out a better representation of the same numerical target.

## 3. Reproduction and scope

    python3 verify_corollaries.py > COROLLARIES_NORMAL.json
    python3 -O verify_corollaries.py > COROLLARIES_OPTIMIZED.json
    cmp COROLLARIES_NORMAL.json COROLLARIES_OPTIMIZED.json

The checker uses exact integer and rational arithmetic, checks the matrix by raw cone rays, checks integer and noninteger budget cases on a declared complete word universe, and verifies the J sharp family and a literal graph. All-order conclusions use the proofs above, not finite extrapolation.
