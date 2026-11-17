# Decimal scales and a recurring nonphysical Pell branch

This sharpens the [least-hole and near-twin checkpoint](README.md). The constructions do not certify arbitrary offsets or prime quotients.

## The integer inverse is single-valued on each candidate

Use its notation (K=uv), (W=4vR-K+4), and (uS-vR+2uv=2). The two-slot counts are even, so the determinant-one condition implies gcd(u,R)=1. Thus reducing the fraction `(W+K−4)/(4K)` to lowest terms recovers R/u uniquely. Then v=K/u and S=((R−2u)v+2)/u are forced. No divisor enumeration of K is needed. Positivity, parity, order constraints and physical realization are still required.

Except at the pair 29→31 both cores have order at least two. If X is a singleton then u=1,R=2 and the switch forces S=2, hence E is also a singleton. Conversely if E is a singleton, v=1,S=2 forces R=4u−2, impossible for nonsingleton X. Indeed P₂≥(n+1)H by splitting each Hamiltonian path at every cut (including the two empty slots), and P₂=4H at order two. Thus R≥4u and S≥4v when neither is a singleton. This improves the complete conic search bound to N≥89K in the nonsingleton cases, by N=2(2R+u)(S+v)−K.

## A repeating branch that cannot be a tournament profile

For any odd N>32 with N≡32 mod 39, the relaxed square equation has the exact solution

    K=(N−32)/39,   W=13K+4.

Direct substitution gives equality. Here K is odd and (W+K−4)/(4K)=7/2. Since u=H(X) must be odd, no profile of a tournament can have R/u=7/2. It also violates W+K≡0 mod 4, and N≥89K for sufficiently large N. This branch recurs at N=10^n+7 whenever n≡3 mod 6, since 10³≡25 mod 39 and 25²≡1 mod 39. It is a genuinely repeated arithmetic pattern, but it is a spurious solution of the relaxed inverse, neither a path-count witness nor a proof that another branch cannot work. At n=9 the previous conic exclusion and the two explicit strong witnesses remain valid.

## Real strong tournaments at every decimal height

Use these literal tournament masks with each row storing its outgoing vertex bits:

    X = (24,1,3,6,14),                        (H,P₂)=(11,100)
    Y = (16,65,779,531,806,335,413,575,651,99), (H,P₂)=(4453,100000).

For every n≥4 write n=2a+5b with b=0 for even n and b=1 for odd n. Take the ordinal product B of a copies of X and b copies of Y. H and P₂ both multiply: in each ordered slot, paths in consecutive ordinal components concatenate uniquely. The cyclic quotient C₃[TT₁,B,TT₁] is strong and has count P₂(B)+H(B) by the two-cover gluing law. Its exact count is

    10^n + 11^a · 4453^b.

This is one literal recurrence on each parity, with increment n→n+2 appending the same block and multiplying the scale by 100 and the error by 11. The error exponent is log₁₀(11)/2, about 0.5207, on both parities. When a is odd the count is divisible by 3. So the construction supplies all sufficiently large decimal scales, not the constant offsets +7 and +9.

Replay the two input profiles by D[S,j]=Σ D[S without j,i] over i→j with singleton initial values 1; set h(S)=Σ D[S,j], h(empty)=1 and sum h(S)h(complement S) over all subsets. This gives the displayed pair of integers for each input.
