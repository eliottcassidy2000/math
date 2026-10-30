# Session 117: provisional coefficient-four audit

## Sources and status

Read-only source inspected on 2026-10-30:

- Rules: https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/2cd3a17f3f296b7f0ac248a63a9b54a0fabdf851/rules.md
- Candidate elimination proof: https://github.com/Layr-Labs/circuitbreaker-llm-verifier/blob/4a61f024a4993592d97b49ca7d00979029b7fece/submission/elimination.md
- Public submission: ArthurOstapenko, `5e9fcf7f-8338-4729-9faa-a18994661af6`, visible in https://www.yukon.org/circuitbreaker . Public metadata gives automated validity, human=null, status=review.

The exact candidate proof is available at the pinned source link above; its local read-only capture is not republished. No submission code was executed. The candidate is NOT an established lower-bound record merely because a leaderboard displays 4. The public note explicitly says this is a new claim not peer reviewed.

The new class claim is: if f is nonconstant on every affine subspace of dimension >d, then C_B2(f)>=4n-9(d+1). The proposed explicit-family assembly claims 4N-O(N^(4/5)) with actual prefix inputs counted. The elimination lemma, not a basis conversion or padding trick, is the essential new claim.

## Targeted inspection verdict

No concrete counterexample or blocking gap established in this pass. This is NOT a complete proof certification. In particular, the separate `disperser.md` and full family assembly were not independently verified here. Source claims and our independently established conclusions must remain separate.

Checks made against the exact elimination text:

1. The affine role exchange substitutes only an UNPAIRED free coordinate. For each assignment of other coordinates, x -> x+ell is bijective. Thus the product of existing three-point pair domains is preserved, even when ell mentions paired coordinates. The matrix-column replacement determinant argument is valid over F_2.
2. Safe fixing of one member of a three-point pair leaves its partner free; m and q both fall by one. Hence a=m-q is preserved. The safe affine core has dimension a, not m.
3. GP removes the gate and its actual children once, not independently priced gadgets. It explicitly handles unrestricted sharing and repeated pins. It restricts the domain rather than claiming a nonlinear affine coordinate change.
4. Terminal-output exclusion uses an affine safe core of dimension a and a constant slice of dimension at least a-1. The loop guard a>=d+2 is essential.
5. AP2's shared partner-child case uses the binary-pin limit to ensure a child of G differs from the shared partner child. AP1's additional third deletion requires its explicit output contradiction. These are the right points for a subsequent mechanized invariant audit, not immediate refutations.

There is no justified smaller circuit construction for the candidate family in this report. Finite SAT on small n cannot settle a bound whose negative error dominates at those n.

## Independently proved result: exact ceiling for the stated potential template

This result is elementary algebra about the source's row inequalities, independent of whether its universal gate-elimination theorem is ultimately correct.

Let W=g*s+i*I+p*q, with g>0 and i,p>=0. Demand that every listed zero-dimension-loss row has nonpositive potential change, and every dimension-loss-one row loses at least t. The row bounds from Section 7 require

    i <= 2g+p                    [P2 / PE]
    t <= 2g+2i-p                 [GP, one child]
    t <= 3g+2p                   [AP2]
    t <= 3g+i+p                  [AP1].

The telescope can certify at most coefficient c=(t-i)/g, with dimension penalty b=t/g. The first two inequalities imply

    c <= 2+(i-p)/g <= 4.

Therefore no retuning of these nonnegative weights, while relying on precisely these worst-case row estimates, can certify c>4. This is a limitation of this proof template, NOT an upper bound on Boolean circuit complexity or all possible potentials.

There is also a sharp tradeoff. From AP2, p>=(t-3g)/2. Combining with GP gives i>=(3t-7g)/4, hence

    c <= (b+7)/4,   and so b >= 4c-7.

At c=4 one must have b>=9. Equality is attained by (g,i,p,t)=(1,5,3,9). Consequently the source weights are simultaneously optimal for leading coefficient and, among coefficient-four weights, for the d+1 penalty within this template.

For 5/2<=c<=4, the tradeoff is attained by

    g=1, i=3c-7, p=2c-5, t=4c-7.

All four displayed inequalities hold: GP and AP2 are equalities; P2 is equivalent to c<=4; AP1 is slack for c>=5/2. In particular c=7/2 gives (i,p,t)=(7/2,2,7), and c=4 gives (5,3,9). These are conditional row-certificate tradeoffs, not independently established circuit lower bounds.

## Concrete next target

To exceed four, improve a structurally binding case, enrich the potential, or prove that worst-case P2/PE and one-child GP behavior cannot accumulate at full density. A useful finite search should enumerate legal restricted-state transitions, record deletions with shared descendants, and seek a GLOBAL amortized incompatibility between those two bottlenecks. Merely solving the same linear program again cannot improve the coefficient. Any proposed improvement must retain the affine safe-core dimension and exact product-domain/affine-system semantics.
