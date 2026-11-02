# Independent audit: polynomial weighted stability

**Verdict: ACCEPT, all orders in the stated domains.** This continues the
separate phase-two stability note and does not alter the first-checkpoint
global determinant proof.

For a weighted tournament A put `V=sum_{i<j} Aij Aji`, `u=det(I+A)`,
`v=(-1)^n det A`, `delta=H-u`, and `Delta=H-u-|v|`. Then

    delta >= c_n V,        n>=2,
    c_n = [24-2^(3-n)(n^2+2n+3)]/[n(n-1)],

    Delta >= beta_n V,     n>=3,
    beta_n = [24-2^(4-n)(n^2-n+3)]/[n(n-1)]
           >= 6/[n(n-1)].

Both domains are essential. The delta base case is the exact identity
delta=V at n=2. At n=2 the strong gap Delta is zero even when V is positive.

## Root recurrences and averaging

The augmented-skew argument in `WEIGHTED_STABILITY_AUDIT.md` proves, for a
nonempty r-vertex remainder, `u+v>=2^(1-r)` and
`u-v>=r 2^(1-r)`. Hence `u>=(r+1)/2^r`. The latter formula also holds for
the empty remainder with r=0 and u=1.

In the weak determinant-gap recurrence, retain only the old one-vertex
paths at a fixed inserted root x. With
`I_x=sum_{i!=x} Axi Aix` this gives

    delta(A) >= delta(A-x) + t_n I_x,
    t_n=(n-1)/2^(n-2),                    n>=2.

Induct from delta_2=V. Average over roots, using that an edge survives in
n-2 deletion cards and `sum_x I_x=2V`:

    c_n = ((n-2)c_(n-1)+2t_n)/n.

Putting `a_n=n(n-1)c_n` gives

    a_2=2,
    a_n-a_(n-1)=2^(3-n)(n-1)^2,
    a_n=24-2^(3-n)(n^2+2n+3).

For the two signed strong margins, retain the parent weak gap and the same
one-vertex paths. Taking their minimum gives, when n>=3,

    Delta(A) >= delta(A-x) + 2^(3-n) I_x.

The remainder is now nonempty, as required by the uniform signed-margin
bound. Averaging and using the already-proved weak coefficient gives

    beta_n=((n-2)c_(n-1)+2^(4-n))/n,

which simplifies to the displayed formula. Its numerator B_n starts at
B_3=6 and satisfies

    B_(n+1)-B_n=2^(3-n)(n^2-3n+3)>0.

Thus the simpler 6/[n(n-1)] bound is valid. The coefficient beta_n is
asymptotic to 24/n^2; no optimality claim at general n is justified here.

An immediate metric consequence is also fully scoped: let
`d=sum_{i<j} min(Aij,1-Aij)`, the entrywise distance on unordered edges to
the nearest deterministic tournament. Since `p(1-p)>=min(p,1-p)/2`,

    d <= 2V <= 2 Delta/beta_n <= n(n-1)Delta/3.

For the full off-diagonal matrix entrywise L1 distance, multiply d by two.
This is closeness to a deterministic tournament; it does not by itself
identify that rounded tournament's isomorphism class or transitivity.

## Exact low-order boundaries

At n=3 put z=det A, the sum of the two nonnegative directed triangle
weights. Direct determinant expansion and the weighted odd-cycle formula
give `u=1-V+z`, `v=-z`, `H=1+2z`. Therefore **Delta=V exactly**.

At n=4, take the ordinal sum of two weighted two-vertex blocks, each with
internal probability 1/2. Then

    H=1, u=9/16, v=1/16, V=1/2,
    delta=7/16=(7/8)V,
    Delta=3/8=(3/4)V.

Consequently c_4=7/8 and beta_4=3/4 are sharp. In particular beta_n is
known sharp at both n=3 and n=4, not just at three. This isolated control
does not imply sharpness at larger orders.

## Independent exact replay

`audit_weighted_polynomial.py` imports only independent local kernels. It
uses weighted Hamilton subset DP and exact Fraction Gaussian determinants
for every full matrix and every deletion card. It checks all ternary-edge
matrices through order four and 24 fixed-seed rational matrices at each
order five through eight: 856 matrices, 3,627 weak root inequalities, 3,621
strong root inequalities, and all 27 ternary order-three identities. Closed
coefficient recurrences are additionally checked for n=3..200. The n=4
sharp example and the n=2 false-extension hostile are explicit.

    python3 -B audit_weighted_polynomial.py > weighted_polynomial.normal.json
    python3 -B -O audit_weighted_polynomial.py > weighted_polynomial.optimized.json
    cmp weighted_polynomial.normal.json weighted_polynomial.optimized.json

Both receipts have SHA256
`5a953fb038f72f19b256e0962706907efdc48f15f4003b243b267c896293fa7b`.
