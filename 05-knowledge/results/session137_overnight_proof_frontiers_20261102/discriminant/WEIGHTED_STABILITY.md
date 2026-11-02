# Polynomial stability of the weighted determinant bound

**PROVED, independently audited.** This is a quantitative continuation of the [global theorem](GLOBAL_PROOF.md), not a new proof of the full integer value-spectrum conjecture.

For a weighted tournament matrix A of order n let H be its weighted Hamiltonian-path sum, u=det(I+A), v=(-1)^n det A, and

    delta=H-u,       Delta=H-u-|v|,
    V=sum_(i<j) A_ij A_ji.

The unordered-edge variance V vanishes exactly for deterministic tournaments.

## 1. Cofactor lower bounds

For a nonempty r-vertex weighted remainder put K=A-A^T, M=I+K and s=1^T M^(-1)1. Every principal minor of a real skew matrix is nonnegative: odd-order minors vanish and even-order minors are Pfaffian squares. The principal-minor expansion therefore gives det(M)>=1 from its empty term.

Augment K by the all-one root to form the skew matrix

    L=[[K,1],[-1^T,0]].

Schur complementation gives det(I+L)=det(M)(1+s). In the principal-minor expansion, terms avoiding the new root sum to det(M); terms containing it include r two-by-two minors, each equal to1, and all other terms are nonnegative. Consequently det(M)s>=r. The determinant coordinates in the global proof give

    u+v >= 2^(1-r),       u-v >= r*2^(1-r),
    u >= (r+1)/2^r,
    u-|v| >= 2^(1-r).                                  (1)

The bound on u also holds at r=0, where u=1. The final signed lower bound must not be applied to the empty remainder, where u-|v|=0.

## 2. Averaging the weak gap over deletion roots

Fix a root x and set I_x=sum_(i!=x) A_xi A_ix. Retain the old one-vertex paths in the global weak-gap recurrence; every discarded summand is nonnegative. This yields, for n>=2,

    delta(A) >= delta(A-x)+t_n I_x,
    t_n=(n-1)/2^(n-2).                                 (2)

The base is the exact identity delta=V at order two. An edge survives in n-2 deletion cards, and sum_x I_x=2V. Averaging (2) and inducting gives

    delta >= c_n V,            n>=2,
    c_n=((n-2)c_(n-1)+2t_n)/n,
    c_n=[24-2^(3-n)(n^2+2n+3)]/[n(n-1)].                (3)

For the closed form put a_n=n(n-1)c_n. Then a_2=2 and

    a_n-a_(n-1)=2^(3-n)(n-1)^2.

Summing this elementary finite power series gives (3).

## 3. The strong gap has polynomial stability

Apply the two signed-margin identities, retaining the parent weak gap and the one-vertex paths. Their remainders have size n-2, so for n>=3 equation (1) gives

    Delta(A) >= delta(A-x)+2^(3-n) I_x.                 (4)

This holds for every x and both signs before taking the minimum. Averaging (4) and using (3) proves

    Delta >= beta_n V,           n>=3,
    beta_n=[24-2^(4-n)(n^2-n+3)]/[n(n-1)].               (5)

The numerator starts at6 at n=3 and its next increment is

    2^(3-n)(n^2-3n+3)>0.

Hence

    beta_n>=6/[n(n-1)],       beta_n~24/n^2.             (6)

The domain is essential. At order two Delta is identically zero, although V can be positive.

## 4. Exact small-order sharpness

At n=3 let z=det A be the sum of the two directed triangle weights. Direct expansion gives

    H=1+2z,       u=1-V+z,       v=-z,
    Delta=V.

Thus beta_3=1 is exact for every weighted order-three matrix.

At n=4 take the ordinal sum of two symmetric weighted two-vertex blocks, all internal arc weights1/2 and all cross-block arcs forward. Then

    H=1, u=9/16, v=1/16, V=1/2,
    delta=7/16, Delta=3/8.

Thus c_4=7/8 and beta_4=3/4 are sharp. This session makes no optimality assertion for beta_n at higher orders.

## 5. Rounding distance

Let d(A)=sum_(i<j) min(A_ij,1-A_ij), the distance on unordered edges to a nearest deterministic tournament. Since min(p,1-p)<=2p(1-p),

    d(A)<=2V<=2Delta/beta_n<=n(n-1)Delta/3.              (7)

For the full off-diagonal entrywise matrix L1 distance, multiply d by two. This bounds distance to a deterministic matrix; it does not by itself assert that a rounded matrix is transitive or identify its isomorphism class.

## Verification and provenance

The [independent proof audit](../independent/WEIGHTED_POLYNOMIAL_AUDIT.md) verifies the root recurrences, normalization, coefficient summation, n=2 hostile, and n=3/n=4 controls. Its exact Fraction/DP checker covers856 rational matrices,3,627 weak and3,621 strong deletion-root inequalities, all27 ternary order-three identities, and coefficient recurrences through n=200. Normal and optimized receipts match.

    cd ../independent
    python3 -B audit_weighted_polynomial.py
    python3 -B -O audit_weighted_polynomial.py

The augmented-skew/rooted-minor mechanism is consistent with inherited THM3729; the global positivity and signed-gap identities are those proved in the first checkpoint. This quantitative result does not provide an efficient exact H algorithm, a full-spectrum realization theorem, or a claim of external priority.
