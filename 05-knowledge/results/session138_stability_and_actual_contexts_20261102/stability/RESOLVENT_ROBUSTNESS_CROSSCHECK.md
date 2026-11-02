# Independent resolvent cross-check for variance robustness

Status: PROVED, independently audited by the Session138 audit lane. This note is
outside both frozen publication manifests. It corroborates the root
Pfaffian robustness theorem by a different mechanism, with weaker constants.
The primary sharp-rate theorem is `../asymptotic_stability/PROOF.md`.

## Dimension-free determinant margin

Let B be a weighted tournament of nonempty order r. Write K=B-B^T,
M=I+K, D=det(M), and s=1^T M^(-1)1. The inherited identities give

    u+v=D/2^(r-1),    u-v=Ds/2^(r-1).

All quantities on the right are positive. Since K is real skew-symmetric,

    ||M^(-1)||_op <= 1,
    s=||(M^(-1))1||_2^2=||(M^(-T))1||_2^2.

For the first inequality use ||(I+K)x||^2=||x||^2+||Kx||^2. The identities
for s follow from sym(M^(-1))=(I-K^2)^(-1) and the commutation of I+K with
I-K.

Vary one skew entry x=K_ij with derivative matrix J=E_ij-E_ji. Its operator
norm is1. Jacobi's determinant derivative and the inverse derivative give

    |d logD/dx|=|tr(M^(-1)J)|<=2,
    s'=-1^T M^(-1) J M^(-1)1,
    |s'|<=||(M^(-T))1|| ||J|| ||M^(-1)1||=s.

Consequently |d log(Ds)/dx|<=3, and both signed margins have logarithmic
Lipschitz constant at most3 per old skew coordinate.

Choose a nearest binary tournament B0, rounding each edge to its preferred
orientation (ties arbitrary). At B0 both u+v and u-v are positive integers,
hence at least1. The total skew-coordinate path length from B0 to B is
2d, where d=sum_e min(p_e,1-p_e). Integrating the logarithmic bounds along
coordinate segments gives

    u-|v| >= exp(-6d) >= exp(-12V),                  (1)

because p(1-p)>=min(p,1-p)/2. No dimension-dependent norm estimate or
assumption of transitivity of B0 is used. The empty tournament is excluded.

## Integrated gap check

In a full independent edge-rounding order at n>=3, let q_t be the current
edge variance and S_t the sum of the later unrounded edge variances. Its
pair-deletion remainder has variance at most S_t, so (1) and the one-edge
gap interpolation give a correction at least q_t exp(-12S_t).
The intervals [S_t,S_t+q_t] partition [0,V], and the decreasing exponential
has integral at most q_t exp(-12S_t) on each. Therefore

    Delta(A) >= E Delta(T)+(1-exp(-12V))/12.          (2)

This is weaker than the independently proved Pfaffian bound with4 in
place of12, but is sufficient for the same matching-extraction asymptotic
rate argument: the right side is bounded away from zero when V>=1.
Thus the root's O(n^(-2)) asymptotic error rate does not depend on a
single robustness mechanism.

## Boundaries

The determinant statement requires r>=1. The gap statement requires n>=3,
so every pair-deletion remainder is nonempty. At n=2 the strong gap is
identically zero even with V>0, and (2) would be false. No claim of optimal
constants in (1) or (2) is made.
