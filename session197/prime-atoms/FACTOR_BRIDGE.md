# A finite factorwise inverse bridge for all unit binary amplifiers

Session 197. **Exact all-order physical reduction, including composite amplitudes.** This is not a proof of the full tournament spectrum or of all prime/7p coverage. It uses the bounded-slack theorem in `PROOF.md` and the complete star synthesis from S196.

## 1. Keep the two pieces separate

Every positive binary transitive-cloning amplifier, possibly after taking a converse, arises from two rooted pieces X,Y sharing only their root c and a transitive chain Z. Write A=X-c, B=Y-c and orient A->Z->B, A->B, Z->c. The S196 theorem gives

    H(F_m(X,Y))=alpha*2^m+beta,
    alpha=U L,
    beta=d g-f t,                                   (1)

where

    d=H(X-c), f=R_c(X)+ends_c(X)-H(X)-d,
    U=R_c(X)-d>0,
    g=H(Y), t=H(Y-c), L=L_c(Y).

All of d,g,t are positive odd integers; t<=g<=L. We retain actual rooted pieces throughout. A numerical solution of (1) without physical piece realizations is not a certificate.

## 2. A stronger source dichotomy than coprimality

The first packet's endpoint reinsertion injection gives, for every source with U>0,

    d+2f <= H(X)+2f <= 2U-1.                        (2)

Suppose beta is +1 or -1.

- If beta=-1, then ft=dg+1>=dt+1, so f>d. Since f is even and d odd, f>=d+1. Equation (2) gives 3d<=2U-3.
- If beta=+1 and f<d, write

      1=(d-f)t+d(g-t).

  Both terms are nonnegative integers, with d-f positive. Hence t=1, d-f=1 and g=t=1. In particular the partner is transitive. Its response L is a power of two.
- In every other beta=+1 case, parity gives f>d and again 3d<=2U-3.

Thus every unit synthesis satisfies the uniform inverse budget

    3d <= 2U+2beta-1,  f<=U-1.                      (3)

A positive source offset d-f can only be +1 and can only couple to a transitive partner. All nontrivial cancellation with a varying partner begins from a negative source offset. This is a restriction on physical inputs, not merely the determinant's coprimality condition.

## 3. Explicit bounded-order certificates for a fixed factorization

Let k(s)=floor((1+sqrt(1+8s))/2), ell(s)=floor(log2(max(1,s))), and v2(n) be the exponent of two dividing a positive integer n. For U>=1 put

    D_beta(U)=floor((2U+2beta-1)/3),
    B_X(U,beta)=1+D_beta(U)+k(U-1)+ell(U-1)+v2(U).    (4)

A factorization with D_beta(U)<1 is impossible. For L put

    B_Y(1)=1,
    B_Y(L)=L+k(L-1)+2ell(L-1), for L>=2.             (5)

**Theorem 5 (finite factorwise bridge).** Fix positive alpha and beta in {+1,-1}. There exists any physical positive binary amplifier with ray alpha*2^m+beta if and only if there is a factorization alpha=U L and actual rooted pieces X,Y satisfying (1), with

    |X|<=B_X(U,beta), |Y|<=B_Y(L).                   (6)

Consequently an existing ray always has a certificate graph F_1 of order at most B_X+B_Y=O(U+L+log alpha), with the original alpha and beta preserved. A balanced factorization, when it physically exists, gives a certificate of order O(sqrt(alpha)); no balanced factorization or physical realization is asserted for arbitrary targets.

**Proof.** The backward implication is the literal S196 graph compiler.

For the forward implication, start with any physical star presentation. Normalize X using the first packet's finite-core theorem. It becomes profile-equivalent to K^j(C,c), with

    |C|<=1+d+k(f)+ell(f), U=2^j U_C.

Since U>0, j<=v2(U). Equations (3) and f<=U-1 prove the source bound (4). Rebuilding the K-tail gives a bounded actual X with the same d,f,U, hence does not change (1).

For Y, take the converse W, rooted at c. Its starting-cover response R_c(W) equals L_c(Y)=L. Its profile has h=g, deletion count t and endpoint count s=starts_c(Y), where s<=t. Thus its slack F obeys

    F=L+s-g-t <= L-g <= L-1,
    U_W=L-t.

If U_W=0, then t<=g<=L=t gives g=t=L. Equation (1) becomes L(d-f)=beta, so L=1. Such a partner is transitive with its root at the source and can be replaced by a singleton.

For L>=2 we therefore have U_W>0 and t<=L-1. Normalize W to a core plus a K-tail of length j. The tail satisfies j<=v2(L-t)<=ell(L-1); the core order is at most 1+t+k(F)+ell(F). Hence its rebuilt order is at most (5). Taking the converse returns a bounded actual Y with exactly the same g,t,L. The star identity preserves alpha and beta. QED.

This yields a complete finite inverse procedure for a specified ray: enumerate the factor pairs U L=alpha, then all rooted graphs within the explicit bounds, and test their actual profiles against (1). Successful output is a literal graph; exhaustion proves absence only within the globally classified unit binary-cloning architecture. The bounds are not intended as practical brute-force estimates.

For a specified count N with a unit cloning root, the amplitude is necessarily (N-beta)/2. Thus two finite factorwise tests, one for each sign, decide the existence of such a root at every graph order. This is stronger than searching a fixed finite template bank. It still does not decide whether N is a tournament Hamilton count by a different construction.

## 4. A complete two-sign obstruction at H=15

The first packet already excludes amplitude7, offset+1. This removes a returning H->2H-1 clone root at H15. The remaining unit possibility at H15 has alpha=8, beta=-1.

In the latter case U divides8. If U=1 or2, the budget 3d<=2U-3 is impossible. If U=4, it forces d=1. If U=8, it forces d=1 or3.

For d=1, the S196 binary incidence recursion gives f=0 whenever U is a power of two. But a negative unit determinant requires f>d, so this case is impossible.

For d=3 and U=8, f>d is even, and (2) bounds f<=6. Thus f=4 or6. The f=4 row is already in the first packet: its U bases are {9,15}. Applying the same exact finite-core enumeration at d=3,f=6 gives normalized bases

    {9,15,27,36,60},

again times arbitrary powers of two. None can equal8. This small new table has the same all-order justification: the deletion core is a triangle, |O|<=4, and at most two neutral singleton components can lie between the first and last out-neighbor components. The verifier enumerates every reduced core under those proved bounds, checks its profile by full subset-DP, and saves witnesses.

**Corollary 5.1.** No tournament with H=15 has any positive unit-offset binary transitive-cloning root, of either sign, at any graph order.

H15 itself is attained, for example by the strong tournament with rows [10,12,9,16,7]. The corollary blocks a universal binary-cloning normalization, not this Hamilton count or another operation. It does not exclude a one-off equality H(T[r:2])=2H(T)±1 at a root with higher clone-run degree, because that does not imply an indefinitely returning binary ray.

## 5. Verification and scope

`verify_factor_bridge.py` checks the general source bounds and the constructive piece normalizations on all unit pairs drawn from complete rooted labelled banks through order4, then directly checks the star graph counts. Larger padded pairs exercise removal and relocation of neutral components. The extra d=3,f=6 classification uses the proved finite-core bound and is independently checked on its literal witnesses. Normal and optimized runs must match.

The bridge isolates a concrete remaining problem: for each required atom N, produce actual factor pieces satisfying (1) within these bounds, or use a different operation. The theorem supplies a sharp finite physical certificate language and an explicit failure at15. It does not turn an integer factorization or Bézout solution into a tournament, and supplies no unsupported all-target closure.
