# Combined baseline and light-mass payment: exact stability for n>=314

**Historical checkpoint: valid but SUPERSEDED as the current frontier by [n>=48](STABILITY_48.md).** Its remaining-range sentences record its publication stage.

**Status: PROVED; same-agent proof audit and exact rational certificates, with separate local-energy checks. No independent-agent or formal-proof claim.** This combines the incoming
[n>=466 proof](../stability_priority_spectrum_20261102/STABILITY_466.md)
with the [mass-paid tail](STABILITY_514.md) and Session139's audited dual.

For a real weighted tournament A define V=sum_e p_e(1-p_e),
Delta=H-det(I+A)-abs(det A), R=floor(n/2), and

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every integer n>=314, Delta>=kappa_n V. This is sharp;
positive-variance equality consists exactly of maximal ordinal fair pairs,
with one singleton at odd order. At V=0 equality is exactly transitivity.
Together with inherited n=3..9, this leaves10..313 OPEN.

## 1. Inputs, common baseline and constants

Use the accepted Session138 rounding/matching inequalities and greedy
group invoices stated in [STABILITY_514.md](STABILITY_514.md). In particular
B=E Delta(T) always uses the same independent binary rounding, W is the
incident pair sum, and greedy group i has selected weight w_i and assigned
discarded mass S_i. Let D=sum S_i. Contributions from different groups
are disjoint and satisfy

    B>=W>=D^2/(2n),
    W_i>=w_i S_i,          W_i>=w_i S_i/2+S_i^2/4.

Assume V>0 and Delta<=kappa V, n>=314, R>=157, kappa=kappa_n. Put

    a=3/4, c=249999/1000000, d=32/157,
    t=d/8=4/157, lambda=d/64=1/314, h=t-4lambda=2/157.

Since kappa<=4/R<=4/157 and kappa<=8/(n-1), the same exp(4)>50
bootstrap gives V>49(n-1)/1600>4. Consequently

    kappa V>B+c,     B<=kappa D,
    V-D>c/kappa-n kappa/2,
    V<=n/8+4n/(n-1).                                  (1)

The bound V-D retains B and completes its square before subtracting D.
No extra copy of B is spent in the following argument.

## 2. Count retained edges with the missing slots kept

Suppose L>=1 selected greedy edges have w_i<=t, and let U be their total
weight. Write r_h for the number of retained heavy edges and
M_h=V-D-U for their mass. The slot constraint is

    L+r_h<=R,        U<=tL<=t(R-r_h).

Thus r_h>=4M_h and(1) give

    r_h(1-4t)>(c-4t)R-16-8/R=:g(R).                    (2)

Here kappa R<=4 and 2n kappa<=16+8/R. The coefficient1-4t=141/157
is positive. At R=157,

    g(157)-8(1-4t)=2225351/157000000>0.

Since g'(R)=c-4t+8/R^2>0, (2) implies r_h>8 and therefore r_h>=9.
In particular the heavy matching is nonempty. Set m=R-L, so

    9<=r_h<=m<=R-1,       n>2m.

This sharper slot count is needed: using only U<=tR would lose one
retained edge at the critical endpoint.

## 3. The retained matching keeps its negative baseline

Apply the matching coefficient a_(r_h), use its monotonicity, and KEEP B:

    (kappa V-B)m >=4(1-a^m)(V-D-U).

Writing m=R-L and expanding yields

    L[kappa V-B]-4(1-a^m)U
      <=[kappa R-4+4a^m]V+4(1-a^m)D-BR.

Since kappa V-B>c and kappa R<=4,

    cL-4U <4D-BR+4V a^m.                              (3)

The term-BR is the improvement inherited from the concurrent466 proof.
The exact U is the improvement retained from Session139 and the514 proof.

## 4. Pay the geometric tail with compatible mass

Let Q_R=R+33/2+8/R; (1) gives4V<=Q_R. We prove

    4V a^m <(c-d)L+2U.                                (4)

**Case m>=13.** Convexity of log(a^(R-L)/L) in L gives

    Q_R a^m/L <=max{[Q_R/(R-13)]a^13, Q_R a^(R-1)}.

Both bounds decrease for integer R>=157. For the first use
Q_R/(R-13)=1+(59/2+8/R)/(R-13); for the second use
Q_(R+1)<Q_R+1 and a(1+1/R)<1. At R=157 both are strictly less than
c-d=7249843/157000000. The larger is about0.028633, below0.046177;
the verifier compares exact fractions. Thus (4) holds without using U.

**Case9<=m<=12.** Since M_h<=m/4, (1) implies

    2U>cR/2-8-8/(2R-1)-m/2.

It suffices to prove positivity of

    F(R,m)=(c-d)(R-m)+cR/2-8-8/(2R-1)-m/2-Q_R a^m.

For each fixed R this is concave in the real variable m, so on[9,12]
its minimum is at an endpoint. Each of F(R,9), F(R,12) increases for
R>=157: its derivative is strictly greater than

    3c/2-d-a^9>0.

The exact endpoint slacks are

    F(157,9)=363349507136067/402563072000000>0,
    F(157,12)=175134030288958737/25764036608000000>0.

This proves(4) in all cases. The proof uses neither independently
maximized incompatible masses nor a finite-n extrapolation.

## 5. One local dual, with both credits retained

For every light group0<w<=t and S>=0, the exact identity

    (wS/2+S^2/4)-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(2t-w)/4

has last term at least tw/4. Heavy groups satisfy W_i>=tS_i. Therefore

    W>=tD-lambda(dL-4U)+(t/4)U.

Use(3) and t/(4lambda)=2 to obtain

    W>hD+lambda RB+lambda[(c-d)L+2U-4V a^m]
     >hD+lambda RB.                                  (5)

As W<=B, this implies(1-lambda R)B>hD. If lambda R>=1 this is
impossible since B,D>=0. Otherwise B<=kappa D implies D>0 and

    h <(1-lambda R)kappa
      <=4/R-4lambda
      <=4/157-4lambda=h,

again impossible. Thus there are no light selected edges. This final
contradiction needs no tiny exponential slack at the critical R=157.

## 6. Equality and boundary

All selected weights exceed t>=kappa. If D>0 the selected-edge invoices
give W>kappa D, contradicting W<=B<=kappa D. Thus D=B=0. Matching-support
stability proves the theorem. Every equality candidate has the same
conclusion, and the inherited binary equality/H=1 classification gives
exactly the stated fair pairs. These attain kappa_n, proving optimality.

The proof newly covers314..465 beyond the concurrent466 result. Its
positive constants fail to certify312: the slot bound no longer forces
nine heavy edges, and the resulting m=8 scalar relaxation fails its tail
test. This is a method hostile, not a counterexample to stability.

The result is about optimal weighted coefficients at vertex orders, not
realization of Hamiltonian-path values. Raw-H spectrum completeness,
all-order weighted stability, and any RH implication remain OPEN.
