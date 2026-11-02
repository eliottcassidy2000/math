# Shared-budget analytic stability for every order n >= 748

Status: PROVED subject to independent review of this new packet. The only
inputs are the accepted Session138 matching-rounding, nonlinear-margin,
and greedy wedge bounds. No new graph census is used. In particular this
proves the previously unresolved orders 748 through 1499; it does not settle
orders 10 through 747.

Write R=floor(n/2), q_e=p_e(1-p_e), V=sum q_e,
Delta=H-det(I+A)-abs(det A), B=E Delta(T), and W=sum q_e q_f over unordered
incident edge pairs. Put

    kappa_(2R)   = (4/R)[1-(3/4)^R-(1/4)^R],
    kappa_(2R+1) = (4/R)[1-(3/4)^R].

The theorem is K_n=kappa_n for all n>=748. Equality with V>0 consists
exactly of ordinal fair pairs, with one singleton in odd order.

## 1. Accepted inputs and the new mechanism

From Session138's `asymptotic_stability/PROOF.md`,
`stability/MULTISTAR_STABILITY.md`, and `threshold_refinement/PROOF_N1500.md`:

* B>=W and W>=2V^2/n-V/4;
* Delta>=B+g(V), where g(V)=(1-exp(-4V))/4;
* any retained matching of m positive edges and mass M gives
  Delta>=B+kappa_n M, and if n>2m, Delta>=B+a_m M,
  where a_m=(4/m)[1-(3/4)^m];
* descending-weight greedy matching with selected weights w_i has
  assigned discarded-edge masses S_i. For D=sum S_i=V-M,
  W>=D^2/(2n). Its disjoint group wedge contributions satisfy both
  W_i>=w_i S_i and W_i>=S_i^2/4.

There are two improvements over the n>=1500 argument. First, retaining B
in the nonlinear inequality before completing a square controls V-D,
rather than separately bounding V below and D above. Second, heavy and
light discarded edges share one wedge budget; the budget is not spent twice.

## 2. Candidate configurations and a square completion

It suffices to consider V>0 and Delta<=kappa V, with kappa=kappa_n. Always

    kappa R<=4,                  kappa<=8/(n-1)<=8/747.

The exact partial exponential sum sum_(j=0)^7 4^j/j!=16319/315>50 gives
c0=(1-exp(-4))/4>49/200. Since g(V)/V decreases, V<=1 would imply
Delta/V>=c0>kappa. Hence V>1. Then kappa V>=g(V)>49/200, and

    V>49(n-1)/1600>4.

Since exp(16)>50^4>250000, g(V)>c, where we fix

    c=249999/1000000.

The variance-wedge bound gives

    V<=n/8+n kappa/2<=n/8+4n/(n-1)<=n/7,                 (1)

where n-1>=224 suffices. Retaining a full greedy matching gives

    B<=kappa D.                                        (2)

More importantly, the nonlinear bound and W>=D^2/(2n) give

    kappa V>=B+g(V)>D^2/(2n)+c,

and therefore, by completing a square in D,

    V-D>c/kappa+D^2/(2n kappa)-D
       >=c/kappa-n kappa/2.                            (3)

This replaces the prior c/kappa-2n kappa lower bound.

## 3. The number of light selected edges

Set eta=1/55. Suppose L>=1 selected greedy edges have weights w_i<=eta,
and remove those edges. Let r_h,M_h be the size and mass of the retained
heavy matching. Using L<=R<=n/2 and (3),

    M_h>=V-D-eta L
        >c/kappa-n kappa/2-eta n/2
        >=c(n-1)/8-4n/(n-1)-n/110.

Every retained weight is at most1/4, so r_h>=4M_h. In fact

    r_h>n/15                                               (4)

uniformly for n>=748: after subtracting n/15 from four times the last
lower bound, the derivative is

    c/2-2/55-1/15+16/(n-1)^2>0,

and its value at748 is428091991/1494000000>0. Thus r_h>0. Since at
least one matching edge was removed, r_h<=R-L and n>2r_h.

Put epsilon=(3/4)^(n/15). For n>=748,

    epsilon<=1/(2000n).                                  (5)

Indeed n(3/4)^(n/15) has logarithmic derivative
1/n-log(4/3)/15<1/n-1/60<0, because log(4/3)>1/4. The endpoint is the
exact integer inequality

    (2000*748)^15 * 3^748 < 4^748.

Consequently a_(r_h)>=4(1-epsilon)/r_h. Applying the retained-matching
inequality (discarding B>=0 only at this step) and using r_h<=R-L gives

    kappa V(R-L)>=4(1-epsilon)(V-D-eta L).

Rearranging and using kappa R<=4,

    L[kappa V-4eta(1-epsilon)]
       <=[kappa R-4+4epsilon]V+4(1-epsilon)D
       <=4epsilon V+4D.

By (1),(5), 4epsilon V<=delta=1/3500. Since L>=1, delta<=delta L.
Also kappa V>c, so

    L(c-4eta-delta)<4D,
    L<C D,       C=4/(c-4eta-delta)>0.                   (6)

In particular D>0.

## 4. One wedge budget, not two

Split discarded mass as D=D_h+D_l according to whether its selected
edge is heavy or light. Heavy groups contribute at least eta D_h.
Light groups contribute at least sum S_i^2/4>=D_l^2/(4L). These groups
have disjoint counted wedge pairs, hence the single simultaneous bound

    W>=eta D_h+D_l^2/(4L)
      >=eta D_h+D_l^2/(4CD).                            (7)

The last weak inequality is valid also if D_l=0. Put x=D_h/D in[0,1].
The right side divided by D is eta x+(1-x)^2/(4C). We have C eta<1/2,
so its minimum on[0,1] occurs at x=1-2C eta and equals

    h=eta-C eta^2.

The exact endpoint calculations are

    1/2-C eta =2427923/27255846>0,
    h-8/747   =190787/186634405485>0.

Thus W>=hD>kappa D>=B>=W, a contradiction. Therefore L=0.

## 5. Closure and equality

All greedy selected weights now satisfy w_i>eta>8/747>=kappa. If D>0,
at least one assigned group has positive mass, and the simultaneous bound
W>=sum w_i S_i>kappa D contradicts (2). Thus D=0 and then B=0.
The uncertainty graph is a matching. Its accepted matching-support bound
proves Delta>=kappa V, while ordinal fair pairs attain the ratio kappa.

For equality with V>0 the same argument gives D=B=0. The accepted binary
equality theorem (n>=4) then makes every supported completion transitive;
the accepted H=1 scalar minimization gives exactly the maximal fair-pair
ordinal configurations, with the singleton in odd order.

## Verification and scope

`verify_n748.py` checks all rational endpoint signs and the exponential
integer comparison with standard-library integer/Fraction arithmetic.
Uniform extension is proved by the displayed derivatives, not by numerical
sampling. This argument does not require the invalid termwise Pfaffian
unit-coefficient claim and makes no claim at orders10..747.
