# Exact stability for n>=514 by paying the tail with light mass

**Status: PROVED; same-session proof audit, exact arithmetic, and a separate local-energy implementation. No independent-agent audit is claimed.** This extends the
independently audited Session139 n>=526 argument. No graph census enters.

Let A be a real weighted tournament, q_e=p_e(1-p_e), V=sum q_e,
Delta=H-det(I+A)-abs(det A), R=floor(n/2), and

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every n>=514, Delta>=kappa_n V, with optimal coefficient.
Equality at V>0 is precisely an ordinal sum of R fair pairs and, when n is
odd, one singleton. At V=0 equality is precisely the transitive case.

## Accepted mechanism and new coordinate

Read Session139's [proof](../session139_incidence_moduli_and_stability_20261102/stability/PROOF_N526.md)
with its [independent audit](../session139_incidence_moduli_and_stability_20261102/independent_threshold526/AUDIT.md).
All Session138 inputs remain as there: B=E Delta(T) for the SAME independent
binary rounding; W is the incident edge-pair variance sum; the greedy
matching has weights w_i and assigned discarded masses S_i. Write D=sum S_i.

    B>=W>=D^2/(2n),    W>=2V^2/n-V/4,
    Delta>=B+(1-exp(-4V))/4,
    Delta>=B+kappa_n(V-D).

A retained r-edge matching of mass M has Delta>=B+a_r M if n>2r,
where a_r=(4/r)(1-(3/4)^r) decreases with r. Disjoint greedy groups have

    W_i>=w_i S_i,       W_i>=w_i S_i/2+S_i^2/4.        (1)

Session139 uses the exact total light mass U in the matching bound, but
discards a positive U term in the final local dual. Here that term pays
the geometric tail precisely when the retained matching could be small.

## 1. Uniform bootstrap and positive retained matching

Assume V>0 and Delta<=kappa V, where kappa=kappa_n, n>=514, R>=257.
Then kappa<=4/257, kappa<=8/(n-1), and kappa R<=4. The inherited bootstrap
from exp(4)>50 gives V>1, kappa V>49/200, and

    V>49(n-1)/1600>4.

Thus with c=249999/1000000 we have kappa V>B+c. Exactly as in Session139,

    B<=kappa D,
    V-D>c/kappa-n kappa/2,
    V<=n/8+4n/(n-1).                                 (2)

Choose

    d=2491/10000, t=d/8, lambda=d/64, h=d/16.

Suppose L>=1 selected edges have w_i<=t; call their total weight U.
The retained heavy matching has size r_h and mass M_h=V-D-U.
Since U<=tL<=tn/2 and each selected weight is at most1/4,

    r_h>=4M_h>c(n-1)/2-16n/(n-1)-dn/4>16.            (3)

At n=514 the last residual over16 is86393731/1026000000>0;
its derivative is125449/2000000+16/(n-1)^2>0. Hence r_h>=17.
Put m=R-L. Then m>=r_h>=17 and m<=R-1; in particular n>2m.
Keeping this heavy matching yields, with a=3/4,

    cL-4U < 4D+4V a^m.                              (4)

This is the same rearrangement as Session139, with no tail estimate yet.

## 2. A tail bound that includes the actual light mass

We prove for every such candidate and every17<=m<=R-1 that

    4V a^m < (c-d)L+2U.                             (5)

Set Q_R=R+33/2+8/R. Equation(2), maximizing over n=2R,2R+1, gives4V<=Q_R.

**Case m>=32.** Since L=R-m, convexity of
log(a^(R-L)/L) in L shows

    Q_R a^m/L <= max{[Q_R/(R-32)]a^32, Q_R a^(R-1)}.

Both terms decrease for integer R>=257: Q_(R+1)<Q_R+1 gives the second
claim, and Q_R/(R-32)=1+(97/2+8/R)/(R-32) gives the first.
At R=257 both are strictly below c-d=899/1000000. The larger is

    [Q_257/225](3/4)^32
      =5789452743369435231/47408132269433547653120
      <899/1000000.

Thus the tail is already smaller than(c-d)L in this case.

**Case17<=m<=31.** By M_h<=m/4 and(2),

    U>c/kappa-n kappa/2-m/4,
    2U>cR/2-8n/(n-1)-m/2
       >=cR/2-8-8/(2R-1)-31/2.                    (6)

The last expression exceeds Q_R a^17 for every R>=257. At257 its surplus
is463842995492257303903/70781597908992000000>0. Its derivative minus the
derivative of Q_R a^17 is

    c/2+16/(2R-1)^2-(1-8/R^2)a^17 > c/2-a^17>0.

Consequently2U>Q_R a^17>=4V a^m. This proves(5) in the remaining case.
The two cases show why a small retained matching must carry enough
light mass to pay its larger tail; treating their worst cases independently
loses this compatibility.

## 3. Retain the positive term in the local dual

For a light group0<w<=t, Session139's exact identity is

    (wS/2+S^2/4)-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(2t-w)/4.

Because w<=t, its final term is at least tw/4. For heavy groups use
W_i>=w_iS_i>=tS_i. Summing the disjoint group contributions gives

    W>=tD-lambda(dL-4U)+(t/4)U.                     (7)

Equation(4) implies

    dL-4U < 4D+4V a^m-(c-d)L.

Substitution into(7), using t/(4lambda)=2, yields the strict bound

    W>(t-4lambda)D
       +lambda[(c-d)L+2U-4V a^m]
     >hD.                                         (8)

Here strictness holds even if D=0. But

    h-4/257=187/41120000>0,

so(8) contradicts W<=B<=kappa D. There are no light selected edges.

## 4. Closure, equality and scope

Every selected weight now exceeds t>kappa. If D>0, (1) gives
W>=sum_i w_iS_i>kappa D, again impossible. Therefore D=B=0.
The uncertainty graph is a matching; the inherited matching theorem
proves the target inequality. The same argument applies to equality.
B=0 means every supported binary completion is transitive, and the
accepted H=1 classification gives precisely the maximal fair-pair
configurations. They attain kappa_n, proving sharpness.

This resolves514..525 beyond the concurrently published n>=526 theorem.
Together with inherited orders3..9, the remaining interval is10..513.
There is no inference about universal raw-H spectrum completeness.

The pointwise dual still has a meaningful boundary: a hypothetical light
group with w tending to0 forces lambda*d>=t^2. Therefore the coefficient
t-4lambda can never exceed d/16 in this particular linear dual. With
d approaching1/4, its ideal ceiling is1/64. This is a relaxation boundary,
not a counterexample at order512 or a reason to abandon richer group data.
