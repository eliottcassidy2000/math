# Tensor reserve, discarded-degree energy and exact slots: sharp stability from66

**Status: PROVED; same-agent proof audit with exact rational certificates and a separate weighted-degree implementation. No independent-agent audit of this new combined theorem is claimed.** This combines the
independently audited incoming tensor matching theorem, the concurrent
discarded-degree estimate, and the retained-baseline/slot method. It does
not modify their frozen proofs or claim an independent-agent audit here.

For a real weighted tournament A, put V=sum_e p_e(1-p_e),
Delta=H-det(I+A)-abs(det A), R=floor(n/2), and

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every integer n>=66, Delta>=kappa_n V with optimal
coefficient. Equality with V>0 is precisely maximal ordinal fair pairs,
with one singleton for odd n. At V=0 it is precisely transitivity.
Together with inherited n=3..9, the remaining interval is10..65 OPEN.

## 1. Inputs and common rounding law

Let B=E Delta(T) for independent binary rounding, W=sum q_e q_f over
unordered incident edge pairs, and let greedy matching select weights w_i
with assigned discarded masses S_i. Put M=sum w_i and D=sum S_i=V-M.
The accepted Session138 bounds are

    B>=W>=2V^2/n-V/4,
    Delta>=B+(1-exp(-4V))/4,
    Delta>=B+kappa_n M.

A retained nonperfect r-edge matching of mass M_r also gives
Delta>=B+a_r M_r, a_r=(4/r)(1-(3/4)^r), decreasing in r.
The disjoint greedy groups have W_i>=w_iS_i and W_i>=S_i^2/4.

The first new input is the [discarded-degree lemma](../stability_priority_moduli_20261102/STABILITY_N500.md),
Section1, valid for every nonnegative weighted graph:

    W>=2D^2/n.                                         (1)

To recall its mechanism, if d_v is the discarded weighted degree, the
discarded/discarded wedges equal(sum d_v^2)/2-sum q_f^2. The selected
edge removing each discarded f has weight at least q_f, and its distinct
selected/discarded wedge pays q_f^2. They cancel the negative term.
Cauchy--Schwarz and sum d_v=2D prove(1).

The second is the independently audited [tensor matching theorem](../session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md):
after rounding outside the greedy matching,

    Delta>=B+1-product_i(1-w_i)-1_(n=2r)product_i w_i
         >=B+1-exp(-M)-4^(-R).                          (2)

The perfect-matching subtraction occurs only when n=2r=2R. All other
rounding corrections are nonnegative and every occurrence of B is the
SAME terminal product distribution. No coordinatewise Pfaffian upper
bound, or false single-edge multiplicative inequality, is used.

## 2. Bootstrap the nonlinear budget to99/100

Assume V>0 and Delta<=kappa V, with kappa=kappa_n and R>=33. Then
kappa<=4/R<=4/33 and kappa<=8/(n-1). The exponential Taylor bounds
exp(4)>50 and exp(8)>1000 are exact finite rational checks. The decreasing
ratio(1-exp(-4V))/(4V) excludes V<=1, since g(1)>49/200>4/33.
Consequently kappa V>49/200, whence

    V>49R/800>=1617/800>2,
    kappa V>B+c_0,            c_0=999/4000.              (3)

Whenever kappa V>B+c, (1) gives

    M=V-D>c/kappa+2D^2/(n kappa)-D
      =c/kappa-n kappa/8
         +[2/(n kappa)](D-n kappa/4)^2
      >=cR/4-1-1/(2R)=:m_c(R).                         (4)

Unlike the weaker degree bound, the asymptotic subtraction here is one
unit of matching mass. Combining(2),(4) bootstraps the constant c.
For every R>=33 the following advances hold:

    999/4000 -> 3/5 -> 19/20 -> 99/100.

At R=33 the mass lower bounds are respectively

    551911/528000,        2597/660,        18011/2640.

For these three values, the Taylor sums of exp(m) through degrees2,3,4
respectively satisfy

    T_k(m_c(33)) [1-c_next-4^(-33)]>1.

All factors are positive and the comparisons are exact rational checks.
Thus1-exp(-m_c(33))-4^(-33)>c_next. Monotonic increase of m_c(R) and
decrease of4^(-R) extend every step uniformly. Hence

    kappa V>B+c,                 c=99/100.              (5)

## 3. Preserve the slots: at least25 heavy edges

Set eta=4/33, and suppose L>=1 selected edges have weight at most eta.
Their mass U satisfies U<=eta L. The retained heavy matching has size r_h
and mass M_h=M-U. Because L+r_h<=R and each selected weight is at most1/4,
(4),(5) yield

    r_h>=4M_h>4[c/kappa-n kappa/8-eta(R-r_h)],
    r_h(1-4eta)>(c-4eta)R-4-2/R=:g(R).                 (6)

Here1-4eta=17/33>0. At R=33,

    g(33)-24(1-4eta)=811/3300>0.

The function g increases for R>=33, because c-4eta>0. Therefore r_h>24,
so r_h>=25. Let m=R-L. Then25<=r_h<=m<=R-1 and n>2m.
No empty retained matching is used. The exact slot relation is what makes
this count compatible with the large light-edge cut eta.

## 4. Matching deletion keeps-BR and pays its tail per light edge

Keep B in the retained matching inequality, replace r_h by m using the
decreasing coefficient, and use M_h>=V-D-eta L. With e=(3/4)^m,

    (kappa V-B)(R-L)>=4(1-e)(V-D-eta L).

Expansion, kappa R<=4 and(5) imply

    L(c-4eta)<4D-RB+4eV.                               (7)

The variance-wedge bound gives

    4V<=n/2+16n/(n-1)<=Q_R:=R+33/2+8/R.

For1<=L<=R-25, convexity of log((3/4)^(R-L)/L) gives

    4eV/L <=max{Q_R(3/4)^(R-1),
                 [Q_R/(R-25)](3/4)^25}<1/100.          (8)

Both endpoint expressions decrease for integer R>=33. The first uses
Q_(R+1)<Q_R+1; the second uses
Q_R/(R-25)=1+(83/2+8/R)/(R-25). At R=33 their values are respectively
less than0.005 and0.0047, and exact rational comparisons with1/100
are supplied by the verifier. Thus, with

    a_0=c-4eta-1/100=817/1650,

equation(7) becomes

    a_0 L<4D-RB.                                      (9)

In particular D>0. Note a_0-4eta=17/1650>0.

## 5. A single shared group budget closes the proof

Write D_l,D_h for discarded masses assigned to light/heavy selected edges.
The disjoint invoices and Cauchy--Schwarz imply

    B>=W>=eta D_h+D_l^2/(4L)>=eta D-L eta^2.             (10)

The last inequality is the square(D_l-2Leta)^2>=0. Combining(9),(10),

    B(a_0-R eta^2)>D eta(a_0-4eta).                    (11)

The right side is positive. If a_0-R eta^2<=0 this is impossible. Otherwise
the full retained matching bound gives B<=kappa D<=4D/R, and

    B(a_0-R eta^2)<=D(4a_0/R-4eta^2)
                   <=D(eta a_0-4eta^2),

because R>=33=4/eta. This contradicts(11), including R=33. Hence L=0.

All selected weights now exceed eta>=kappa. If D>0 then
B>=W>=sum w_iS_i>kappa D, again impossible. Thus D=B=0. The uncertainty
graph is a matching and accepted matching-support stability proves the
desired bound. In an equality candidate B=0 makes every supported binary
completion transitive; the accepted H=1 scalar classification yields
exactly maximal fair pairs, which attain kappa_n.

## Scope and next boundary

The theorem newly covers66..139 beyond the incoming140 result. It uses
the tensor input, stronger discarded-degree input, exact slot accounting,
and the negative baseline together; none alone gives this statement.
At n65 the chosen R>=33 comparison fails; this is a certificate boundary,
not an actual counterexample. Any further reduction must be proved with
new constants/compatibility, not extrapolated from the finite checks.

The separate guarded path and all-modulus theories remain informative but
are not inputs to this cutoff. Raw-H completeness, orders10..65, Collatz,
RH and every other unclosed problem retain their stated OPEN status.
