# Explicit prime-power kernels and CRT consequences

Status: PROVED. These formulas and derivations were supplied and freshly checked by the independent reviewer of the reconstructed packet.

Let K_ab(q,r)=E[g h 1(g|D+a)1(h|r+b d)], where d,U,D are independent uniform modulo q, g=gcd(d,U,q), h=gcd(D,U,q), and a,b are 0 or 1. Every branch is nonnegative and sums to q over r: summing over r cancels h, and the remaining D-average of g1(g|D+a) is 1. Thus K_ab/q is a probability law. CRT factors each branch separately.

## Derivation at p^k

Let alpha=1-1/p, t=min(v_p(U),k), and v=min(v_p(r),k), with v=k for r=0. The weights of t are alpha p^(-t) for t<k and p^(-k) for t=k.

For K_00, put j=min(v_p(D),t). If j<t, the D probability times h is alpha, while the d-average of g1(g|D) is alpha(j+1). If j=t, that average is 1+alpha t. Requiring j<=v gives the conditional kernel

A_t = 1+t alpha+t(t+1)alpha^2/2 if t<=v;
A_t = (v+1)(v+2)alpha^2/2 if t>v.

For K_10, split d into units and nonunits. The nonunit contribution is 1/p. The unit contribution is alpha E_D[h1(h|r)]. Consequently

B_t = 1+t alpha^2 if t<=v;
B_t = (v+1)alpha^2+1/p if t>v.

For K_01, the contribution at j<t is alpha^2 if j<=v, and alpha p^(v-j) if j>v. The j=t contribution is p^(min(v,t)-t). When v<t the geometric tail equals 1/p, proving K_01=K_10 exactly.

For K_11 and t>=1, nonunit d forces D to be a unit under g|(D+1), contributing 1/p. Unit d contributes alpha^2+1/p if r is a unit, or alpha^2 otherwise. At t=0 the kernel is 1. Since Pr(t>=1)=1/p, its averaged value is 1+p^(-3) on units and 1-p^(-2)+p^(-3) otherwise.

## Closed forms

Write A=K_00, B=K_01=K_10, C=K_11. For 0<=v<k and z=p^(-v-1), summing the preceding finite geometric series gives

A=1+2/p-z[1+2/p+(v+1)(1-p^(-2))],
B=1+1/p-1/p^2-z(1-p^(-2)).

At r=0 modulo p^k,

A_0=1+2/p-p^(-k-1)[2+k alpha],
B_0=1+1/p-1/p^2-alpha p^(-k-1).

For all k>=1,

C=1+p^(-3) if p does not divide r;
C=1-p^(-2)+p^(-3) if p divides r.

Therefore

L_q(r)=(9/16) product A(p^k,r)+(6/16) product B(p^k,r)+(1/16) product C(p^k,r),

where products range over p^k exactly dividing q. The common branch label must be retained; multiplying averaged local factors is incorrect.

For units the prime-power factor is 1-3/(2p^2)+1/p^3. Put

F_infinity=1+3/(2p)-7/(16p^2)+1/(16p^3).

For 1<=v<k,

L_(p^k)(r)=F_infinity-p^(-v-1)/16*[9+18/p+(9v+15)(1-p^(-2))].

For r=0,

L_(p^k)(0)=F_infinity-p^(-k-1)/16*[18+(9k+6)alpha].

At p=3 these reduce to unit factor 47/54, positive-valuation factor 157/108-3^(-v)(v/6+85/144), and zero-tail factor 157/108-3^(-k)(k/8+11/24).

## Strict positive divisibility correlation

At every prime power A_0>B_0>C_0. Indeed,

A_0-B_0=E[t alpha+t(t-1)alpha^2/2]>0,
B_0-C_0=alpha/p*(1-p^(-k))+alpha/p^2>0.

The normalized three branches are product probability laws across distinct prime powers. For divisibility by distinct prime powers p^k and ell^m, let x_i and y_i be the conditional probabilities in branches A,B,C, with weights w=(9/16,6/16,1/16). Both triples strictly decrease, so

Cov(1[p^k|F],1[ell^m|F])
=sum_(i<j) w_i w_j(x_i-x_j)(y_i-y_j)>0.

This is the limiting fixed-length law. For each fixed pair it also implies positivity for sufficiently large length, without a uniform threshold.

## Uniform nonvanishing

The conditional formulas make A and B strictly increasing with the valuation v, including the terminal step to the zero residue. For A, a conditional increment below its cutoff is alpha^2(v+2), and the cutoff increment is 1+t alpha. For B these increments are alpha^2 and alpha. Positive mass above the current valuation gives strict increase after averaging.

Thus

A>=1-2/p^2+p^(-3)>(1-p^(-2))^2,
B>=1-1/p^2+p^(-3)>1-p^(-2),
C>=1-1/p^2+p^(-3)>1-p^(-2).

Using product_(p odd)(1-p^(-2))=8/pi^2 in the three-branch formula gives, for every odd q and every residue r,

L_q(r)>36/pi^4+7/(2pi^2)=0.7241995039168183... .

For q=1 this follows from L_1=1. This is a uniform bound on limiting modular densities only. It provides neither uniform convergence in q, integer coverage, nor a height-truncated result.

## Fresh independent checks

Independent exact rational enumeration verified all four branch formulas at every residue for q=3,5,7,9,25 and branchwise CRT at q=15; normalization was checked for each branch. The fresh q=15 example is L_15(0)=5179/3375 versus L_3(0)L_5(0)=5134/3375, difference 1/75. These new calculations and their scope are documented in AUDIT.md.
