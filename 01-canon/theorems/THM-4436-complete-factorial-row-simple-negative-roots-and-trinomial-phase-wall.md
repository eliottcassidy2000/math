---
id: THM-4436
title: "Complete factorial-row simple negative roots and trinomial phase wall"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 third research wave
---

# THM-4436 -- Complete factorial-row simple negative roots and trinomial phase wall

**PROVED ANALYTICALLY RELATIVE TO A CITED FINITE-PRESERVER THEOREM;
FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof, precise primary-source hypothesis, Laurent consumer and controls](../../05-knowledge/results/nc2_channel_realrooted_overnight_hexagon_sep05.md)
are part of this theorem.

Let A,B,h,r,z,x be integers with
0<A<B,h>=0,x>=0,0<=r<B,0<=z<A. Write delta=B-A and
m=x+B*h+r+z. The polynomial

```text
P(t)=sum_(j=0)^h m!*t^j/
       [(x+delta*j)! (B*h+r-B*j)! (z+A*j)!]
```

is nonzero constant for h=0, and for h>=1 has exactly h distinct strictly
negative real roots. Neither gcd(A,B)=1 nor a Laurent charge equation is
required.

The proof splits its normalized coefficients into a binomial row times
B-1 falling-factor and B-1 inverse-rising-factor sequences. Their finite
symbols are reflected or reversed Laguerre polynomials with parameters>-1.
The precise cited input is Borcea--Branden,
[*Polya--Schur master theorems*, arXiv:math/0607416v6, Theorem2(b)](https://arxiv.org/pdf/math/0607416#page=5):
each bivariate symbol is a product of positive linear forms and hence real
stable, so its diagonal operator preserves real-rootedness on the **full**
degree-bounded real coefficient space. Invertibility and openness preserve
simplicity after the first strictly rooted Laguerre factor. The factorial
identities, Laguerre zero proof, simplicity step and consumers are derived
in the linked proof; no additional external multiplicity result is assumed.

For f(u)=alpha*u^(-a)+beta*u^b+gamma*u^c,
a>=1,1<=b<c,alpha*beta*gamma!=0, without content normalization, put

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g,
tau=alpha^(B-A)*gamma^A/beta^B.
```

At **every positive mass** with nonempty support-return fibre, CT(f^m)
is a nonzero coefficient monomial times one of these complete P(tau).
Its cancellations occur only at finitely many strictly negative real
values of tau, each a simple scalar root. Each individual moment's
zero set on the coefficient torus is reduced and smooth, and all scalar
roots are attainable by actual nonzero coefficients.

Off that negative ray, every nonempty moment is nonzero; the first nonzero
moment is exactly the first support return, at most
(a+b)/gcd(a,b)<a+c. The tau quotient preserves common coefficient scaling
and Laurent-variable scaling. It retains every carry channel.

The same weak factorial theorem applies to the proportional central-resonance
polynomials of [THM-2021, Legendre finite recurrence](THM-2021-gmc2-legendre-finite-recurrence-closure.md)
with arbitrary positive p,q:

```text
sum_j m!*kappa^j / [(p*j)! (q*j)! (m-(p+q)*j)!].
```

Use A=q,B=p+q,x=z=0,h=floor(m/B),r=m mod B. This extends
the old symmetric zero-location mechanism, not its recurrence or common-root
claims.

Individual simple negative roots do **not** imply coprimality between
different masses. The sharp general trinomial return problem on the
remaining negative ray is OPEN. Deleting a middle channel can destroy
real-rootedness:5+30t+10t^2 passes, while5+10t^2 fails.

The independent exact checker retains1,530 unfiltered factorial rows,
1,045 charge-derived nonempty rows,112 Laguerre symbol pairs,30 wide
rows and304 proportional-resonance controls;27,847 gates pass normally and
with optimization. Root and independent
agent proof/source audits passed. The recovered finite composition operation
comes from [THM-2760's Faber-flux thread](THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion.md); it is the actual coefficient
map, not an asserted equivalence between the original mathematical objects.
