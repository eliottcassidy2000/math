---
id: THM-4432
title: "Trinomial two-channel two-rung noncancellation with carries"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 second research wave
---

# THM-4432 -- Trinomial two-channel two-rung noncancellation with carries

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full proof, polynomial certificate and independent controls](../../05-knowledge/results/nc2_two_rung_overnight_hexagon_sep05.md)
are part of this theorem.

Let f(u)=alpha*u^(-a)+beta*u^b+gamma*u^c with a>=1,1<=b<c,
gcd(a,b,c)=1 and all three coefficients nonzero. Suppose its first
nonempty support-return fibre has exactly two channels. Put

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g, delta=B-A.
```

The first support return is g, and its two count vectors are
v=(x,B+r,z), w=(x+delta,r,z+A), where
0<=r<B,0<=z<A,x=g-B-r-z>0 and a=A(B+r)+Bz.
With d=(delta,-B,A), the **complete** doubled fibre is

```text
2v+j*d,   -floor(2z/A)<=j<=2+floor(2r/B).
```

Thus there may be three, four or five channels. No no-carry assumption or
free all-level semigroup is used.

Write p=multinomial(v),q=multinomial(w), and C,D,E for the multinomials
at2v,v+w,2w. Elementary central-binomial ratios give

```text
Delta=C*q^2-D*p*q+E*p^2<0.
```

At cancellation of the first moment, the coefficient-monomial ratio is
Y/X=-p/q<0. Both optional doubled carry terms have negative sign there.
Consequently the second moment divided by X^2 is strictly negative and
cannot vanish. The full proof gives the literal polynomial identity

```text
q^4*X*Y*T2-Q3*T1=K*X^4,
K=-p*q*Delta+Lminus*q^4+Lplus*p^4>0.
```

Here Lminus,Lplus are the actual optional carry multinomials (zero if
absent), and Q3 is displayed explicitly in the proof. This certifies
noncancellation over the entire coefficient torus without numerical roots.

The first nonzero moment is therefore exactly g or2g, both attainable,
and2g<=a+c. A family with A=n,B=n+1,r=1,z=0,g=(n+1)^2,n>=3
has unbounded smaller endpoint and unavoidable later carries; the conclusion
is not merely a bounded-endpoint or free-semigroup result.

Normal/optimized checks pass76,923 gates on3,665 parameter rows, retaining
all four carry classes and3,530 nonfree cases. The hostile(-13,1,8) has
two first and five second channels. Root and independent agent proof reviews
passed. Concurrent empty-core independently obtained the same two-channel
advance; the overlap is credited, not counted as separate novelty.
The first fibre with three or more channels remains outside this theorem.
