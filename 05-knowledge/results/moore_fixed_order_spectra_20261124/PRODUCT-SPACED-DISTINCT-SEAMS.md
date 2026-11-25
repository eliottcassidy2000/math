# Many distinct physical +2 pairs, and an all-Mersenne obstruction in a decodable subfamily

This strengthens the [two-tail description count](TAIL-LOCAL-LAWS.md) to a theorem about **distinct** counts. It also isolates an exact arithmetic obstruction that is invisible to every fixed modulus.

Use either established variable-tail guard, including its transpose. Write
\[
 G=(H-k)/2=a vu+bv+cu-h(q(v)+q(u)),\qquad
 (h,a,b,c,k)=(29,640,191,158,111)\ \text{or}\ (651,17616,6768,2691,2313).
\]
Transposition exchanges \(b,c\), and the actual guarded arc changes \(H\) to \(H+2\). Here \(0\le q(n)\le n/2\) and \(q(2n)=q(n)\).

## 1. Force the product to be an exact readout

Fix two dyadic scales \(X_0,Y\ge1\). Choose arbitrary
\[
 X_0\le x<2X_0,\qquad Y\le u<2Y,\qquad
 S=h(X_0+Y),\qquad
 T=\text{least power of two strictly larger than }2cY+S,
\]
and set \(v=Tx\). Then \(q(v)=q(x)\) and
\[
 G+S=T\,x(au+b)+\bigl(cu-h(q(x)+q(u))+S\bigr),\qquad
 0\le cu-h(q(x)+q(u))+S<T.
\]
Therefore the height reveals **exactly** \(\lfloor(G+S)/T\rfloor=x(au+b)\). It is essential to include the known offset \(S\); the earlier use of a shorter positive-remainder spacing is a different subfamily. Different product quotients cannot give equal heights in the same box. Given the target and a box, enumerate divisors \(x\) of this product, set \(u=(P/x-b)/a\), and check the actual integer identity. No congruence solution is counted as a height.

## 2. A distinct-output lower bound, even with a fixed congruence filter

Take \(X_0=Y\), so \(T\) is a fixed constant times \(Y\) and the heights are \(O(Y^3)\). There are \(Y^2\) descriptions. Write \(y=au+b\); both factors now have size comparable to \(Y\). In an equality \(xy=x'y'\) put \(x=d r,\ x'=d s\) with \(\gcd(r,s)=1\), so \(y=s t,\ y'=r t\). For \(m=\max(r,s)\), there are at most \(2m\) ordered pairs \(r,s\), and \(O(Y/m)\) possibilities each for \(d,t\). Summing gives \(O(Y^2\log Y)\) ordered collisions, even if the second factor is allowed to be any integer in its comparison interval instead of in the progression. Cauchy–Schwarz gives \(\Omega(Y^2/\log Y)\) distinct product quotients and hence distinct lower heights. Choosing a dyadic \(Y\) comparable to a sufficiently small constant times \(X^{1/3}\) gives
\[
 \boxed{\#\{\text{distinct literal pairs }H,H+2\text{ with }H+2\le X\}
       =\Omega\!\left(X^{2/3}/\log X\right)}
\]
as a lower bound, separately for each guard.

This lower bound persists inside **every fixed admissible lower-height residue class**. The \(x,u\) bits mix independently from their dyadic leading prefix. At primes dividing \(h\), \(T\) is a unit (these primes are odd), so the limiting flank residues are still uniform. At all other primes the conditional uniformity of \(q(x)+q(u)\) absorbs the target, including the normalized dyadic part. There are finitely many relevant residues of \(T\) for a fixed modulus, so the positive survival fraction follows uniformly from the same proof as the local law. Apply Cauchy–Schwarz to the surviving descriptions. This is not positive density of all odd integers and does not establish a growing-cutoff prime sieve.

## 3. Why this product-decodable subfamily misses Mersenne perfect-number neighbors

In this precise spaced subfamily, **neither** \(E_p-1\) nor \(E_p+1\) occurs at either endpoint for any \(p\ge15\) for which \(P_p=2^p-1\) is prime, with \(E_p=2^{p-1}P_p\). This includes \(p=61\) and holds with unequal dyadic scales and either chart marking. It does **not** exclude these targets in the full two-tail families, or in other tournaments.

Suppose a lower query \(H=E_p+\delta\), \(\delta\in\{-3,-1,1\}\), occurs. Then \(G=2^{p-2}P_p+(\delta-k)/2\). Put \(D=\) the least power of two strictly larger than \(2(c+h)\). In the four chart orderings \(D\) equals \(512,512,8192,16384\), and always \(a>D\). The defining argument for \(T\) is at most \(2(c+h)X_0 Y<D X_0Y\) since both scales are at least one. It is at least \(T/2\) by minimality of \(T\), so \(X_0Y>T/(2D)\).

If \(T\ge2^{p-1}\), then \(T\ge D\), and the dyadic scales force \(X_0Y\ge T/D\). The exact product quotient then yields \(G+S\ge T a X_0Y\ge(a/D)T^2\). But \(G<2^{2p-2}\le T^2\) and \(S<T\), so \(G+S<T^2+T\), contradiction since \(a/D-1>1/T\) at the stipulated exponents. Thus \(T\le2^{p-2}\) and divides the main half-perfect term. In all charts \(S\ge2h\ge(k-\delta)/2\). The added remainder \(S+(\delta-k)/2\) is now in \([0,T)\). The exact quotient becomes
\[
 x(au+b)=\frac{2^{p-2}}{T}\,P_p.
\]

This product cannot contain the prime \(P_p\). In fact \(x<2X_0<2T/h<P_p\). For the second factor, \(Y\le T/D'\) with \(D'\) the least power of two strictly larger than \(2c+h\), because \(T\) and \(Y\) are dyadic. Thus \(au+b<2aT/D'+b\). Here \(D'\) equals \(512,512,8192,16384\) in the same order; the bound is \(<P_p\) in both six-site orderings and in the transposed eight-site ordering for \(p\ge15\). In the remaining order it is \(<2P_p\), but \(a,b\) are both even, so \(au+b\) cannot equal the odd prime \(P_p\) or any of its positive multiples below \(2P_p\). The contradiction proves the exclusion.

A fixed modulus sees no difficulty here, by §2, but the actual height fixes an impossible product divisor. This is a specific global failure of the construction, not evidence that the target is absent from the broader chart. The replay tests 5 dyadic boxes per guard and checks the exponent-61 box bounds with an exact primitive-root certificate for the large Mersenne prime:

    python3 product-spaced.py
