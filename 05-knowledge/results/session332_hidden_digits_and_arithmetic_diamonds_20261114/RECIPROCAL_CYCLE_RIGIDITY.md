# A cycle shape at the intersection of two arithmetic descriptions

The [odd-coordinate predecessor note](ODD_COORDINATE_PREDECESSORS.md) identifies the actual minus-sheet \(\{4,11\}\) cycle with \(\{7,21\}\). There is a useful stricter arithmetic question: when can a pair \(H,pH\) be simultaneously an ordinary multiplicative step and a reciprocal two-step affine loop? Here \(p\ge3\) is any odd integer, not necessarily prime.

Use the normalized affine letters
\[
C_{p,j}(X)=\frac{2^j(X+1)-1}{p}.
\]
The \(j=0\) letter is exact division by \(p\) on its integrality domain, so one reciprocal loop is \(H\xrightarrow{j}pH\xrightarrow{0}H\) exactly when
\[
(p^2-2^j)H=2^j-1.
\]
We impose positive odd \(H\) and \(j\ge1\). This generalizes the arithmetic form; it does not claim that a generalized greedy algorithm chooses these letters. For \(p=3\), positivity forces \(j\le3\) and only \(j=3\) gives an integral head \(H=7\). More generally, writing \(D=p^2-2^j>0\), the condition is \(D\mid2^j-1\) and \(H=(2^j-1)/D\). If \(H\) is odd then \(H\equiv7\bmod8\): \(j=1,2\) make \(D\ge5>2^j-1\), and for \(j\ge3\) one has odd \(D\equiv1\bmod8\).

## Even exponents can never make this reciprocal loop

If \(j\) is even, let \(u=2^{j/2}\) and put \(a=p-u,b=p+u\). They are positive odd integers with \(a\le b\), and \(ab=D\). The proposed loop would force
\[
a^2+b^2-4=(4H+2)ab .
\]
Indeed \(b-a=2u\) and \(u^2-1=Hab\). Now fix the even integer \(K=4H+2\ge6\) and apply Vieta jumping. When \(a\ge3\), replace \(b\) by the other root of this quadratic, \(b'=Ka-b=(a^2-4)/b\). It is positive odd and strictly smaller than \(a\). Swapping the variables and repeating must reach \(a=1\). There \(b\mid3\), hence \(b=1\) or \(3\), and the resulting \(K\) is \(-2\) or \(2\), contradiction. Thus **every such reciprocal loop uses an odd exponent**, without a bound on \(p\).

## The head seven is rigid even when the multiplier is free

At \(H=7\) the equation becomes \(7p^2+1=2^{j+3}\). Since \(j\) is odd, write \(x=2^{(j+3)/2}\), obtaining \(x^2-7p^2=1\). The least positive nontrivial unit is \(8+3\sqrt7\): the trials \(p=1,2\) give nonsquares and \(p=3\) gives \(x=8\). Every positive unit solving the equation is a power of this least unit by the elementary unit reduction (divide by the largest power leaving a real value at least one).

The \(x\)-coordinate of an even power is \(2y^2-1>1\) for a previous \(x\)-coordinate \(y\), so it is odd and cannot be a power of two. For an odd power \(r\), the coordinate is \(T_r(8)\), where \(T_r\) is the Chebyshev polynomial with \(T_r((t+t^{-1})/2)=(t^r+t^{-r})/2\). The polynomial \(T_r(X)/X\) is integral for odd \(r\) and its constant term is odd: \(T_{2s+1}'(0)=(2s+1)(-1)^s\), from the recurrence or the cosine identity. Therefore \(T_r(8)/8\) is odd and is greater than one if \(r>1\). Only \(r=1\) works, giving \(p=3,j=3\). This gives two different rigidity statements: for multiplier 3 the only head is 7, and for head 7 the only multiplier is 3 in this affine family.

The same Pell argument gives at most one eligible multiplier/exponent for each fixed head \(H=2^a-1\) with odd \(a\ge3\). For such a head \(j+a\) is even and \(x=2^{(j+a)/2}\) is a coordinate of a Pell unit for \(x^2-Hp^2=1\) (the discriminant \(H\equiv7\bmod8\) is nonsquare). If its fundamental \(x\)-coordinate is odd then all the power coordinates are odd by the recurrence. If it is even, all even powers are odd and every odd-index coordinate is its fundamental coordinate times an odd integer; only the fundamental unit can have \(x\) a power of two. No existence of a suitable fundamental unit is claimed for the other heads.

The [finite divisor replay](replay_reciprocal.py) found no other reciprocal loop through exponent 61, but the odd-exponent, arbitrary-head classification remains open in this note. Even complete classification of this arithmetic geometry would not prove tournament realizability; here forward arithmetic closure is already refuted by \(63\to21\).
