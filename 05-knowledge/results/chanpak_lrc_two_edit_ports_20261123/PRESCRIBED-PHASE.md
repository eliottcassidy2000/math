# Prescribed cofactor phases and the common content of the fork

**PROVED**, with the same literal filler and five-edit route as [the original wheel](WHEEL.md). The cofactor is prescribed modulo the wheel, not as an absolute integer.

Let \(D=H(M_{\rm fill})\), \(Q=Q_2(M_{\rm fill})\), and \(c=40Q+D\). Keep the wheel modulus \(M=9\prod_{5\le\ell\le y}\ell\) and prescribe any unit \(r\bmod M\). To choose the physical filler residues, set \(D=r-40Q\). In addition to \(c\), keep these six values free of the wheel primes:
\[
 D,\quad 3r-68Q,\quad (5r-148Q)/3,\quad
 67r-2368Q,\quad 19r+80Q,\quad 23r-80Q .
\]
The displayed expressions with \(r\) are local evaluations: the actual integer forms use \(c\), not a chosen integer representative of \(r\). We may also require \(Q\) to avoid all **odd** wheel primes; \(Q_2\) is even, so no condition at two is imposed on \(Q\).

At an odd \(\ell\ge11\), the six numerator forms and \(Q\) forbid at most seven points: each form has a coprime coefficient pair and \(r\) is a unit. The roots dilate by \(r\). At seven the six forms exclude \(Q/r=2,3,5\); adding zero leaves three choices. At five the excluded set with or without the condition on \(Q\) is \(\{0,1,4\}\). At three impose \(Q\equiv2r\bmod3\), so \(Q\equiv D\bmod3\) in the even-order phase and both are units. Among the three lifts to mod nine, delete the unique one with \(5r-148Q=0\bmod9\). The other reads are units mod three. The normalized physical reads and \(D\) are automatically odd.

These choices combine by CRT, with a uniform correction length as in [the original growing wheel](WHEEL.md#a-growing-rigorous-wheel-in-every-large-scale-window). For every sufficiently large \(X\), uniformly for \(r\in(\mathbb Z/M\mathbb Z)^\times\), the physical fork parent has
\[
 h=7c\equiv7r\bmod 7M,\qquad
 h=X+O(X^{1-\frac12\log_3(3/2)+o(1)}).
\]
In particular \(r=p\bmod M\) is legal for a prime \(p>y\) at scale \(X=7p\). For any fixed \(\epsilon>0\) less than \(\frac12\log_3(3/2)\), applying the same construction at \(X=7p\pm(7p)^{1-\frac12\log_3(3/2)+\epsilon}\) gives two such parents on strictly opposite sides of \(7p\). Use the cutoff defined at the smaller scale in both constructions (the correction bound stays valid). Both satisfy the same congruence and wheel conditions. Straddling is at the error-window scale, not at one modulus or one arc.

## Extract one gcd, and all eight readings separate

View the eight distinct readings as linear forms in \((c,Q)\) with rows
\[
 (1,-40),(3,-68),(5,-148),(67,-2368),(19,80),(23,-80),(1,0),(0,1),
\]
with the third row divided by 3. No pair of these rows is proportional. The prime divisors of the pair determinants of the integer rows are precisely
\[
 2,3,5,7,11,13,17,19,23,31,37,67,73,331,383,751,1049 .
\]
In particular outside these primes the seven forbidden clock residues are all distinct. For the exceptional primes \(\ge5\), the counts are 3 at 5; 4 at 7 and 13; 5 at 11 and 37; and 6 at the remaining exceptions.

If \(y\ge1049\), **every pair of the eight integer readings has gcd exactly \(G=\gcd(D,Q)\)**. As \(c=D+40Q\), this is also \(\gcd(c,Q)\). The wheel makes \(G\) a unit at 3; hence it divides the integer quotient in the third form as well. The eight values divided by \(G\) have no common prime in any pair: at a wheel prime this follows from local avoidance (at 2 only \(Q\) may be even), and outside the wheel an invertible two-row matrix would force that prime to divide both \(c/G,Q/G\). Thus any shared obstruction across the normalized route is exactly common filler content, not a pair-specific prime.

This does not prove \(G=1\). For example at any fixed wheel depth we can add a prime \(P>\max(y,1049)\) to the independent controller and prescribe \(D=Q=0\bmod P\), without changing the wheel phase. Nor does dividing by \(G\) itself construct tournaments for the divided heights.
