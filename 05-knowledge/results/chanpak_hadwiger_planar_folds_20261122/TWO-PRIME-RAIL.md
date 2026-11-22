# A crossed unit rail with simultaneous dyadic and ternary phases

Use the [adjacent-opposite construction](TWO-APEX.md) with the far tags
\(W=11010\,1^k\) on the left and the reverse complement of \(W\)
on the right; both nearer common tags are zero. The address
difference is one, so this is a strong literal unit seam, of
order \(2k+16\), for every \(k\ge1\). Its exact midpoint with
\(X=2^k,\ Z=3^k\) is
\[
 M_k=2910-59246X+18304Z+459162X^2-1498068XZ+3312738Z^2.
\]
To certify the polynomial it suffices to count \(k=1,\ldots,6\)
by the 12-state absolute-height transfer. The two growing
blocks each occupy at most three runs because every return
crosses an apex. Thus the midpoint is a linear combination of
\(1^k,2^k,3^k,4^k,6^k,9^k\) by the ordered-run filling formula;
the six nonzero distinct bases give an invertible Vandermonde
system. The companion replay counts both orientations at each
of these widths.

**Simultaneous phase theorem.** For every \(s\ge5,\ t\ge1\),
an even \(k\) can be chosen to make \(M_k\) equal any prescribed
residue modulo \(2^s3^t\) that is divisible by 32 and equal to 1
modulo 3. In particular \(M_k\equiv E_{61}\pmod{2^s3^t}\) at
every such pair of precisions. This is an existence theorem
for simultaneous local phases, not equality with the target.

Indeed write \(k=2m\). At \(X=0\) the midpoint polynomial \(P\)
satisfies
\[
 P(1+8u)/32=104186+1660945u+6625476u^2 .
\]
Its difference at any two arguments is their difference times
an odd integer, so it permutes each modulus \(2^j\).
The base \(Z=9^m\) generates the residues 1 modulo 8 at every
dyadic precision; one solves the equation for \(m\) modulo
\(2^{s-3}\). At \(Z=0\) the ternary polynomial is
\[
 Q(X)=2910-59246X+459162X^2.
\]
We have \(Q(1)\equiv Q'(1)\equiv1\pmod3\).
The difference \(Q(x)-Q(y)\) for \(x,y\equiv1\pmod3\) is
\((x-y)(-59246+459162(x+y))\), whose second factor is a unit
at 3. Hence \(Q\) permutes the residues 1 modulo 3 for every
modulus \(3^t\), with \(X=4^m\) generating that subgroup:
\(v_3(4^{3^j}-1)=j+1\). This determines \(m\) modulo \(3^{t-1}\).
Use the Chinese remainder theorem on these two coprime moduli,
then add a common period to arrange \(k\ge\max(s,t)\).
Now \(X=0\pmod{2^s}\) and \(Z=0\pmod{3^t}\), so both phases hold
for the full midpoint. The construction gives arbitrarily
large such \(k\) for each fixed phase.

**A third prime rules out actual equality.** All five nonconstant
coefficients of \(M_k\) are divisible by 11, and \(2910\equiv6
\pmod{11}\). Consequently \(M_k\equiv6\pmod{11}\) for every width,
including odd ones. On the other hand \(2^{10}\equiv1\pmod{11}\)
gives \(E_{61}\equiv1\pmod{11}\). Thus this one explicit family
matches the target to arbitrary and independent 2-adic and 3-adic
precision, but a single fixed obstruction at 11 prevents an
exact exponent-61 seam. Other tournament constructions are not
restricted by this obstruction.

**No Euclid midpoint at any exponent.** The obstruction at 11 is
specific to particular exponent classes; 3 and 5 give a uniform
obstruction by instead exploiting the shared width. For \(k\ge1\),
\(M_k\equiv X\equiv(-1)^k\pmod3\), whereas \(E_p\equiv1\pmod3\)
for odd \(p\) and \(0\pmod3\) for even \(p\). Equality therefore
requires both \(k\) and \(p\) to have the following parities: \(k\)
even, \(p\) odd. For even \(k\) we have
\(X=Z=(-1)^{k/2}\pmod5\), and substitution gives
\(M_k\equiv2+3(-1)^{k/2}\in\{0,4\}\pmod5\).
For odd \(p\) the Euclid midpoint is \(1\) or \(3\pmod5\), according
as \(p\equiv1\) or \(3\pmod4\). These sets are disjoint.
Thus the whole family is ruled out modulo 15, despite its full
simultaneous coverage at 2 and 3 within the stated residue class.
