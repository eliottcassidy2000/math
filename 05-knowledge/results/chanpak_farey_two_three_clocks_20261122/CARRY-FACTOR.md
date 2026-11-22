# Carry-to-factor inverse for two prescribed low words

**PROVED, with finite E61 sectors.** This inverse uses the binary carry to reduce the original digital equation to finitely many exact factorizations. It concerns the same actual Farey marked unit flip, not a modular surrogate. It does not settle unrestricted E61.

Write the lower chart as
\[
 H=k+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+h(f(x)+f(y)),
\]
where \(h\) is odd and \(f(2u+1)=f(u)\), \(f(2u)=f(u)+u+[u\text{ odd}]\). Fix padded suffix words \(i,j\in[0,B)\), \(B=2^\ell\), and write \(x=BX+i,y=BY+j\). For \(\epsilon,\eta\in\{0,1\}\), put
\[
\begin{aligned}
 C_i&=B-1-\operatorname{rev}_\ell(i)-B(1-\epsilon),&
 \phi_i(X)&=f(X)+[i<B/2][X\text{ odd}]+C_iX,\\
 C_j&=B-1-\operatorname{rev}_\ell(j)-B(1-\eta),&
 \phi_j(Y)&=f(Y)+[j<B/2][Y\text{ odd}]+C_jY,\\
 K&=k+qij+ti+cj+r[i\text{ odd}]+e[j\text{ odd}]+h(f(i)+f(j)),\\
 A&=t+qj+h(1-\epsilon),&
 C&=c+qi+h(1-\eta),& P&=qB.
\end{aligned}
\]
The standard exact split of \(f(BX+i)\) now gives, without roundoff,
\[
 H=K+B(PXY+AX+CY)+hE,\qquad E=\phi_i(X)+\phi_j(Y).
\]
Thus every candidate for target \(N\) must have a **single integer carry residue**
\[
 E\equiv h^{-1}(N-K)\pmod B.
\]
Once an actual integer \(E\) is chosen, set \(V=(N-K-hE)/B\). The remaining bilinear equality is equivalent to the factorization
\[
             (PX+C)(PY+A)=PV+AC .
\]
These two factors are positive. Enumerate the positive divisors \(d\) of the right side, put \(X=(d-C)/P\), \(Y=((PV+AC)/d-A)/P\), require nonnegative integrality, and check the original \(f\)-reads. Any nonpositive right side is impossible. The factorization is not a relaxation of the integer part of the equation and no bounded search over the high binary words remains. One can also replace \(B\) in the carry congruence by \(B\gcd(P,A,C)\) before solving the linear congruence, since that larger modulus divides the entire bilinear term. The dyadic residue/division step plays the same role as a specified branch of an affine Collatz inverse: the valuation guard is checked before the integer reconstruction. Here the reconstructed object is a divisor pair, not a Collatz predecessor; no assertion about Collatz orbits is needed.

## Explicit carry intervals

One has \(0\le f(X)+[X\text{ odd}]\le X\), hence
\[
 C_i X\le\phi_i(X)\le(C_i+1)X
\]
and similarly for \(Y\). Given a bound \(X+Y\le Z\), all carries are therefore in the exact finite interval
\[
 Z\min(0,C_i,C_j)\le E\le Z\max(0,C_i+1,C_j+1).
\]
For homogeneous suffixes choose \(\epsilon,\eta\) to match zero versus one: when \(i=0\), \(C_i=-1\); when \(i=B-1\), \(C_i=0\). Then \(\phi_i\) is respectively nonpositive or nonnegative and at most \(X\) in magnitude. The same method allows arbitrary prescribed leading defects among the top \(d\) bits of the low block: then \(|C_i|\le2^d\), so the inverse stays cheap if those defects are confined high in the block. It is not an efficient inverse for an arbitrary unrestricted suffix bank.

For positive \(x,y\) and an aspect ratio within \(R_0\), use \(xy\ge R_0(x+y)^2/(R_0+1)^2\) and \(qxy\le N-k\). Without any ratio bound, let \(D=\lfloor(N-k)/q\rfloor\), \(x_{\min}=i\) if \(i>0\), otherwise \(B\), and similarly for \(y_{\min}\). Reject immediately if these minima have product greater than \(D\). Convexity of \(x+D/x\) on the interval \([x_{\min},D/y_{\min}]\) otherwise gives
\[
 x+y\le\max\{x_{\min}+\lfloor D/x_{\min}\rfloor,\ y_{\min}+\lfloor D/y_{\min}\rfloor\}.
\]
In each case subtract \(i+j\), then divide by \(B\) and take the floor to obtain \(Z\). For pure suffixes with no aspect restriction the number of carries is \(O(1+D/B^3)\), whereas in a fixed aspect sector it is \(O_{R_0}(1+\sqrt D/B^2)\). Thus the valuations can bring an unbounded-aspect sector down to only one or two factorizations.

## E61 replay

At \(N=2^{121}-2^{60}-1\), the script uses the canonical charts for spines 4 and 8. In each row it tests all four independent choices \(i,j\in\{0,B-1\}\), all positive addresses and all widths of the remaining high words:

| Spine | Each address ends in a constant block of at least \(\ell\) bits | Aspect cap | Carries to factor | Divisor pairs in both required classes | Exact hits |
|---:|---:|---:|---:|---:|---:|
|4|26|8|3|0|0|
|8|26|8|0|0|0|
|4|32|none|87|0|0|
|8|30|none|89|0|0|

These are additional absolute sectors, including patterns where both tails are zeros or both ones, for which the sum need not end in ones. The no-aspect rows exclude collapse to highly unequal positive flanks, not only the balanced windows treated earlier; axes are separate. All primality claims used to list divisors are independently verified by recursive Lucas certificates on exact factorizations of \(p-1\). Failure to construct any certificate fails the replay rather than silently treating a probable prime as proven.

For the standard-library-only replay, run the three commands below with the supplied exact factor/witness table; this path does not invoke a factorization package.

    python farey-carry-factor.py 26 8 --certificate farey-carry-certificates.json
    python farey-carry-factor.py 32 0 --spine 4 --certificate farey-carry-certificates.json
    python farey-carry-factor.py 30 0 --spine 8 --certificate farey-carry-certificates.json
