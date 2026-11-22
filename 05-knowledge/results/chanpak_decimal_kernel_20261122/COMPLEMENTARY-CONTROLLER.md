# One suffix controller for all primes at a fixed sum of the two addresses

**PROVED.** In any of the eight canonical Farey charts let \(H(x,y)\) be the lower reading and \(H(x,y)+2\) its marked neighbor. Fix \(D=2^aM\), \(M\) odd, and an odd target \(T\). There is a common suffix length \(\ell=\ell(D)\) with the following property. Choose a positive total \(S\) satisfying
\[
 S+1\equiv0\pmod {2^{\ell+a+1}},\qquad
 tS\equiv T-k-r\pmod{\gcd(h,M)} .
\]
For **every** pair of nonnegative prefix addresses \(x_0,y_0\) with \(x_0+y_0=(S+1)/2^\ell-1\), one can append an \(\ell\)-bit word to \(x_0\) and its bitwise complement to \(y_0\) so that
\[
       x+y=S,\qquad x\text{ odd},\qquad
       H(x,y)\equiv T\pmod D.
\]
The congruences on \(S\) specify a nonempty progression because \(\gcd(t,h)=1\) on each chart. This controls the physical flip at a fixed address sum, including powers of 2, repeated odd primes and arbitrary products. The downstream cube-root placement uses the same length for all prefixes. [Single-address controller and canonical chart conditions](LOCAL-WORDS.md) · [Absolute placement](CUBE-ROOT-PHASE.md)

## Complementary digits and reversed words

Let \(F(x)=f(x)/2\), \(Z=F(x)+F(y)\), \(s=x+y\), and
\[
           C=Z+\lfloor x/2\rfloor-2F(x).
\]
When a bit \(b\) is appended to \(x\) and \(1-b\) to \(y\), the new sum is \(2s+1\). The updates of \(F(x),F(y)\) add \((1-b)\lceil x/2\rceil\) and \(b\lceil y/2\rceil\), respectively. For odd \(s\),
\[
                      C'-C=b(s+1)/2 .
\]
Indeed \(\lfloor(2x+b)/2\rfloor-\lfloor x/2\rfloor=\lceil x/2\rceil\), and \(\lceil x/2\rceil+\lceil y/2\rceil=(s+1)/2\). Across a word \(W\) of fixed padded width \(w\), this adds \((s+1)\operatorname{rev}_w(W)/2\), since the factor \(s+1\) doubles at each step. The initial sum is \(-1\bmod 2^{a+1}\), so \(C\) is fixed modulo \(2^a\) throughout the construction. The low \(\ell\) bits of \(S\) are all ones, giving \(x+y=S\) for every allowed choice of \(W\).

## Switches with equal addresses in both reading directions

For \(M>1\) choose an even \(L\ge4\) with \(2^L=1\bmod M\); for \(M=1\) take \(L=4\). Write \(B=2^L\). A switch appends \(UV\) or \(VU\) to \(x\), where \(U,V\) are padded \(L\)-bit words, followed in either case by \(a+1\) zero bits. Their complements go to \(y\). Initialize the entire suffix with \(a+1\) zero bits on the left. Thus the incoming \(x\) at each switch is \(0\bmod 2^{a+1}\). In each switch below \(U,V\) end in the same bit. The two choices therefore have identical terminal ordinary address residues modulo \(2M\) on both tapes, even before the common tail. Their reversed-word residues agree modulo \(M\) as well, since reversal changes the order of the two reversed blocks and \(B=1\bmod M\).

Take the difference of register readings in the order \(UV\) minus \(VU\). At fixed incoming state the extra increment in \(C\) agrees modulo \(M\) by equality of the reversed residues. Also \(\lfloor x/2\rfloor\) agrees modulo \(M\) by equality modulo \(2M\). Hence the difference in terminal \(Z\) equals the difference in \(2F(x)=f(x)\) modulo \(M\). To obtain the latter, append one \(L\)-bit word \(V\) to an address \(u\) using
\[
       f(Bu+V)=f(u)+f(V)+(B-1-\operatorname{rev}_L(V))u
                   +[\text{first bit of }V=0][u\text{ odd}].
\]
The two-word dependence on the incoming address cancels modulo \(M\) by the equal reversed residues, and the first incoming parity is even. Consequently the difference before the common tail is
\[
 V\operatorname{rev}_L(U)-U\operatorname{rev}_L(V)
   +[\text{first bit of }V=0][U\text{ odd}]
   -[\text{first bit of }U=0][V\text{ odd}]\pmod M .
\]
The common tail preserves the difference because the preceding addresses agree modulo \(2M\).

* With \(U=2,V=4\) the change in \(Z\) is \(3/4\bmod M\).
* With \(U=1,V=B/2+1\), the bilinear part is \(B^2/4-1\) and the boundary term is \(-1\). The change is \(-7/4\bmod M\).

Here division by 4 uses its inverse at the odd modulus. Join two switches of the first type and one of the second type, toggling all three together. The difference in \(Z\) is \(-1/4\bmod M\), independent of both prefixes. Between switches an existing \(Z\)-difference is transmitted with coefficient one because their address residues agree.

For control at \(2^a\) as well, toggle \(2^a\) copies of this triple as one batch. Each switch begins at \(x=0\bmod 2^{a+1}\), and its increment to \(F(x)\bmod 2^a\) depends only on that address residue and the chosen word. Thus each triple has the same difference in \(F(x)\bmod 2^a\). A whole batch contributes zero difference there and \(-2^a/4\bmod M\) to \(Z\); the latter is a unit when \(M>1\). Install \(M-1\) batches. Relative to the baseline of every \(VU\) choice, toggling \(j=0,\ldots,M-1\) of them supplies every final \(Z\bmod M\), keeping both address residues modulo \(2M\), the left address modulo \(2^{a+1}\), and \(F(x)\bmod 2^a\) unchanged at the end of the batches. For \(M=1\) the bank is empty.

## The final common address word

In the chart \(H=k+tx+cy+qxy+r[x\text{ odd}]+e[y\text{ odd}]+2hZ\), the integers \(c,t,q\) are even, \(h,k\) odd, and \(r\) even (evaluate at \(x=1,y=0\)). These chart properties and \(\gcd(t,h)=1\) hold at each of the eight spines. [Canonical valuations](RUN-COMPLEXITY.md) · [Eight chart gcds](LOCAL-WORDS.md)

For \(a\ge1\), impose \(F(x)=0\bmod 2^a\). For final odd \(x\), even \(y=S-x\) the invariant then gives
\[
 H\equiv k+r+cS+2hC+h +(t-c+qS-h)x-qx^2\pmod {2^a}.
\]
The constant term is even, the linear coefficient odd and \(q\) even. The difference at \(x',x\) is \(x'-x\) times an odd number. Hence the polynomial permutes the residue ring and has a unique solution for \(T\bmod 2^a\), which is odd. For \(a=0\) prescribe odd \(x\) directly. At the odd modulus prescribe \(x=S\bmod M\), hence \(y=0\bmod M\). The remaining height condition is
\[
                     2hZ\equiv T-k-r-tS\pmod M .
\]
Our hypothesis makes it solvable at the full odd modulus.

Let \(N=2^aM\). Use the inherited single-word controller on the baseline state to prescribe \(x\bmod 2N\) by CRT with the just-selected residues and parity, choosing any lift, and \(F(x)=0\bmod N\). Its common length \(L_N\) is independent of the state and the requested residues. Append this same controller word after every possible choice of batches. The terminal \(x\bmod 2N\) and \(F(x)\bmod 2^a\) are the same as on the baseline. Although the batches can change \(F(x)\bmod M\), the two tapes have the same incoming address residues modulo \(2M\) on every variant, so the fixed final word adds the same \(Z\)-increment modulo \(M\) to each. Select the number of toggled batches to solve the last displayed congruence. The common length in the theorem can be taken to be
\[
       \ell=(a+1)+3(M-1)2^a(2L+a+1)+L_N .
\]
This is an existence bound depending only on the requested precision. It also applies to the usual odd phase modulo \(2M\), and to all powers at the supersingular primes. The actual suffix addresses and physical flips depend on the absolute prefix pair.

**FINITE-EXACT.** On spine 4 the replay tests every odd phase modulo each of \(1,2,4,8,9,15,20,24,25,36,45,60\) on three distinct prefix pairs, including zero on either side. For all 516 cases it constructs the literal complementary words, recomputes both full integer chart readings, and checks their difference, sum of addresses and lower phase. The case \(M=25\) includes an odd prime dividing the chart height to lower valuation than the requested modulus. [Integer replay](complementary_controller.py)

~~~sh
python3 -B complementary_controller.py
~~~
