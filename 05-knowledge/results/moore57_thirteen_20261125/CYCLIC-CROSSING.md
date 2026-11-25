# A compulsory quarter-turn on crossing differences that the crossing set cannot have

**CONDITIONAL THEOREM.** Suppose a degree-57 Moore graph has an equitable equal split preserved by a free automorphism of order 13. Assume its crossing matrix on the two panels of 125 orbit labels is circulant after relabeling. Write \(S\subset\mathbb Z/125\) for the set of crossing shifts. The earlier [square-reflection crossing theorem](../moore_fixed_order_spectra_20261124/SQUARE-REFLECTION-CROSSING.md) proves that \(S\) contains exactly one lift of each residue modulo 25. Then:

1. The difference multiset \(S-S\), counting ordered differences, is invariant under multiplication by **57**, whose square is \(-1\bmod125\). The set \(57S\) itself is inequivalent to \(S\) under both translations and translated reflections. Thus a candidate must produce a genuinely different set with the same pair differences.
2. The previous coarse internal trace pair \((20,140)\) is impossible. With sides interchangeable, the remaining trace pairs are \((50,110)\), with respectively \((0,30)\) active order-13 orbits, and \((80,80)\), with \((15,15)\). The 100-dimensional fine spectra of the two internal matrices agree, multiplicities included. In the balanced case the two full 125-dimensional matrices are cospectral; in the other case changing panels trades four copies of \(-8\) for four of 7.
3. The shift polynomial \(f(X)=\sum_{s\in S}X^s\) must vanish at **both** \(\alpha\) and \(\alpha^{57}\) for a primitive 125th root in some prime field \(\mathbb F_q\), with \(q\equiv1\pmod{125}\). One may take \(q\le5^{38}\). A candidate whose primitive resultant has no such paired root is excluded.

4. Given the crossing set and both 25-dimensional coarse internal matrices, **two adaptively chosen rows of one 125-dimensional internal matrix suffice to recover both panels in full**. After row 0, the next row is 1 or 2 according to an explicit test. Three fixed consecutive rows also suffice.

The quarter-turn is an automorphism of the *difference multiset*. No graph automorphism of order four is inferred. The crossing rank and trace inputs concern this additional cyclic representation; the [pentagon crossing constraints](SQUARE-CROSSINGS.md) apply to any equitable equal split without those hypotheses.

## The singular-value field and its trace

Write \(C\) for the circulating crossing matrix and \(U,V\) for the integer symmetric internal orbit matrices. On the orthocomplement \(E\) of functions constant on residue classes modulo 25, the eigenvalues of \(CC^t\) are the values \(g_k=|f(\zeta^k)|^2\) for the 100 primitive frequencies of \(\zeta=\exp(2\pi i/125)\). The previous theorem gives \(0<g_k\le225/4\) and \((2U+I)^2=225I-4CC^t\) on \(E\). The same identity holds for \(V\) because \(C\) is normal. Both internal matrices preserve the coarse and fine spaces. If \(t\) is the trace of the chosen side's 25-dimensional coarse matrix and \(c\) counts its active order-13 orbits (internal degree two within the orbit), then the full trace is \(2c\). The two sides together have 30 active orbits; the earlier candidates for \(t\) were 20, 50, 80, 110, 140. Every coarse *diagonal entry* is one of 0, 2, 4, 6.

Use \(F=\mathbb Q(\zeta)\), \(\pi=1-\zeta\), with \(v_\pi(5)=100\). The ring of integers is \(\mathbb Z[\zeta]\) and the prime over 5 is totally ramified. These are the standard prime-power cyclotomic integral-basis facts. [Sharifi, *Algebraic Number Theory*, cyclotomic fields](https://www.math.ucla.edu/~sharifi/notes/algnum-ch03.html) Since \(S\) is a transversal, the coefficient polynomial of degree at most 124 can be written
\[
 f(X)=(1+X+\cdots+X^{24}) h(X),\qquad h(1)=1 .
\]
Indeed \(f(X)-(1+\cdots+X^{24})\) is divisible by \(X^{25}-1\). Modulo \(\pi^{25}\) this gives \(f(\zeta)\equiv\pi^{24}\), by the binomial expansion of \(1-(1-\pi)^{25}\). The conjugate has the same leading coefficient since 24 is even. Thus \(g=g_1\) has valuation 48 and leading coefficient one. The field \(L=\mathbb Q(g)\) is real, and \(m=[F:L]\) divides both 100 and 48 by total ramification. Therefore \(m=2\) or 4. The distinct Galois conjugates of \(g\) have multiplicity \(m\) on \(E\). In particular \(225-4g\ne0\) by valuation.

Put \(R=(2U+I)|_E\). Over the field \(L\), the eigenspace of \(CC^t\) for \(g\) has dimension \(m\), and \(R^2=(225-4g)I\) there. If this scalar is a nonsquare in \(L\), both roots have equal multiplicities, conjugately at every primitive eigenvalue. If it is a square, write \(\beta^2=225-4g\) in \(L\). The two multiplicities have difference \(d\in\{-m,-m+2,\ldots,m\}\), and
\[
                     \operatorname{tr} R=d\operatorname{Tr}_{L/\mathbb Q}\beta.
\]
This follows by conjugating the eigenspace and both signs around all embeddings of \(L\).

Here \(\beta\) is integral, and \(\beta\equiv\pm\pi^{24}\bmod\pi^{25}\) since \(-4\equiv1\bmod5\). Moreover
\[
                 m\operatorname{Tr}_{L/\mathbb Q}\beta
                 \equiv\pm100\pmod{125}.
\]
Indeed the integral ideal \((\pi^{25})\) equals \((1-\zeta^{25})\). Its trace into \(\mathbb Z\) is contained in \(125\mathbb Z\): on multiples by the integral basis powers this is the difference of consecutive shifts by 25 in the Ramanujan sum \(c_{125}(j)\), whose values are 100 at multiples of 125, \(-25\) at the other multiples of 25, and zero otherwise. Expanding \(\pi^{24}\) gives trace 100 because the remaining powers have zero trace.

In either square class we now have
\[
           \operatorname{tr}R=4c-2t+100\equiv0\pmod{25}.
\]
With \(0\le c\le30\) the possibilities are as follows.

| \(t\) | Possible \(c\) | Corresponding \(\operatorname{tr}R\) |
|---:|---:|---:|
|20|10|100|
|50|0, 25|0, 100|
|80|15|0|
|110|5, 30|−100, 0|
|140|20|−100|

A nonzero trace would force \(|d|=m\). For \(m=2\) it is the only nonzero possibility. At \(m=4\), \(\operatorname{Tr}_L\beta/25\equiv\pm1\bmod5\), whereas \(\pm100/25\equiv\pm4\); multiplying by \(d=\pm2\) cannot give this residue. With \(|d|=m\), the real symmetric \(R\) acts as a scalar on each entire primitive eigenspace. After extension by zero on the coarse space it is diagonal in the Fourier basis and hence circulant, as is the fine part of \(U\).

The diagonal of the fine part of \(U\) is therefore constant. It equals 0 for trace \(R=100\) and \(-4/5\) for trace \(R=-100\). The diagonal contribution of the coarse part at any orbit is its cell diagonal divided by five. Since the orbit diagonal belongs to \(\{0,2\}\) and the cell diagonal to \(\{0,2,4,6\}\), the first case would require *all* cell diagonals zero, giving \(t=0\). The second requires every diagonal four, giving \(t=100\). Both violate the table. Hence the trace of \(R\) is zero. In the square case \(\operatorname{Tr}_L\beta\ne0\) by the same trace congruence, so the two signs have equal multiplicities in every primitive eigenspace in all cases. This proves the joint trace and fine-spectral conclusions. The coarse spectrum consists of 32 and copies of 7 and \(-8\); a trace change of 60 exchanges four copies.

## Why the quarter-turn is compulsory on differences

If \(m=2\), each primitive eigenspace is the pair of Fourier directions \(k,-k\). The corresponding two diagonal entries of \(R\) in the unitary Fourier basis are equal: the eigenvectors are conjugate and the entries are real by symmetry. Their sum is zero, so each is zero. In the integer circulant obtained by summing the 125 cyclic conjugates of \(U\), each primitive Fourier eigenvalue would consequently be \(-125/2\). A rational eigenvalue of an integer matrix must be integral. This excludes \(m=2\).

Thus \([F:L]=4\). The Galois group of \(F\) is cyclic of order 100 and the subgroup fixing \(g\) has order four, generated by multiplication by 57 since \(57^2+1=26\cdot125\). All primitive Fourier coefficients of the autocorrelation are fixed by this multiplier. Its other nonconstant coefficients are zero by transversality. Inversion of the 125-point Fourier transform proves that the actual multiplicities in \(S-S\) are invariant under 57.

**FINITE-EXACT exclusion on the shifts themselves.** If \(57S\) is a translate of \(S\), shifting the origin makes \(S\) invariant under 57 because 56 is a unit modulo 125. A translated reflected image gives the same implication: the affine action with linear part \(-57\) has cube with linear part 57. In normalized coordinates the unique representative of residue zero mod 25 must be fixed, hence zero. The other 24 residue classes form six orbits of length four, and for each there are five choices of lift generating the entire orbit. For **every** one of these \(5^6=15625\) normalized crossing sets the largest centered singular value is greater than \(39/5=7.8\), contradicting the necessary upper bound \(15/2\).

For an exact certificate select a centered integer vector \(v\) for each primitive frequency up to the quarter-turn, using rounded cosine values as trial vectors. For each of the 15,625 sets at least one vector satisfies the integer inequality \(5|v^tCv|>39v^tv\). Its Rayleigh quotient bounds the centered operator norm below. The floating-point choice of trial vectors does not enter any acceptance or rejection inequality. The replay checks them all using integer arithmetic. [Quarter-turn replay](quarter_cross.py)

## Recover both fine matrices from two adaptively selected rows on one side

Once \(C\) is given, the primitive directions break into the 25 distinct real four-dimensional eigenspaces of \(CC^t\), indexed by the quarter-turn orbits \(\{k,-k,57k,-57k\}\). Every one of those spaces reduces \(U\), and its trace there is \(-2\) by equal multiplicities of the two signs of \(2U+I\). For a particular eigenspace \(W\) let \(P\) be its real orthogonal projection. A complete physical row at label \(j\) gives \(P U e_j=U P e_j\) by symmetry and reduction. Take \(j=0,1,2\), or any three consecutive labels. The corresponding three vectors \(P e_j\) are linearly independent. In the Fourier basis their coordinates at the four distinct frequencies form the first three rows of a Vandermonde matrix, multiplied by nonzero scales.

Hence those rows determine the symmetric operator \(U|_W\) on a three-dimensional subspace, including its images in all four coordinates. Two symmetric completions can differ only by a multiple of the orthogonal projection onto the remaining line. The prescribed trace forces that multiple to vanish. Thus every primitive block is reconstructed. On the full fine space the earlier crossing identity gives \(V=-I-C^{-1}UC\), since \(C\) is invertible there. The coarse blocks lift as matrices constant on pairs of five-label fibers, with each coefficient divided by five. Together this reconstructs the two 125-by-125 internal matrices. Outputs from arbitrary proposed three rows must still satisfy the square equations, nonnegative integral entries, and the physical order-13 voltage conditions.

In fact the row budget falls to two. On one four-dimensional eigenspace \(W\) continue to write \(R=2U+I\); its square is the known positive scalar \(225-4g\). First **the projected row probe \(P e_0\) cannot be an eigenvector of \(R\)**. By rational Galois conjugacy it suffices to test the block containing the frequency 1. If it were an eigenvector, the Fourier coefficient at frequency 1 of \((2U+I)e_0\) would be the eigenvalue times that of \(e_0\). Equivalently the polynomial
\[
                        h(X)=1+2\sum_{j=0}^{124}U_{j0}X^j
\]
would satisfy \(h(\zeta)^2=225-4g\) (using either Fourier sign). Hence \(h(\zeta)\equiv\pm(1-\zeta)^{24}\bmod\pi^{25}\). The quotient \(\mathcal O_F/(\pi^{25})\) is \(\mathbb F_5[X]/((X-1)^{25})\) under \(X\mapsto\zeta\), since it is generated by \(1-\zeta\) and has cardinality \(5^{25}\). In characteristic five, \((X-1)^{25}=X^{25}-1\) and \((X-1)^{24}=1+X+\cdots+X^{24}\). For every nonzero residue \(j\bmod25\) this gives \(2b_j\equiv\pm1\pmod5\), where \(b_j\) counts the 0th orbit's internal neighbors in that residue class. For either sign every such \(b_j\) is at least two. This would require at least 48 internal neighbors in a row of sum 32.

Consequently \(P e_0, R P e_0\) span a two-plane known after the first row. At least one of \(P e_1,P e_2\) lies outside this plane by the Vandermonde independence of \(P e_0,P e_1,P e_2\). Select the least such label \(j\) on the block containing frequency 1. The selection works on *all* primitive blocks: the projections and the tested dependence relations are Galois conjugate, and the row vectors and matrices before projection are rational. After reading row \(j\) we know \(R\) on \(Z=\operatorname{span}(P e_0,P e_j)\). The equality for \(R^2\) gives its action on \(RZ\) as well. Now \(Z+RZ\) has dimension at least three because \(R P e_0\notin Z\); fill any remaining scalar by the trace as above. This proves the stated two-row decoder without a generic-position assumption.

## A prime-field witness for the lost phase

Let \(\sigma(\zeta)=\zeta^{57}\). If the ideals of \(f(\zeta)\) and \(\sigma f(\zeta)\) coincided, their quotient would be an algebraic unit all of whose complex conjugates have modulus one by the autocorrelation equality. Kronecker's theorem makes it a root of unity, hence \(\pm\zeta^j\). The quotient is one modulo \(\pi\) after canceling its valuation 24 since \(57^{24}=1\bmod5\), excluding the minus sign. Equality at all primitive roots up to a translate gives equality of the two shift sets: the difference polynomial mod \(X^{125}-1\) is a multiple of \(\Phi_{125}(X)=1+X^{25}+\cdots+X^{100}\), with coefficients constant in each residue class modulo 25, while each class has coefficient sum zero by transversality. This contradicts the finite exclusion. Some prime ideal therefore distinguishes the two factors.

At the unique prime over 5 their valuations agree. For \(q\ne5\) the decomposition group is generated by \(q\bmod125\). If its order is even, complex conjugation fixes the prime, and \(f\bar f=\sigma(f\bar f)\) forces equal valuations at the prime and its \(\sigma\)-image. Thus the distinguishing prime has \(q\equiv1\bmod5\) and multiplicative order 1, 5 or 25 modulo 125. At its four \(\sigma\)-conjugate primes the valuations \(e_0,\ldots,e_3\) satisfy \(e_0+e_2=e_1+e_3>0\). Two adjacent primes must therefore divide \(f\), giving a primitive root \(\alpha\) in the residue field with \(f(\alpha)=f(\alpha^{57})=0\).

Order 25 is impossible. Here \(\alpha^{25}\in\mathbb F_q\) and \(\alpha\) has degree 25. Writing \(f(\alpha)=\sum_{i=0}^{24}\alpha^i(\alpha^{25})^{t_i}\) in its power basis yields 25 nonzero coefficients.

At order 5 write \(\beta=\alpha^5\in\mathbb F_q\) of order 25. Separating the five residues modulo 5 gives \(f(\alpha)=\sum_{i=0}^4\alpha^i g_i(\beta)\), where each \(g_i(X)\) has five exponents, one in each residue modulo 5, read modulo 25. Both primitive roots give bases; hence every \(g_i\) vanishes at \(\beta\) and \(\beta^7\). For a primitive complex 25th root \(\xi\), \(g_i(\xi)\) has valuation four at \(1-\xi\) by the same transversal calculation. Parseval gives mean squared absolute value five over its 20 conjugates: on all 25 characters the sum is 125, with 25 from the constant and zero from the four other nonprimitive characters. Consequently its norm in absolute value is at most \(5^{10}\), is nonzero by valuation, and is divisible by \(5^4q^2\), since the two distinct roots mod \(q\) give distinct degree-one primes. Thus \(q\le125\). The only prime in this range congruent to 1 modulo 25 is 101.

Any primitive 25th root modulo 101 is a power of 16 by a unit exponent. Multiplying the five shift exponents by that unit preserves transversality, so one can fix \(\beta=16\). A further translate by a multiple of five fixes the representative of zero to zero. Of the \(5^4\) remaining polynomials ten vanish at 16. Their residues at \(16^7\) are 62, 83, 56, 74, 19, 73, 95, 15, 80, 31. This excludes order five. [Exact modular test](quarter_cross.py)

The remaining multiplicative order is one, giving the claimed prime-field witness. Parseval now gives mean square 25 on the 100 primitive characters at order 125. Thus \(|\operatorname{Norm}_F f|\le5^{100}\). This norm is divisible by \(5^{24}q^2\), giving \(q\le5^{38}\). The filter concerns a proposed transversal, so its own resultant supplies a finite list of primes to test. None of these necessary conditions supplies the missing graph.

Run `python3 -B quarter_cross.py` in this directory for the two finite exclusions. The spectral and norm deductions apply to all sizes of coefficients satisfying the stated transversal and quotient conditions.
