# Balance versus abundance: five prime families, seven exceptions, and650

Date:2026-10-28. **PROVED elementary all-size classification; FINITE-EXACT independent controls.** No literature-priority claim. This note connects the inherited F=U+S shape theorem with concurrent factor-pair compression work; it does not change the literal three-proper-divisor definition. The exact650 uniqueness is already proved in the incoming overlap/dyadic notes cited in §5. It is recovered here through the broader five-family/seven-exception abundance classification, with an independent proof; no priority is claimed for either statement.

## 1. The result

Use F,U,S as counts of proper nontrivial divisors, prime proper divisors and squarefree proper nontrivial divisors. Write δ=F-U-S and Delta=s(N)-N, keeping these two defects distinct.

The complete abundant set among balanced sources δ=0 is the following **disjoint union**:

\[
\begin{split}
&\{12r:r\ge5\text{ prime}\}\ \cup\ \{20r:r\ge7\text{ prime}\}\\
&\cup\ \{28r:r\ge11\text{ prime}\}\ \cup\ \{18r:r\ge5\text{ prime}\}\\
&\cup\ \{6p^2:p\ge5\text{ prime}\}\\
&\cup\ \{350,490,550,572,650,748,836\}.
\end{split}
\]

**No F=U+S-balanced source is perfect. Every other balanced source is deficient.** Moreover:

\[
\boxed{\delta(N)=0\ \text{and}\ s(N)=N+2\quad\Longleftrightarrow\quad N=650.}
\]

This is a complete theorem, not merely the outcome of a finite search. It concerns the conjunction of two specified predicates, not every excess-two source.

## 2. Reduce to the mixed balanced shape

The inherited [FUS classification](FUS_DIVISOR_FRONTIER.md) says balanced sources are primes, prime cubes, and p²qr for three distinct primes. Prime powers are deficient. We therefore need only classify the abundance of p²qr.

If all primes are odd, the abundancy ratio is maximized by squaring the smallest prime and taking the prime set{3,5,7}. To see the first assertion, the gain from increasing a prime exponent from1 to2 is

\[
\frac{1+1/p+1/p^2}{1+1/p}=1+\frac1{p(p+1)},
\]

which decreases with p. Each abundancy factor also decreases with its prime. Hence

\[
\frac{\sigma(p^2qr)}{p^2qr}\le
\frac{13}{9}\frac65\frac87=\frac{624}{315}<2.
\]

So all odd balanced sources are deficient. This assertion is restricted to the inherited balanced shapes; it is not an odd-perfect-number theorem in general.

## 3. Squared prime2

Write N=4qr with odd primes q<r. Then

\[
\Delta(N)=7q+7r+7-qr=56-(q-7)(r-7).
\]

- q=3 gives12r and excess28+4r, for every prime r>=5.
- q=5 gives20r and excess42+2r, for every prime r>=7.
- q=7 gives28r and constant excess56, for every prime r>=11.
- q=11 gives excess84-4r, positive exactly at r=13,17,19. These are572,748,836 with excess32,16,8. The zero would be at r=21, not prime.
- q>=13 forces r>=17 and(q-7)(r-7)>=60, so the source is deficient.

These cases are exhaustive. There is no perfect source in this shape. The excess is always divisible by4, so no such source can have excess2.

## 4. Squared odd prime, with the even factor to exponent1

Write N=2p²q, p,q odd distinct primes. Put A=p²+p+1 and D=p²-3p-3. Then

\[
\Delta(N)=3A-Dq.
\]

- p=3: Delta=39+3q>0, giving18q for every prime q>=5.
- p=5: Delta=93-7q. Positive at q=3,7,11,13. The q=3 source belongs to6p²; the other three are350,550,650, with excess44,16,2. No integral q gives zero.
- p=7: Delta=171-25q. Positive at q=3,5. The first belongs to6p²; the second is490, with excess46. No integral q gives zero.
- p>=11: at q=3 the excess is12p+12>0, giving6p². At q>=5 it is at most-2p²+18p+18<0, so no other source is abundant or perfect.

Together with the preceding section this proves the complete classification and the absence of perfect balanced sources. Unique factorization shows the five families and seven exceptions are disjoint: the squared prime and the other two primes identify the displayed case.

Only the exceptional650 has excess2. The five infinite families have excess at least a larger positive constant; the seven exceptional excesses were computed explicitly. Conversely

\[
650=2\cdot5^2\cdot13,\quad(F,U,S)=(10,3,7),\quad\sigma(650)=1302.
\]

Thus its count balance and excess two are both exact.

## 5. Why this repairs the20 connection without changing definitions

Concurrent repository work formalizes factor-pair compression as

\[
C(N):\quad\exists a,b\ge1:\ ab=N,\ s(N)=1+a+b.
\]

It proves C(N) iff τ(N)=4 or s(N)=N+2, and its final incoming-results section already gives the exact mixed-balanced650 uniqueness proof. See the [incoming overlap theorem](../aliquot_divisor_overlap_20261028/THEORY.md), which credits the independent [divisor-pair proof](../dyadic_module_spectrum_20261026/ALIQUOT_OVERLAP.md).

For20 the literal proper-divisor list has five terms, but C uses(a,b)=(1,20). For650 the same unit-pair mechanism yields1+1+650=652, while its actual proper-divisor list has eleven terms and satisfies F=U+S. Thus650 is the **unique balanced mixed source** whose compression uses the excess-two branch.

Prime cubes already satisfy both balance and the literal τ=4 branch, so650 is not the only balanced source satisfying C. More precisely, the balanced compressed sources are exactly prime cubes together with650. Balanced primes themselves do not satisfy C. These distinctions prevent one representation from silently replacing another.

The overlap is not a square-sigma source: sigma(650)=1302 is not square. Balanced square-sigma sources such as364 are a different intersection. The two defect coordinates and the square predicate must remain separate.

## 6. Exact counting and verification

Let B(X) count abundant balanced sources at most X. The disjoint classification gives

\[
\begin{split}
B(X)={}&\max(\pi(X/12)-2,0)+\max(\pi(X/20)-3,0)\\
&+\max(\pi(X/28)-4,0)+\max(\pi(X/18)-2,0)\\
&+\max(\pi(\sqrt{X/6})-2,0)+\#\{e\le X:e\in E\},
\end{split}
\]

where E is the seven-element exceptional set. For X>=836 this is the sum of the five displayed prime-counting functions **minus6**. This counts sources, not distinct aliquot targets. It needs no prime-density conjecture.

Run `python3 check_balance_excess.py`, also with `-O`. An additive divisor sieve independently obtains sigma, divisor count and proper-prime/squarefree counts on every source through100000. A prime-based generator produces the predicted abundant set. The lists, zero-excess exclusion, unique excess-two intersection, and counting formula at every integer cutoff agree. The all-size conclusion rests on the exhaustive prime case split above, not the finite cutoff.
