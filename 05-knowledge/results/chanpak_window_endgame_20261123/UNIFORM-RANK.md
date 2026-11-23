# A uniform contradiction on multiplicative windows, and what windows cannot prove

**Conditional reduction, not Chan–Pak and not an RH/GRH claim.** A quantitative “lower witness somewhere in every power window versus upper bound everywhere” argument can be made exact for tournaments. The needed statistic measures the number of nontrivial strong components at a fixed count, not ordinary density of attained integers. The novelty here is a uniform finite observation scale and a countermodel for replacing that statistic by ordinary windows.

Let \(K(H)\) denote the largest number of nonsingleton strongly connected components in a tournament with \(H\) Hamilton paths; call \(H\) absent when there is no such tournament. Let \(a_H(n)\) count isomorphism classes of \(n\)-vertex tournaments having \(H\) paths. We use the inherited SCC padding formula and the established values 49, 63, 343 and holes 7,21. In an ordinal sum Hamilton counts multiply. Strip singleton components from a fixed ordered list of \(k\) nonsingleton components of total length \(m\). The exact padding multiplicity at order \(n\) is
\[
 \binom{n-m+k}{k}\qquad(n\ge m),
\]
as the \(n-m\) singletons are distributed among \(k+1\) gaps. Distinct lists and positions are distinguishable by unique ordered SCC decomposition.

## 1. Read the grade from one uniformly specified order

For an odd \(H>1\) put \(L=\lceil\log_2(2H)\rceil\) and
\[
 n_H=2^{16L^2},\qquad
 \Lambda(H)=\frac{\log_2 a_H(n_H)}{\log_2 n_H},
\]
with \(\Lambda(H)=-\infty\) for a zero class count. Then for every attainable \(H>1\),
\[
                      |\Lambda(H)-K(H)|<1/4 .
\]
Indeed the [Busch strong-order lower bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) bounds the order \(m_i\) of a component of count \(h_i\) by \(1+3\log_5 h_i\). Thus the whole stripped list has order
\[
 m\le k+3\log_5 H \le
 \left(\frac1{\log_2 3}+\frac3{\log_2 5}\right)\log_2 H <2L,
\]
and \(k\le L\). There are at most \((2L+1)2^{2L^2}\) stripped isomorphism classes, using all labelled tournaments at each permitted order as an upper bound. A single top-grade class contributes at least \((n_H/2)^{K(H)}/K(H)!\) padded copies. All stripped classes together contribute at most \((2L+1)2^{2L^2}(2n_H)^{K(H)}\), using the binomial upper bound and \(K(H)\le L\). Taking logarithms and dividing by \(16L^2\) gives the claimed strict error. Unlike a pointwise polynomial asymptotic, the error bound is simultaneous in \(H\), even when \(H\) runs through a growing window. The order specified is large; this is a finite theoretical criterion, not an efficient census proposal.

## 2. Two exact window contradictions

For a prime \(p\ne7\), if \(p\) is missing, every nontrivial SCC of count a power of \(p\) uses at least two exponents. Hence \(K(p^r)\le\lfloor r/2\rfloor\) whenever \(p^r\) is attained. If \(p\) is attained, \(K(p^r)=r\) by \(r\) ordinal copies. Thus an exact lower witness
\[
                    \Lambda(p^r)>r/2+1/4
\]
at even one \(r\ge1\) implies that the ordinary prime lane at \(p\) is attained; conversely it holds at all \(r\ge1\) when \(p\) is attained. A missing \(p\) gives the opposite upper bound at **every** point \(p^r\). In particular one can phrase the test as finding a witness along this ray in each sufficiently large window \([Y,Y^{3/2}]\), since those windows contain a ray point for large \(Y\). Having both \(p^2,p^3\) does not help this test: if \(p\) is missing but they are attained, their products reach the upper grade \(\lfloor r/2\rfloor\) for all \(r\ge2\).

For a now-attained prime \(p\ne7\), define
\[
 \eta_p=\min\{s\ge1:7p^s\text{ is attained}\},
\]
with \(\eta_p=\infty\) if the ray is empty. In a witness for \(7p^r\) there is exactly one 7-bearing component, consuming at least \(\eta_p\) of the \(p\)-exponents. At the smallest exponent \(\eta_p\), some realization has a single nonsingleton component: any additional one would make the 7-bearing exponent smaller. Ordinally add \(r-\eta_p\) copies of a \(p\)-witness. Therefore for all \(r\ge\eta_p\),
\[
                  K(7p^r)=r+1-\eta_p.
\]
For \(r<\eta_p\) the target is absent. Hence the second lower witness
\[
                     \Lambda(7p^r)>r-1/2
\]
at even one \(r\) implies \(7p\) is attained, and conversely holds at every \(r\ge1\) if \(7p\) is attained. It too admits the exact “somewhere in each power window versus uniform upper” form. But the lower witness assertions are open, not consequences of profinite control, ordinary density or the analogy to an analytic argument.

If these two window witnesses hold for each odd \(p\ne3,7\), the desired spectrum follows by elementary arithmetic: products of attained non-7 primes; for at least two 7-adic factors use 49 and 343; for exactly one 7 and a factor \(p\ne3,7\) use \(7p\); for exactly one 7 and \(3^j\), \(j\ge2\), use 63 and more 3's. The only leftover cases are 7 and21. Ordinary prime 3 is a directed triangle.

## 3. Additional 7-adic powers erase the marker

Let \(m\) be supported on attained odd primes other than 7, and write \(\Omega(m)\) for their total multiplicity. In fact for **every** \(e\ge2\),
\[
                      K(7^e m)=\Omega(m)+\lfloor e/2\rfloor ,
\]
whether or not any \(7p\) count is attained. For an SCC using \(e_i\ge2\) factors of 7, its contribution to the factorization deficit \(\Omega(H)-k\) is at least \(e_i-1\ge e_i/2\). If \(e_i=1\) it must contain another prime (7 is missing), so the deficit is at least one. Sum and round up: deficit at least \(\lceil e/2\rceil\). The lower construction uses 49 for pairs of 7's, one 343 if \(e\) is odd, and separately each factor of \(m\). Thus unlike ordinary prime holes, the sevenfold hole cannot be amplified merely by adding further powers of the missing prime 7. The one-7 ray really is the appropriate test.

More generally at \(H=7m\) the exact deficit is
\[
 K(7m)=1+\Omega(m)-\min_{d\mid m:\ 7d\ {\rm attained}}\Omega(d)
\]
when there is such a divisor. The same minimality argument shows a least carrier \(7d\) has a strong realization after stripping singletons, and the other prime factors can be attained separately. Consequently replacing a particular \(p\)-ray by a window of mixed products can hide a failed \(7p\) lane behind another already successful carrier.

## 4. Why ordinary windows cannot give this contradiction

Abstract multiplicativity plus finite seeds, bounded gaps, and all congruence phases cannot force either primitive lane. Let \(P,Q\) be arbitrary sets of odd primes above any chosen finite verification threshold, excluding 7. Delete from the positive odds exactly
\[
                  \{7,21\}\;\cup\; P\;\cup\;\{7q:q\in Q\}.
\]
The surviving set is multiplicatively closed: an excluded prime has no nontrivial factorization; an excluded \(7q\) or 21 has only a factorization using the excluded factor 7. All finite verified seeds below the threshold can be retained while choosing either omitted lane independently. The set still has relative density one (the possible exclusions are at most the ordinary and sevenfold primes) and has no more than two consecutive missing odd places past 21, because among three consecutive odd places one is a multiple of 3 and is not excluded. In any fixed arithmetic progression containing odd integers, choose an odd prime \(\ell\ne7\) not dividing the modulus and impose divisibility by \(\ell^2\) by CRT; the resulting subprogression is retained past 21, giving bounded gaps even in the prescribed phase. Thus arbitrarily strong ordinary window coverage and all fixed phase tests can coexist with infinitely many primitive holes. The obstruction is against those abstract implications, not against a construction using actual strong-component data.

The two thresholds are themselves exact integer comparisons:
on the ordinary ray compare \(a_H(n_H)\) to \(2^{4L^2(2r+1)}\),
and on the one-seven ray compare it to \(2^{8L^2(2r-1)}\).
No numerical logarithm or rounding decision is required.

One can also evaluate the entire window at the same order. Given an
upper cap \(X\), replace \(L\) throughout the padding bound by
\(L_X=\lceil\log_2(2X)\rceil\) and take \(n_X=2^{16L_X^2}\).
Then the uniform error \(<1/4\) holds simultaneously for every
attainable \(H\le X\). This uses one common order on
\([Y,Y^{3/2}]\).
