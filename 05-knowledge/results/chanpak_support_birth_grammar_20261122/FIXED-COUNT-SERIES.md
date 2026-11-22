# Fixed Hamilton count: exact transients, symmetry by parity, and limiting padding geometry

**Status: PROVED.** The ordinary and converse rational generating functions themselves are inherited from [Session 24](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md) and the [Session 324 fixed-count fiber](../session324_busch_inverse_and_gap_structure_20261113/FIXED_COUNT_FIBER.md). The extensions here determine the exact onset of their tails, the missing parity leading terms and pole order, the first-witness coefficients, and the geometry in a large fixed-count fiber. They do not decide whether an arbitrary odd count occurs.

## 1. The finite core, with grades kept

Write \(h=H(T)\) for the ordinary directed Hamilton path count, not \(H/|\operatorname{Aut}|\). Every tournament has a unique ordered list of strongly connected components. Uniform forward arcs between components force Hamilton paths to concatenate and give product \(H\). Every isomorphism preserves component order, so it is an ordered word of isomorphism classes, not a multiset. For \(h=1\), use the empty word with arbitrary singleton padding; all nonempty tournaments in this fiber are transitive.

For a strong count \(d>1\), let \(A_d(z)\) enumerate nontrivial strong isomorphism classes with \(H=d\) by order. This is finite: [Busch](https://doi.org/10.37236/1141) gives \(d\ge 5^{(|C|-1)/3}\). Delete all singleton components of a tournament to obtain its ordered nontrivial core \(S\), of grade \(k=k(S)\) components and order \(m=m(S)\). Put \(b_{h,k,m}=\#\{S:H(S)=h,k(S)=k,m(S)=m\}\). These are finite and can be computed by ordered divisor products of the \(A_d\). Retain the empty \(h=1\) core, \(k=m=0\). Then the exact inherited ordinary series is
\[
 F_h(z)=\sum_{k,m}b_{h,k,m}\frac{z^m}{(1-z)^{k+1}},
 \qquad
 N_n(h)=\sum_{\substack{k,m\\m\le n}}b_{h,k,m}\binom{n-m+k}{k}.
\tag{1}
\]
There are exactly \(k+1\) scalar singleton-gap variables. Formula (1) imposes their sum \(n-m\).

For attained \(h\), the largest grade \(r_h\) is the largest number of factors above one in a factorization of \(h\) into attainable counts: any attained factor can be refined into its own nontrivial strong factors. The leading coefficient is \(b_h/r_h!\), with \(b_h=\sum_m b_{h,r_h,m}\), and \(r_h\le\Omega(h)\). If all the prime factors of \(h\) occur in the spectrum, equality holds.

## 2. Exact onset, including transient cancellation

For fixed nonempty grade \(k\), let \(M_k=\max\{m:b_{h,k,m}>0\}\). Then the actual grade-\(k\) count equals its eventual polynomial at **every** \(n\ge M_k-k\), and this onset is optimal. A summand at order \(m\) is already zero as a binomial polynomial for \(m-k\le n<m\). At \(n=M_k-k-1\), however, the maximal-order terms contribute \((-1)^k b_{h,k,M_k}\) to the polynomial and are absent from the actual count. All smaller \(m\) are already correct. Thus within a grade the obstruction cannot cancel.

For the ungraded count let
\[
 P_h(n)=\sum_{k,m} b_{h,k,m}\binom{n-m+k}{k},
\]
now viewing the binomial as a polynomial for every integer argument. The exact error at any integer \(n\ge0\) is
\[
 N_n(h)-P_h(n)
 =\sum_{\substack{k,m\\n<m-k}}
 (-1)^{k+1}b_{h,k,m}\binom{m-n-1}{k}.
\tag{2}
\]
This is the negative-argument binomial identity on precisely the terms not yet correct. So \(\max(m-k)\) is an upper bound, but different grades can cancel. Formula (2) determines the exact onset from the finite core without sampling the infinite tail: take one more than the last nonnegative \(n\) for which it is nonzero, or zero if it is never nonzero. Equivalently, divide the rational function into its polynomial at infinity plus its principal part at \(z=1\); the finite discrepancy is exactly that polynomial at infinity.

For instance the inherited complete inventory at \(h=15\) gives the core terms \((k,m,b)=(1,5,2),(1,6,2),(2,7,2)\). The grade-based upper bound is five, but the errors at orders four and three cancel in (2), while the order-two error is nonzero. Thus the exact ordinary polynomial starts at **three**, even though the first actual witness has order five. These are different notions of onset. Busch bounds every \(m-k\le3\log_5 h\), so the search remains logarithmic in the *number of candidate orders* once the finite inventory is known.

A secondary statistic in the tail is visible: for \(r_h\ge1\), the coefficient of \(n^{r_h-1}\) is
\[
 \frac{\sum_m b_{h,r_h-1,m}
  +\tfrac{r_h+1}{2}b_h-\sum_m m\,b_{h,r_h,m}}
 {(r_h-1)!}.
\tag{3}
\]
It is not given by the leading-dimensional core histogram alone; it mixes the immediately lower grade with the mean size of a top-grade core.

## 3. Parity-resolved self-converse tails and the exact negative pole

Let \(C_d(z)\) enumerate nontrivial strong self-converse centers at count \(d\). Set \(C_1^+(z)=1+z\) for empty or singleton center and \(C_d^+=C_d\) for \(d>1\). A self-converse tournament has a unique ordered half \(U\), a possible center \(C\), and the reflected half \(U^{\rm op}\). Here whole-half converse reverses component order and converses its factors; a choice of ordinary half, not a converse orbit, determines the full class. Thus
\[
 J_h(z)=\sum_{a^2\mid h}F_a(z^2)C^+_{h/a^2}(z).
\tag{4}
\]
The converse-merged series is \((F_h+J_h)/2\).

For order parity \(e\in\{0,1\}\), write \(n=2j+e\) and let \(c\) be a center order of parity \(e\). A half core of order \(m\) and grade \(k\) has \(j-m-(c-e)/2\) singleton padding vertices to distribute. Consequently the coefficient for this parity is again an ordinary positive-binomial mixture with effective shift \(m+(c-e)/2\), its own exact grade onset, and its exact discrepancy given by (2) with that shift. There is no cancellation in a fixed grade of a fixed parity. In particular, let
\[
 t_{h,e}=\max\{r_a: a^2\mid h,\ a\text{ attained},\
 C^+_{h/a^2}\text{ has an order of parity }e\},
\]
with absent parity if the set is empty, and define
\[
 b_{h,e}^{\rm sc}
 =\sum_{\substack{a^2\mid h\\r_a=t_{h,e}}}
 b_a\sum_{c\equiv e\ (2)}[z^c]C^+_{h/a^2}(z).
\]
For this parity,
\[
 S_n(h)\sim \frac{b_{h,e}^{\rm sc}}{2^{t_{h,e}}t_{h,e}!}\,n^{t_{h,e}},
 \qquad 2t_{h,e}\le r_h.
\tag{5}
\]
An absent parity has identically zero coefficient, not just zero leading coefficient. Formula (5) gives both part of the missing symmetry coefficient and its actual parity dependence; the degree bound with a single pooled \(t_h\) may hide an entire absent parity.

Let \(P_0(n),P_1(n)\) be the two eventual parity polynomials in the **same order variable** \(n\), including zero for an absent parity. Then the actual pole order of \(J_h(z)\) at \(-1\) is
\[
 0\ \text{if }P_0=P_1,\qquad
 1+\deg(P_0-P_1)\ \text{otherwise}.
\tag{6}
\]
Indeed the infinite tail is \((P_0+P_1)/2+(-1)^n(P_0-P_1)/2\), and finite errors cannot alter a pole. The same negative pole holds for the merged series because the ordinary series has none there. This specifies exactly when the two parity leading terms cancel and, if they do, which next degree controls the oscillation; no unproved asymptotic for strong centers is involved.

## 4. Reconstructing inventories and first witnesses

Use divisor convolution \((f*g)_h=\sum_{d\mid h}f_d g_{h/d}\), with values in rational functions of \(z\); all sums are finite. The ordinary word identity gives
\[
 A=(1-z)\delta-F^{-1},
\]
where \(A_1=0\) and \(F^{-1}\) is the divisor-convolution inverse. This inverse exists because \(F_1=(1-z)^{-1}\). Thus the full ordinary count series uniquely and recursively determine the strong ordinary inventory, not just its degree or total. For centers define \(D_h(z)=F_a(z^2)\) when \(h=a^2\), and zero otherwise. Then \(C^+=D^{-1}*J\). This gives an exact coefficient inverse from scalar census series, provided the full series are retained; it is not a claim about reconstructing individual tournaments or recovering all strong sizes from a finite leading polynomial.

There is also a minimum-order version that avoids full unlabelled inventory multiplication. Let \(a_d\) be the least order of a nontrivial strong witness of count \(d>1\), with \(a_d=\infty\) if absent, and \(c_d\) its number of strong isomorphism classes at that order. Set \(\nu(1)=0,f(1)=1\) for the auxiliary empty witness. The least ordinary order and its number of ordinary first witnesses follow recursively:
\[
 \nu(h)=\min_{\substack{d\mid h\\d>1}}\bigl(a_d+\nu(h/d)\bigr),\qquad
 f(h)=\sum_{\substack{d\mid h,\ d>1\\a_d+\nu(h/d)=\nu(h)}}c_d f(h/d).
\tag{7}
\]
An unattained count has \(\nu=\infty,f=0\). At minimum order there can be no singleton padding, and the first nontrivial strong component is unique; replacing either factor by a shorter witness would shorten the whole. This proves both formulas without rearrangement overcount.

For defined self-converse parity \(e\), let \(a^{\rm sc}_{d,e}\), \(c^{\rm sc}_{d,e}\) be the least order and number of centers of parity \(e\) of count \(d\), allowing the empty and singleton centers for \(d=1\). Then
\[
 \nu_e^{\rm sc}(h)=\min_{a^2\mid h}\bigl(2\nu(a)+a^{\rm sc}_{h/a^2,e}\bigr),\qquad
 f_e^{\rm sc}(h)=
 \sum_{\substack{a^2\mid h\\2\nu(a)+a^{\rm sc}_{h/a^2,e}=\nu_e^{\rm sc}(h)}}
 f(a)c^{\rm sc}_{h/a^2,e}.
\tag{8}
\]
Again ignore absent summands. The empty convention at \(h=1,n=0\) can be removed if only nonempty tournaments are requested. For \(h>1\), the least converse-merged order is \(\nu(h)\) and its number of classes is \((f(h)+S_{\nu(h)}(h))/2\). The overall finite strong inventory bounds come from Busch; these recursions are not a fast recognition theorem in the bit length of an arbitrary \(h\).

## 5. Random fixed-count classes: the actual geometry of the neutral dimensions

Choose a uniform ordinary class of order \(n\) with fixed attained \(h\). As \(n\to\infty\), its nontrivial ordered core converges to the uniform law on the \(b_h\) top-grade cores, with error \(O_h(1/n)\) in total variation. This is stronger than merely having top grade with high probability: each such named core contributes \(\binom{n-m+r_h}{r_h}=n^{r_h}/r_h!+O_h(n^{r_h-1})\), independent of its size at leading order.

Conditional on any fixed core \(S\) with \(k\) components and \(M=n-|S|\) padding vertices, the gap vector is **exactly uniform** on the weak compositions of \(M\) into \(k+1\) slots. Equivalently, the stars-and-bars positions
\[
 t_i=i+\sum_{j=0}^{i-1}g_j,\qquad 1\le i\le k,
\]
are a uniform \(k\)-subset of \(\{1,\ldots,M+k\}\). Hence normalized gaps converge to the uniform simplex law \(\operatorname{Dirichlet}(1,\ldots,1)\); the scaled positions of the nontrivial strong components converge to the order statistics of \(k\) independent uniform points. One can couple the discrete cuts to the continuous ones by first discretizing iid uniforms, then conditioning on distinct cells. The collision probability is \(O(k^2/(M+k))\) and cell error \(O(1/M)\), giving \(O_h(1/n)\) for bounded Lipschitz tests. Under the ordinary law these limiting spacings are independent of the limiting top core, since every top core has the same limiting simplex law.

For self-converse classes on a fixed attained parity, the same argument applies to the *half* with the weights in (5): among maximizing choices of a named half-core and named center the limiting law is uniform, and the normalized half-gap vector is uniform on its simplex. The other half is its exact reflected copy, not an independent simplex. Any real central component stays \(O_h(1)\) from the middle vertex position. For converse-merged classes, choose uniformly one orientation from the selected converse orbit; the same ordinary limiting law follows, since the proportion of self-converse classes tends to zero for \(h>1\).

Finally, under uniform **labelled** sampling at fixed \((n,h)\), core weights change to \(1/|\operatorname{Aut}(S)|\) and the conditional gap law is unchanged. Exactly the same core weights arise from uniform **Hamiltonian-path-normalized tiling presentations** at fixed \((n,h)\): a class contributes \(n!/|\operatorname{Aut}|\) labelled graphs and \(h/|\operatorname{Aut}|\) such presentations, so the conditional class distributions coincide already at finite \(n\). Automorphisms preserve SCC order and multiply over nontrivial factors; singleton padding contributes no weight. Consequently the limiting labelled/tiling law uses the same top cores but with normalized reciprocal automorphism weights. The distinction between an ordinary class, a labelled graph, and a normalized presentation is real even where all three see the same growth degree.

## 6. A growing prime-power fiber: uniform logarithmic-size stability

The previous estimates fix \(h\). There is nevertheless a uniform range in which the exact rank law survives as the target grows. Fix any *attained prime* \(p\), and let \(M_p\) be the total number of nontrivial strong isomorphism classes with \(H=p\). These exist: a witness for prime \(p\) has exactly one nontrivial strong component. Let \(W_p=\sum_C1/|\operatorname{Aut}(C)|>0\) over the same classes. There is a constant \(\gamma_p>0\) such that, uniformly for \(1\le a\le\gamma_p\log n\),
\[
 \begin{split}
 N_n(p^a)&=\frac{M_p^a n^a}{a!}\bigl(1+O_p(a^2/n)\bigr),\\
 \frac{\#\{\text{labelled order-}n\text{ tournaments with }H=p^a\}}{n!}
 &=\frac{W_p^a n^a}{a!}\bigl(1+O_p(a^2/n)\bigr).
 \end{split}
\tag{9}
\]
The same weighted formula is the normalized number of path presentations divided by \(p^a\).

To prove it, write a core count as \(p^{e_1}\cdots p^{e_k}\), with \(e_i\ge1\) and \(\sum_i e_i=a\). At top grade \(k=a\), all \(e_i=1\), so every word of \(a\) independent \(p\)-cores appears; the sum of its sizes is \(O_p(a)\). Its number of pads is
\[
 \binom{n-m+a}{a}=\frac{n^a}{a!}(1+O_p(a^2/n))
\]
uniformly over these words. This gives both leading terms, using ordinary weights 1 or labelled weights \(1/\operatorname{Aut}\) on each component.

For deficit \(d=a-k\ge1\), there are \(\binom{a-1}{d}\le a^d\) ordered exponent compositions. The strong inventory of count \(p^e\) is at most \(\exp(C_p e^2)\), since Busch bounds its order by \(O_p(e)\), and the number of labelled orientations up to that order is at most exponential in its square. The same bound holds for reciprocal-automorphism mass. Among \(e_i\ge2\), their sum is at most \(2d\), so their total squared exponent is at most \(4d^2\). Their contribution, relative to the \(M_p^a\) or \(W_p^a\) inventory factor, is therefore at most \(\exp(C'_p d^2)\); in the weighted case the missing prime factors cost at most \(\max(1,W_p^{-2})^d\), absorbed in that exponential. At grade \(k\), the pad count is at most \(n^k/k!\). Relative to the top scale, these inventories and pads cost
\[
 \exp(C'_p d^2)(a^2/n)^d.
\]
The \(d=1\) contribution is \(O_p(a^2/n)\). Choose \(\gamma_p\) so that \(C'_p\gamma_p\le1/4\). Summing \(d\ge2\) then costs at most \(O((a^2/n^{3/4})^2)=O_p(a^2/n)\). This proves (9).

Thus for this growing range the probability that **every nontrivial component has count \(p\)** is \(1-O_p(a^2/n)\), and the ordered top core is within \(O_p(a^2/n)\) in total variation of \(a\) iid strong \(p\)-cores (uniform, or reciprocal-automorphism weighted). This uses the same uniform pad asymptotic to compare each named word with its product-law probability. The sort-and-couple implementation in Section 5 then couples their normalized positions to \(a\) independent uniform order statistics with failure probability \(O_p(a^2/n)\) and coordinate error \(O_p(a/n)\) on success. For \(p=3\) there is just the directed triangle, recovering the exact triangle top core but now over a target range that increases with \(n\). No assertion is made for primes not known to be attained or for exponents beyond this small logarithmic range.

**Finite-prime-bank extension.** In (9), one may fix any finite set \(\mathcal P\) of attained primes and take
\(h=\prod_{p\in\mathcal P}p^{a_p}\), \(a=\sum a_p\le\gamma_{\mathcal P}\log n\).
The ordinary and weighted leading terms are, respectively,
\[
 n^a\prod_{p\in\mathcal P}\frac{M_p^{a_p}}{a_p!},
 \qquad
 n^a\prod_{p\in\mathcal P}\frac{W_p^{a_p}}{a_p!},
\]
with relative error \(O_{\mathcal P}(a^2/n)\). Indeed the top prime-count word has
\(a!/\prod a_p!\) possible prime-label arrangements, all with the same component inventory weights and the same leading pad count. At deficit \(d\), expand each ordered factor into a canonically sorted sequence of its prime factors. The resulting prime word and its \(d\) missing separators recover the original ordered factorization. Thus there are at most \((a!/\prod a_p!)\binom{a-1}{d}\) such ordered factorizations. A composite factor with \(e\) prime factors has logarithm \(O_{\mathcal P}(e)\); the same Busch square-exponent estimate and defect proof apply. At top grade the limiting comparison model first uniformly shuffles the multiset of prime labels, then independently chooses a strong core conditional on each label with the stated weights. The total variation and position coupling bounds remain \(O_{\mathcal P}(a^2/n)\) and \(O_{\mathcal P}(a/n)\), respectively.

**Concrete negative-pole controls.** In the same complete Session 24 inventory, \(H=25\) has one strong self-converse center on seven vertices, together with the half-count-5 option. Thus \(J_{25}=z^7/(1-z^2)+(1+z)z^8/(1-z^2)^2\). Its two parity polynomials have equal degree-one coefficient but constants differing by \(1/2\); the negative pole has order one rather than two. At \(H=33\) the only centers have one on each of orders six and seven, and \(J_{33}=(z^6+z^7)/(1-z^2)=z^6/(1-z)\), so the negative pole cancels entirely. The positivity at the positive pole is no contradiction.

## 7. The seven-slot degree defect as a closed optimization law

The repository lists the indecomposables of the conjectured monoid in THM-5258 and proves several individual \(7\)-rays in Session 24. Their maximum SCC lengths admit a unified closed rule. Suppose every prime other than \(7\) dividing the particular target \(h\) is attained; no claim about all primes is needed. Put \(v=v_7(h)\), and let \(a=\Omega(h)\).

- If \(v=0\), then \(r_h=a\).
- If \(v\ge2\), then \(r_h=a-\lceil v/2\rceil\). For even \(v\), the top cores have only \(49\) blocks in their \(7\)-bearing components. For odd \(v\), they have pairs \(49\) and exactly one exceptional block, either \(343\) or some attained \(7p\) with prime \(p\ne3,7\) dividing \(h\).
- If \(v=1\), write \(h=7u\) with \(7\nmid u\) and put
  \[
    b=\min_{d\mid u:\ 7d\text{ attained}}\Omega(d).
  \]
  If this set is empty, the fiber is absent. Otherwise \(r_h=a-b\), and its top cores have exactly one \(7d\) block with \(d\) attaining this minimum, with all other factors prime. In particular \(h=7\cdot3^e\) with \(e\ge2\) loses two degrees (one \(63\) block), and \(v=1\) with some attained \(7p\) dividing \(h\) loses exactly one.

For proof, a multiplicative block with \(e\) factors of \(7\) and \(f\) other prime factors spends degree deficit \(e+f-1\). When \(e>0\), its deficit is at least \(e/2\): \(e=1\) needs \(f\ge1\) since \(7\) is absent; \(e\ge2\) is immediate. The *only* equality is \(49\); the only excess exactly \(1/2\) is \(343\) or \(7p\), where \(p=3\) is excluded by the gap at \(21\). Any non-\(7\) block spends no deficit only when prime. The known strong realizers of \(49\) and \(343\) attain the lower bound on total deficit for every \(v\ge2\). This also classifies all cases of equality. When \(v=1\), exactly one component carries \(7\), so its other factors set the cost; any minimizing attained \(7d\) must itself be strongly realized, else splitting off a nontrivial component would lower \(\Omega(d)\). This proves the last bullet.

Under the full Chan–Pak conjecture this gives the uniform formula
\[
 r_h=\Omega(h)-\lceil v_7(h)/2\rceil
       -[h=7\cdot3^e,\ e\ge2]
\]
for every attained \(h\). The conditional scope enters only in making all the needed individual primes and \(7p\) available. The preceding cases are unconditional whenever their stated inputs are known.

The leading coefficient can be read from the classified top **multisets**, then restored to ordered words: if a multiset has \(m_d\) blocks of count \(d\) and \(a_d^{\rm tot}=A_d(1)\) strong choices, its coefficient in the top polynomial is
\(\prod_d (a_d^{\rm tot})^{m_d}/m_d!\).
Indeed ordering those labels contributes \(r_h!/\prod m_d!\), exactly cancelling the \(r_h!\) from spacing. Sum over the classified top multisets. Zero choices contribute nothing. The identical formula with reciprocal-automorphism mass gives the labelled/presentation version.
