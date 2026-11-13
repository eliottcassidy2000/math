# The ternary count chart has a real and an 11-adic limit

This follows from the [physical count chart](TERNARY_COUNT_CHART.md), with the two-sink palette and its first site fixed. Let $r=33/2880=11/960$. Normalize an $n$-site count by $2880^{n-2}$ and read its address from the last site. Write $a_\ell\in\{2465,2669,2675\}$ for the two-slot end response at digit $\ell$. The overlap part of the normalized count is exactly
$$
508\sum_{\ell=0}^{n-2} a_\ell r^\ell\sum_{u=0}^{n-\ell-2}r^u .
$$
The other terms are a bounded sum of $O(n)$ entries times $33r^{n-2}$; the adjacent contribution vanishes because $s=0$. Extending the digit word to any infinite word gives uniform convergence, in both the real and 11-adic metrics, to the same series
$$
F(a)=\frac{508}{1-r}\sum_{\ell\ge0}a_\ell r^\ell .
$$
The discarded overlap terms have total degree at least $n-1$ in $r$. The real error is $O(nr^{n-2})$, and the 11-adic error is at most $11^{-(n-1)}$ because $33r^{n-2}$ has that valuation and all local entries are integral.

In the real field the three pieces are disjoint: after translating the digits to $0,204,210$, the tail width is $210r/(1-r)<6$, the smallest digit gap. By the separated self-similar dimension formula ([Hutchinson](https://maths-people.anu.edu.au/~john/Assets/Research%20Papers/fractals_self-similarity.pdf)), the real Hausdorff dimension is $\log3/\log(2880/33)\approx0.246$. Over the 11-adics the same addresses have distance exactly $11^{-\ell}$ if the first differing digit is $\ell$: the leading digit differences and the prefactor are units. The equal-weight cylinder measure gives dimension $\log3/\log11\approx0.458$. Both differ from the exponent $\log3/\log2880$ for the number of distinct unnormalized integers through $X$.

## What happens if a future compiler can add copies?

Up to an affine change the real limit is $K=\{\sum_{\ell\ge0}d_\ell r^\ell:d_\ell\in\{0,34,35\}\}$. The $k$-fold numerical sumset has the same contraction and digit set
$$
D_k=\bigcup_{m=0}^k\{34m,34m+1,\ldots,35m\}.
$$
Each digit sum independently chooses $m$ nonzeros and then the number of 35s. This formula is about a numerical sumset; it is not a construction of a tournament operation that adds counts.

* **Up to 11 equal-scale copies have measure zero.** There are at most $\binom{k+2}2$ distinct digits. Depth-$j$ cylinders have total length at most a constant times $(\binom{k+2}2 r)^j$, and $78r<1$.
* **35 equal-scale copies contain a nondegenerate interval.** At $k=35$, the last three digit blocks, $m=33,34,35$, contain every integer from 1122 through 1225. Using just this subalphabet at every depth gives an interval: the candidate hull has tail width $103r/(1-r)>1$, so its images under the consecutive digits overlap and cover it. The minimum $k$ for any interval is not determined here.
* **Exactly 84 equal-scale copies are necessary and sufficient for the whole hull.** The largest digit gap is 34 (attained between 0 and 34). The complete candidate hull is $[0,35k/(1-r)]$. Its first-level translates cover it iff $35kr/(1-r)\ge34$, equivalent to $k\ge84$. If this fails the first digit gap is already outside the attractor; if it holds, the hull is invariant and hence is the attractor.
* **Seven unequally scaled copies already contain an interval.** In each copy keep only digits 34 and 35, and scale the seven copies by $1,2,4,\ldots,64$. After translation this gives all effective digits $0,\ldots,127$, and $127r/(1-r)>1$. The same invariant-hull argument applies. This is a numerical design target for an unequal-scale compiler, not evidence of a legal physical lift.

So one can identify exactly which gaps survive normalization, and one can make numerical ranges overlap by changing the composition rule. None of these limit intervals gives an integer count construction: a physical lift must implement the weights, preserve the continuation ports, and bound finite-stage integer errors at unit scale.
