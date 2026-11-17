# Period 3p: no equally inflated homogeneous skeleton at width p or larger

**Proved:** no fixed tournament skeleton, with any number of marked vertices simultaneously replaced by transitive modules of common width \(m(p)\ge p\), realizes \(Q(p)=4^p+2\cdot2^p+3\) infinitely often. For a single marked slot the conclusion holds at *any* width schedule. Two further prescribed period witnesses, 57 and 69, are exact. The wider problem remains open.

## 1. Unshifted equal width: why more identical modules cannot help

Marked slots of a fixed skeleton must induce a transitive tournament if common width \(m=p\) works infinitely often. Otherwise three marked slots form a directed triangle. Repeat its slot cycle \(p\) times, with independently permuted vertices in each slot, to get \((p!)^3\) Hamiltonian paths in that induced union. Induced monotonicity (or the previous [cut theorem](../figurate_mandelbrot_shuffle_gap_20261117/README.md)) gives the same lower bound in the full tournament. It eventually exceeds \(Q(p)\).

When marked slots are transitive, a path visits each slot at most once between consecutive fixed exterior vertices. Compress each contiguous slot run to one letter. There are finitely many possible compressed traces; a trace with run counts \(j_i\ge1\) has exactly \(\prod_i j_i!S(m,j_i)\) lifts, where \(S(m,j)\) is a Stirling number of the second kind. Indeed the \(m\) vertices of a transitive slot are assigned onto its ordered visits, and within each visit their order is forced. Each trace is a polynomial exponential in \(m\), with leading base \(\prod_i j_i\) and leading coefficient one. Therefore the largest base is the largest such product and its coefficient is a positive integer.

At \(m=p\), matching \(Q\) on an unbounded sequence forces largest base 4, coefficient 1. A trace with four visits in one slot yields at least six distinct traces of the same profile by permuting the three intervening nonempty gaps. Every gap contains an exterior vertex, since a path inside the marked transitive quotient cannot return to a slot. Interchanging gaps preserves all arc constraints at that slot and keeps the other visit counts. It follows that the leading contribution is exactly one trace with two doubled slots. After ruling out four visits in one slot, a three-visit trace contributes positive leading base 3 with no cancelling base 3 from the two-doubled trace. Hence there is no three-visit trace. If \(h\) is the one-visit coefficient and \(s\) the sum of the single-double coefficients, the result is
\[
h+s(2^p-2)+(2^p-2)^2.
\]
The target forces \(h=11,s=6\). Contracting all slots gives the skeleton itself, of count 11.

Since 11 is prime, this skeleton has exactly one nonsingleton strong component, and [Busch's bound](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md) puts it at order at most five. The supplied independent permutation census gives its unique isomorphism class, with outgoing rows \((24,1,3,6,14)\). Both doubled slots must be in this component; a slot outside it is visited only once. Contracting other marked slots to singletons preserves all relevant traces with their visit counts of one. In the five-vertex core, every pair involving vertex 0 permits a trace with three visits in one slot. Each other pair has double-double coefficient 3 or 4, never 1. Padding by neutral strong components of order one cannot change these counts. This contradiction holds for all fixed skeleton sizes.

## 2. Arbitrary common width: matching the leading power is necessary, but not sufficient

For transitive marked slots and unbounded width the same finite expansion has form \(H(m)=cK^m+O((K-1)^m)\) with positive integers \(c,K\) and \(K\ge2\) (otherwise it is constant). On infinitely many exact hits \(H(m)=Q(p)\), necessarily \(m=\Theta(p)\). The logarithmic form \(\log c+m\log K-p\log4\) would be exponentially small in \(p\). By [Matveev's lower bound for linear forms in logarithms](https://www.mathnet.ru/eng/im314), a nonzero such form in fixed logarithms is bounded below by an inverse power of \(p\); hence \(cK^m=4^p\) on all sufficiently large hits. Thus \(K=2^r,c=2^j\) with \(rm+j=2p\). Eliminating the top mode reveals a necessary second mode: the largest remaining base must be \(2^{r/2}\), since the required remainder is \(2\cdot2^{(rm+j)/2}+3\). Consequently \(r\) is even. The integer coefficient of that mode must be \(2^{1+j/2}\), making \(j\) even as well.

If \(m\ge p\), a triple of cyclic marked slots is excluded as above. The relation \(rm+j=2p\), with \(r>0\) even and \(j\ge0\), now forces \(r=2, j=0,m=p\) at every sufficiently large hit. Section 1 rules this out. Thus longer equal modules, including nonlinear schedules with \(m(p)\ge p\), do not repair the construction. The argument does not exclude equal schedules below \(p\), nor unequal growth rates or nonuniform contacts.

For a single slot the same log-form comparison gives \(cK^m=4^p\) regardless of the width schedule. Its \(c\) is divisible by \((K-1)!\) from permuting the gaps, so \(K\ge4\) is ruled out by the factor 3; \(K=3\) is ruled out by its base. At \(K=2\) the remaining term is constant, not \(2\cdot2^p+3\). This is an all-width obstruction for that one-slot architecture.

## 3. The two further physical targets

Using the [exact binary inverse](../prime_packing_inverse_20261116/README.md), the next two prime inputs after 17 both have a normalized rooted source. Address bits are physical births in low-to-high order, with 0 an ordinary sink and 1 a near-sink.

| p | Q(p) | Outgoing masks | Root | Address | Hamiltonian witness for period 3p |
|---:|---:|---|---:|---:|---:|
| 19 | 274878955523 | 2,4,1 | 0 | 95011519192 | 72057594037665789 |
| 23 | 70368760954883 | 2,4,25,19,35,15 | 2 | 2330279290103 | 295147905179348631549 |

The final column joins the literal Fibonacci fan of count \(2^{p-1}-1\). The second row reuses the earlier \(p=13\) source at a different address. Each row is separately checked by literal small-core permutation counting and integer recurrence; they give no all-\(p\) conclusion.

The lesson for the stagger analogy is more precise now. A homogeneous schedule must match **two** exponential modes, not just the leading one, and at width \(p\) the two-mode requirement collapses the skeleton to count 11 with the wrong seam coefficient. Any proposed longer equal-width decoder can therefore be discarded before a large search.
