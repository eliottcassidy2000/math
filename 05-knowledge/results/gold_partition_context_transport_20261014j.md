# Gold certificates under nonuniform extension contexts

Date: 2026-10-14. Status: **PROVED** finite-poset transport and buffer theorems; **FINITE-EXACT** implementation controls. The gold partition conjecture and the 1/3–2/3 conjecture remain **OPEN**. This packet supplies certificates for newly constructed buffered posets; it does not remove the buffer from the original source.

## 1. Inheritance and exact target

Peczarski's gold condition is the weak inequality $t_0\ge t_1+t_2$, for every possible outcome of two comparisons. Comparable elements and a repeated comparison are allowed. See Peczarski, [*Comments on the Golden Partition Conjecture* (2017), pp. 106–109](https://cdm.ucalgary.ca/article/view/62349), and the [2010 workshop formulation, first page](https://dwest.web.illinois.edu/regs/UL2010.pdf). The two-element antichain forces equality $2=1+1$, so replacing $\ge$ by $>$ changes the conjecture.

The closest external structural mechanism is the equal-fibre decomposition for lexicographic substitution: [Dolores-Cuenca–Guzmán-Sáenz–Kim, arXiv:2410.12494v1, Lemmas 2.3 and 2.5](https://arxiv.org/html/2410.12494v1). Reordering the substituted module in its occupied slots preserves every exterior relation. Their displayed Conjecture 1.2 uses a strict sign; here the equal-fibre argument is used with the weak original inequality, proved directly below. No recent census is needed.

The closest repository warning is [crossroads_poset_20260926_board.md](crossroads_poset_20260926_board.md), with the [associated barrier note](crossroads_poset_20260926_barrier.md): conditioning or deleting vertices can change the uniform linear-extension law into a nonuniform law on a smaller poset. The retained sidecar here is the exact number of extensions over each smaller ordering. The live concepts are a fixed comparison certificate, its four deficits, extension fibres, a marked order ideal, and a buffer chain.

Adaptive and fixed second comparisons are equivalent for this weak target. If the first comparison has branches of sizes $m\ge n$, the light branch satisfies $n+t_2\le2n\le m+n$ for every second comparison. Therefore one may choose the successful heavy-branch comparison in advance. This does not justify choosing the most balanced first comparison: that is a separate stronger assertion.

## 2. Four deficits and a sharp robustness radius

Fix a finite poset $Q$ and two comparisons $e,f$. For $i,j\in\{0,1\}$ let $C_{ij}$ be the set of its extensions having those outcomes, $N_{ij}=|C_{ij}|$, and $N_i=N_{i0}+N_{i1}$. The source total is $N_0+N_1$; a leaf on row $i$ satisfies gold exactly when

$$\Delta_{ij}=N_{1-i}-N_{ij}\ge0. \tag{1}$$

Thus a fixed certificate is four linear inequalities. Impossible cells have count zero and impose no additional obstruction. Define

$$\rho=\min_{N_{ij}>0}\frac{N_{1-i}}{N_{ij}}. \tag{2}$$

The certificate is gold iff $\rho\ge1$ and is strict on every actual outcome iff $\rho>1$.

Now let $Q$ be an induced subposet of a larger $P$. Every linear extension $\sigma$ of $Q$ extends to $P$: adding its total order to $P$ creates no cycle, since a contrary path between two core vertices would already be a relation of the induced $Q$. Let $w(\sigma)>0$ count this fibre. Restriction preserves each marked comparison, so the same certificate is valid in $P$ iff

$$\sum_{\sigma\in C_{1-i,0}\cup C_{1-i,1}}w(\sigma)
 \;\ge\;\sum_{\sigma\in C_{ij}}w(\sigma)
 \quad(i,j\in\{0,1\}). \tag{3}$$

This is an iff with the actual weights, not an independence approximation. If $\ell\le w(\sigma)\le u$, a sufficient condition is

$$u/\ell\le\rho. \tag{4}$$

It is sharp for unrestricted positive weights in that box: when (4) fails, put weight $u$ on a cell attaining (2) and $\ell$ on its opposite first row. These sets are disjoint, and their inequality fails. This sharpness does **not** claim that every artificial box assignment is realizable by a poset. Equal module fibres give $u/\ell=1$ and recover the inherited substitution mechanism.

## 3. Exact one-vertex context compiler

Let $Q$ have $n$ vertices, let $D$ be an order ideal of $Q$, and fix a gold certificate. Construct $P_m$ by adding a chain $c_1<\cdots<c_m$ above all of $Q$, and a vertex $z$ above exactly $D$. The vertex $z$ is incomparable with $Q\setminus D$ and the chain. All labels and the original comparison pair are retained.

For $\sigma\in L(Q)$, let $d(\sigma)$ be the position of its last member of $D$, or zero if $D$ is empty. Then

$$w_m(\sigma)=m+\ell(\sigma),\qquad
\ell(\sigma)=n-d(\sigma)+1. \tag{5}$$

Indeed, after fixing $\sigma$, the core followed by the new chain is a single ordered sequence. The allowed gaps for $z$ are precisely those after the last member of $D$. This is a bijection, including the two exterior gaps when permitted.

Put $W_{ij}=\sum_{\sigma\in C_{ij}}\ell(\sigma)$ and $\Gamma_{ij}=W_{1-i,0}+W_{1-i,1}-W_{ij}$. Equation (3) becomes exactly

$$m\Delta_{ij}+\Gamma_{ij}\ge0. \tag{6}$$

Since the supplied certificate has $\Delta_{ij}\ge0$, there are only two possibilities:

* If $\Delta_{ij}=0>\Gamma_{ij}$ for some cell, this certificate fails for every $m\ge0$.
* Otherwise its valid parameters are exactly all integers $m\ge m_*$, where

$$m_* = \max\left(0,\max_{\Delta_{ij}>0}
\left\lceil-\Gamma_{ij}/\Delta_{ij}\right\rceil\right). \tag{7}$$

This is a complete, terminating all-parameter classifier for a supplied marked source, ideal, and certificate. It makes no assertion about other comparison choices in a rejected source.

**Strict positive example.** Let $Q$ be the antichain $\{a,b,c\}$, compare $(a,b)$ then $(a,c)$, and take $D=\{a\}$. In row-major outcome order,

$$N=(2,1,1,2),\quad\Delta=(1,2,2,1),\quad
\Gamma=(-2,2,6,6),\quad\rho=3/2.$$

The certificate therefore works exactly for $m\ge2$. The coarse bound $(m+3)/(m+1)\le3/2$ only gives $m\ge3$. Retaining the cell weights improves the conclusion.

**Zero-slack hostile.** Let $Q=\{a,b\}$ be an antichain and compare $(a,b)$ twice, with $D=\{a\}$. The two extension weights are $m+2$ and $m+1$. The heavy actual leaf has deficit $-1$ for every $m$. The normalized law tends to the uniform law, but the certificate never becomes valid. This refutes an induction step based only on approximate uniformity and a weak core certificate. It is not a counterexample to gold: the enlarged source can use different comparisons.

## 4. Strict certificates survive arbitrary upper contexts

Let $Q$ be an order ideal of a finite $P$, with $|Q|=n$ and $|P\setminus Q|=r$. Write $R=P\setminus Q$ with its induced order. Add a chain of $m$ new vertices above every vertex of $Q$ and incomparable with $R$. Denote the resulting poset by $P_m$. For each $\sigma\in L(Q)$,

$$e(R)\binom{m+r}{r}\le w_m(\sigma)
 \le e(R)\binom{m+n+r}{r}. \tag{8}$$

For the lower bound, place the whole core first, then interleave an extension of $R$ with the new chain. The order-ideal premise ensures that every such order is legal. For the upper bound, erase the cross-relations and interleave an extension of $R$ with the fixed ordered core followed by the chain. Every actual extension has exactly one such presentation.

Consequently, a strict core certificate transports whenever

$$\frac{\binom{m+n+r}{r}}{\binom{m+r}{r}}
=\prod_{j=1}^{r}\left(1+\frac{n}{m+j}\right)\le\rho. \tag{9}$$

For $r>0$ this ratio decreases to one, so a finite sufficient buffer always exists. An exact doubling/binary search finds the least integer satisfying (9); $r=0$ needs no buffer. This applies to arbitrary cross-relations from the core into the exterior, not only a module. The weak hostile in §3 proves that the strict premise cannot simply be dropped for the same-certificate conclusion.

There is also a finite exact refinement. For each $\tau\in L(P)$ restricting to $\sigma$, let $h(\tau)$ be the number of exterior vertices after the last core vertex. Removing the new chain from an extension of $P_m$ gives a unique $\tau$. Reinserting it among those $h$ exterior vertices gives

$$w_m(\sigma)=\sum_{\substack{\tau\in L(P)\\\tau|_Q=\sigma}}
\binom{m+h(\tau)}{h(\tau)}. \tag{10}$$

Thus the marked fibres are exact degree-$r$ polynomials in a binomial basis. There are exactly $e(R)$ summands with $h=r$, independently of $\sigma$; every fibre has leading coefficient $e(R)/r!$. A cell deficit has leading coefficient $e(R)\Delta_{ij}/r!$. Strict slack is the positive leading term that permits the buffer argument. When the leading term vanishes, a lower coefficient can remain negative forever.

Equations (9) and (10) give a complete all-$m$ algorithm for any supplied strict certificate and upper context: find the finite bound from (9), evaluate (10) at every smaller nonnegative integer, and retain every exceptional invalid parameter. Unlike the one-vertex case, monotonicity of the exact deficits is not assumed. The implementation enumerates the original source only, with a declared nine-vertex literal cap; it never enumerates arbitrarily long buffered sources.

**Nonmodule example.** Take the same antichain $a,b,c$ and add $r_1<r_2<r_3$ with $a<r_1$, $b<r_2$, and $c<r_3$. The core is not a module: $r_1$ sees $a$ but neither $b$ nor $c$. Formula (9) gives the sufficient bound $m=19$. Formula (10) gives the exact critical deficit

$$\binom{m+3}{3}-3\binom{m+2}{2}-(m+1)
=\frac{(m+1)(m^2-4m-18)}6. \tag{11}$$

The other three deficits have nonnegative binomial coefficients. Since the positive root of $m^2-4m-18$ lies between six and seven, this certificate works exactly for $m\ge7$. It fails on the original six-vertex source at $m=0$.

## 5. What the proof attempt achieves and still needs

The source is a labelled poset plus two comparison pairs. The map is restriction of full extensions to an induced core. It preserves comparison outcomes and the exact weighted inequality (3). It destroys uniformity unless the fibre counts are retained. A four-cell weighted table suffices for this particular certificate; it does not encode every possible alternative comparison. The cheapest hostile is the two-element core with one asymmetric extra upper vertex.

The universal remaining obligation is to find, in every nonchain source, comparisons satisfying the four actual source inequalities. A deletion-based proof would need to select a core and a certificate whose genuine context weights obey (3) **at the original source**, or prove a different invariant that yields such a certificate. Gold on the smaller core alone does not establish this. Adding a buffer proves a new family's theorem but changes the source. Removing that buffer is exactly the unproved step, and the zero-slack example shows why limiting uniformity cannot replace it.

No claim is made that these elementary transport lemmas are a new proof of any known poset class, that all weights are poset-realizable, or that a failed fixed certificate is a counterexample to the conjecture.

## 6. Exact implementation and reproduction

Owned files are [the script](../../04-computation/experiments/gold_partition_context_transport_20261014j.py) and [the JSON receipt](gold_partition_context_transport_20261014j.json). `compile_gold` authenticates four cell counts by ideal dynamic programming; `insertion_plan` retains the labelled source, ideal, slopes and intercepts; `accepts` recomputes the source-bound plan before consumption. `context_buffer` checks the order-ideal premise; `context_schedule` retains the full binomial fibre histograms. Exact integer/type checks precede cached counting and dataclass equality.

The declared finite universe is all 49 naturally labelled posets of sizes two through four; 356 fixed gold certificates; 3,514 ideal-insertion plans; 14,152 independently counted buffered instances; all 729 assignments of weights in $\{1,2,3\}$ to the six antichain extensions; and the explicit nonmodule context. Malformed counts, bool/float aliases, incorrect source labels, nonideals, altered deficit records and omitted fibre weights are rejected. This is a transport-identity audit, not a new conjecture census.

Reproduce from the repository root:

```text
python -B -X utf8 04-computation/experiments/gold_partition_context_transport_20261014j.py
python -B -O -X utf8 04-computation/experiments/gold_partition_context_transport_20261014j.py
```

Normal and optimized runs agree exactly at **29,251 checks**. The saved JSON,
normalized to LF, has SHA256
`406275e2f9ffdb3c64a3d6fe1e432d74378114f95c81596f047835759c637034`.
Independent peer proof/API audit and optimized-versus-saved replay: **PASS**.
No global conjecture claim depends on the finite controls.
