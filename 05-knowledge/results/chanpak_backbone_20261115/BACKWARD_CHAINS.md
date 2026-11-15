# Backbone chains: exact path-count transfer, full dyadic residue coverage, sparse integer support

Status: **proved in this orientation family**. No new exclusion for unrestricted tournaments and no proof of Chan–Pak. This is an unbounded-gap example where a simple positive two-state geometry is not a universal exact-target decoder.

## Orientation and the interval lemma

Order the anchors $a_0<\cdots<a_k$, with $m_i\ge0$ ordinary vertices between $a_{i-1}$ and $a_i$. Orient every pair forward except the $k$ arcs $a_i\to a_{i-1}$. For $k\ge2$, the anchors alone are strongly connected: $a_k\to\cdots\to a_0$, and $a_0\to a_k$. Every ordinary vertex has $a_0\to v\to a_k$, so the full tournament is strong. For $k=1$ it is strong iff $m_1>0$; for $k=0$ it is a singleton.

An ordinal cut anywhere between consecutive anchors has exactly one backward arc. A directed simple cycle that meets both sides of the cut must therefore use that arc. In particular, its least and greatest vertices in the order must be anchors: if the least were ordinary, the cut just after it forces the preceding anchor into the cycle, and the dual argument handles the greatest. A cycle spanning $r\$ consecutive gaps is thus a fixed backwards chain on $r+1\$ anchors followed by an increasing path through an arbitrary subset of the $M\$ ordinary vertices in those gaps. When $r=1$, the empty subset is illegal, because its closing edge is the reversed pair, but this would in any case have even length. For $M>0$, parity requires subset-size congruent to $r\pmod 2$, giving $2^{M-1}$ odd cycles. For $M=0$, there is one odd cycle if $r\$ is even and at least 2, and none otherwise. Distinct odd cycles are disjoint iff their occupied gap intervals are not overlapping or touching, because selected intervals use all respective anchors.

By the odd directed cycle packing formula of [Grinberg–Stanley](https://arxiv.org/abs/2307.05569), the Hamiltonian path count is the sum over binary strings on the gaps. A maximal run of selected gaps containing a positive $m_i$ contributes $2^M=\prod 2^{m_i}$; an all-zero run of even length contributes 2, and an all-zero run of odd length contributes 0. Unselected gaps contribute 1. This proves an exact effective nearest-neighbor model for arbitrary reversed arc span, rather than approximating the intervals by mutually excluding each other's neighboring single cycles. Two consecutive positive gaps have an odd cycle spanning both; omitting it would incorrectly give counts including the forbidden 7 from putative arithmetic recurrences. The inherited constant-gap-zero specialization is the [ribbon tower of THM-4952](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/ternary_ribbon_metagraph_lines_20261108.md). The positive-gap factorization below includes the one- and two-reversal constructions discussed in [About the number of directed paths in tournaments](https://doi.org/10.1016/j.dam.2019.09.005); this note does not claim literature priority for these special cases.

## Transfer and exact internal factorization

Let $B_r=\bigl(\begin{smallmatrix}\alpha_r&\beta_r\\\beta_r&\gamma_r\end{smallmatrix}\bigr)$ be the weighted sum over a run of $r\$ zero gaps between two positive gaps, for the four selected/unselected boundary conditions. A selected zero run touching a selected positive boundary has weight 1 (the ordinary vertices have already supplied the weight of the positive selected block). An internal all-zero run still has weight 2 in even lengths and 0 in odd lengths. Then
$$
b_0=(1,1,1),\quad b_1=(1,2,2),\quad b_2=(3,3,4),\quad b_{r+3}=b_{r+2}+b_{r+1}+b_r,
$$
where $b_r=(\alpha_r,\beta_r,\gamma_r)$. Indeed
$$
\alpha_r=\alpha_{r-1}+2\sum_{j\ge1:\,2j\le r}\alpha_{r-2j-1}\qquad(r\ge1),
$$
with $\alpha_{-1}=\alpha_0=1$; sum classification by the length of a final selected run, and the last unused gap if present. Subtract the same recurrence two indices earlier to get tribonacci at $r\ge3$ $directly at $r=3$$. Classify by the first unselected gap to obtain
$$
\beta_r=1+\sum_{j<r}\alpha_j,\qquad \gamma_r=1+\sum_{j<r}\beta_j.
$$
Taking first differences propagates the tribonacci recurrence through these sums once the $r=3$ values are checked. In particular $0<\alpha_r\le\beta_r\le\gamma_r$ by the initial triples and forward induction. For the initial zero tail $r\$ and first positive gap $m$, put $(x,y)=(\alpha_r,2^m\beta_r)$. Every later positive gap $m$, preceded by $r$ zeros, updates this as
$$
(x,y)\mapsto(\alpha_rx+\beta_ry,\;2^m(\beta_rx+\gamma_ry)).
$$
A final zero tail $t$ gives $H=\alpha_t x+\beta_t y$; a word with no positive gap has $H=\alpha_t$.

At separator $r=0$, the update gives $(x+y)(1,2^m)$. Therefore one can split **exactly multiplicatively** across any adjacent positive gaps, although their two tournaments share an anchor in the original strong orientation. All words with no zeros have $H=\prod_i(1+2^{m_i})$. In particular gaps $m_i=1,2$ give $3^a5^b$. The powers of 5 run through all odd residues $\equiv1\pmod4$ modulo every $2^j$ (since $v_2(5^{2^d}-1)=d+2$ by successive squaring); one extra 3 toggles the mod-4 class. Every odd residue thus has a strong witness for every fixed dyadic modulus, including residue 1 by a positive period of 5.

## Exact spectrum bound and scope

For $s=49/50$, the half-plane partition sum over the counts of **distinct** members of the multiplicative spectrum converges:
$$
\sum_{h\in \mathcal H}h^{-s}<\infty.
$$
In particular $|\mathcal H\cap[1,X]|=O(X^{49/50})$, simultaneously for unbounded order and unbounded reversal span. To prove it, call a word a piece if it has no adjacent positive gaps, so all its internal separators have $r\ge1$. Every chain count is a product of piece counts by the reset, so we can majorize by the Euler product over distinct piece counts $p>1$.

After any positive block $y\ge2x>0$, since $\beta_r\ge\alpha_r$ and $\gamma_r\ge\beta_r$. The change factor for $x+y$ is a weighted mean of $A=\alpha_r+2^m\beta_r$ and $D=\beta_r+2^m\gamma_r\ge A$, with weight at least (2/3) on $D$. Thus it is at least
$$
\lambda_{r,m}=\frac{\alpha_r+2\beta_r+2^m(\beta_r+2\gamma_r)}3.
$$
The contraction $\sum_{r,m\ge1}\lambda_{r,m}^{-49/50}<24/25$ has an exact rational certificate in `replay_backward_chain.py`. In brief, round each of the $23\times19$ terms with $r\le23,m\le19$ **up** to the nearest millionth by comparing integer fiftieth powers. Their sum is below 0.932. For remaining terms use $\lambda_{r,m}\ge2^m\beta_r$, $\beta_r\ge2^{\lfloor r/2\rfloor}$ by tribonacci, and $2^{-49/50}<3/5=z$. The union bound for $r\ge24$ or $m\ge20$ is at most
$$
\frac{2z^{12}}{1-z}\frac{z}{1-z}+\frac2{1-z}\frac{z^{20}}{1-z}.
$$
The first-positive factor $x+y\ge 2^m\beta_r$ and final tail factor $H\ge\alpha_t(x+y)$ contribute finite sums, again by exponential growth. All-zero words alone also contribute a finite sum. Hence the sum of $H^{-s}$ over all piece words is finite, so the Euler product $\prod_{p>1}$1-p^{-s}$^{-1}$ converges. This also bounds the distinct support of arbitrary products, if additional ordinal sums of these tournaments are allowed. It does not bound counts after arbitrary new outlet arcs or mixed contact cuts, which break the single reverse arc hypothesis.

## Constructive lesson

Here arbitrarily long strong graphs add no power to the scalar channel whenever the gap-zero transfer loses rank. The lost coordinate is the distinction between the two boundary occupation states. Touching a variable-width marked block changes both contact weights, not just one term in a continuant. Before using chain-based clocks for prescribed *absolute* targets, one must find a way to transport both coordinates through a reset, or introduce ordinal cuts with more than one reverse arc so that the interval description does not apply. The dyadic coverage result supplies a concrete hostile control for any criterion certifying absolute coverage from all dyadic residues alone.
