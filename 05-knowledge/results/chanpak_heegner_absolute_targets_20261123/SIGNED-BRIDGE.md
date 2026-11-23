# From two slot readings to absolute signed heights

The modularly rich controller does not by itself produce arbitrary absolute heights. Nevertheless its same two readings combine an exact signed scaffold and the dyadic fillers into a target-directed family, including a literal prime bridge.

## A gap-six bridge with both endpoints prime

Take ordinary filler rows $(244,81,210,39,8,150,56,90)$, and reverse $(0,5)$. The exact reads $(D,Q_2)$ change from $(383,5960)$ to $(389,5960)$, as in the earlier slot-neutral switch. Prepend three ordinal singletons so that $Q$ is now $47680$ in either orientation. Put this filler in two matched homogeneous shells, with absolute height reads $128Q+D$ and $7(128Q+D)$. They give, respectively,
$$ 6103423\leftrightarrow6103429,\qquad 42723961\leftrightarrow42724003. $$
Each displayed pair differs on exactly the same filler-internal arc. Both unscaled heights are prime (exact trial division), and all four output tournaments are strong. These are actual absolute counts, not the small residues of giant unknown counts.

The two exterior shells come from the earlier signed-six rail packet; for self-contained reconstruction the core rows and tag masks are below. A shell with core $A$ and marked set $I$ uses a new singleton $b$ and the entire filler $X$, with $b\to X\to A$ and $I\to b\to A\setminus I$. Enumerate ordered splits $A=R\sqcup S$ into directed paths, with $R$ nonempty ending in $I$; let $G$ be their count and $B$ the subcount for which $S$ is empty or starts outside $I$. A complete Hamilton path uses either one or two filler runs. The two-run exterior skeleton has $G$ choices, and the one-run skeleton has $H(A)+G+B$ choices (start at $b$, reach $b$ before the run from a nonempty prefix, or start in the run then reach $b$). Thus the exact functional is $GQ+(H(A)+B-G)D$. The subset-DP replay recovers these coefficients and verifies the singleton-filler quotients are strong. Replacing its singleton by any nonempty filler preserves strength: arcs between parts lift quotient paths, and a quotient cycle through the filler connects distinct vertices within that part.

| role | core rows | $I$ mask | $(G,H(A)+B-G)$ |
|---|---|---:|---:|
| plus ordinary | 8,1,3,6,15,31,63 | 3 | $(128,1)$ |
| plus sevenfold | 2,132,9,35,15,23,63,125 | 38 | $(896,7)$ |
| minus ordinary | 32,5,81,7,75,30,171,63 | 100 | $(1024,-1)$ |
| minus sevenfold | 8,1,3,166,15,23,63,119 | 80 | $(224,-7)$ |

In the minus-sevenfold shell prepend five extra ordinal singletons to the filler, which multiplies its $Q$ by 32 and leaves $D$ fixed. This aligns it with $7(1024Q-D)$.

## Absolute smooth-offset selector

The certified dyadic atoms have $(D,\log_2Q)$ equal to $(19,8),(217,12),(223,12),(599,13)$, with row masks in the replay. The ordinal join multiplies both reads. For all nonnegative atom multiplicities $a,b,c,d$ put
$$D=19^a217^b223^c599^d,\qquad w=8a+12b+12c+13d.$$
Add $e$ singleton positions to make the filler; require at least one atom or singleton. The shell now gives $p=2^K+D$ and $7p$ when $K=w+e+7$, and $p=2^K-D$ and $7p$ when $K=w+e+10$. These exact phrases are proved for every legal exponent, with no primality assumption.

The inverse on a given positive odd $p$ is especially short. Since $Q\ge2D$, in the plus case $2^K<p\le(1+1/256)2^K$, and in the minus case $(1-1/2048)2^K\le p<2^K$. Therefore each sign forces $K$ immediately from the bit length, and their allowable intervals cannot overlap: the plus interval near the lower power and the minus interval near the next power are disjoint. Compute $D=|p-2^K|$. Repeated exact division by the four fixed pairwise coprime bases decodes the atom multiplicities iff the remainder is one; the remaining exponent determines the number of singleton positions. Check nonnegativity and nonempty filler. At most one signed encoding passes this test, so no absolute search or integer factoring oracle is needed.

For example the replay certifies prime cofactors $1048793=2^{20}+217$ and $1048799=2^{20}+223$, and $549755809651=2^{39}-19\cdot223$, with their sevenfold witnesses. These are new exact instances of this chart, not a claim of first appearance in the full tournament spectrum.

## Why the modular controller cannot be an absolute solver on its own

Fix any nonempty seed tournament of height $H_0$. Consider words using only $P$ (two ordinal singleton sources) and $W$ (cyclic sandwich). Removing final $P$'s changes no height, leaving either the seed or a reduced word $P^{a_1}W\cdots P^{a_k}W$. Put $A=\sum a_i$. A tournament of order $n$ has $Q_2\ge(n+1)H$: cut each spanning path at each of its $n+1$ positions including endpoints, counting an injection into two-slot covers. Also $P^a$ multiplies $Q_2$ by $4^a$ and preserves $H$, and $W$ sends $H$ to $H+Q_2$. Just before the $i$-th stage the existing order is at least $2i-1$. Hence at any such reduced word
$$H\ge H_0\,4^A2^k k!.$$
For $H\le X$, writing $L=\log X$, this imposes $A=O(L)$ and $k=O(L/\log L)$. The number of nonnegative $k$-vectors with total at most $CL$ is $\binom{\lfloor CL\rfloor+k}{k}$. Using the usual binomial bound gives, uniformly for this range, $\exp(O(L\log\log L/\log L))=X^{o(1)}$ reduced words, after summing over depths. The same holds for any fixed finite bank of seeds. Thus within this physical grammar absolute heights are extremely sparse despite the proved local reachability; even covering a positive proportion through absolute searches in this one grammar is impossible.

The same proof yields a genuine exhaustive target inverse for this restricted language. At a node with polynomial $F$, if its height is the target return the current word; otherwise its next wrap after $a$ pads contributes successor height $F(1)+4^aF(2)$, monotone in $a$. Enumerate only the finitely many choices not exceeding the target, transform the polynomial by the exact sandwich operator, and recurse. Every reduced admissible history is included, and the bound applies to the visited nodes. This does not address architectures with other physical edits, or arbitrary target-dependent seeds.

Run `python3 absolute-signed-bridge.py` for the independent row, coefficient, exact inverse, and prime checks.
