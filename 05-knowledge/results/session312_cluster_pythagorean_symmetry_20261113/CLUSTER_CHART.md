# A cluster chart that remembers the ternary addresses

[Session 311](../session311_recursive_symmetry_20261113/WREATH_ATLAS.md) gives equal-polynomial tournaments with a hierarchy of rotations. This note constructs an exact rank-two cluster chart for the *hierarchy*. It also explains why it does not generally give a chart for the interiors of the tournament modules.

## An integer chart with two readouts

Number the leaves $0,\ldots,N-1$ in root-first lexicographic order, where $N=3^d$. Let $r(i)$ reverse the $d$ ternary digits of $i$. Use vectors $(1,x_i)$ with $x_i=Ni+r(i)$. Their positive minors are $p_{ij}=x_j-x_i$ for $i<j$. They satisfy the type-A Plücker exchange
$$
p_{ik}p_{jl}=p_{ij}p_{kl}+p_{il}p_{jk}\qquad (i<j<k<l).
$$
See [Fomin and Zelevinsky, Cluster Algebras I](https://arxiv.org/abs/math/0104151). This adapts the older [positive-minor chart on Collatz carries](../frieze_collatz_reset2_20261007.md), which preserves a word under a chart flip only after retaining the needed frozen data. The old [tropical Hamiltonian experiments](../tropical_cluster_H.out) instead look for dominant contributions to $H$; the chart here uses an exact discrete valuation, not a dominant-term approximation.

The valuation $v_3(p_{ij})$ equals the length of the common address prefix. In fact the first different root digit becomes the first different 3-adic digit; the term $Ni$ does not contribute at that precision. Moreover the leading unit of $x_j-x_i$ modulo three is the digit of $j$ minus that of $i$ at that split. Interpreting unit $1$ as an arc from $i$ to $j$, and unit $-1$ as the reverse arc, **recovers the cyclic tournament on the addresses themselves**. The valuation gives where two leaves diverge and the leading unit gives the direction of their arc. This statement is about the quotient with singleton leaves, before adding the colored modules.

Exchange tropicalizes *exactly* here. With $a=v_3(p_{ij}), b=v_3(p_{jk}), c=v_3(p_{kl})$, the three outer interval valuations are respectively $\min(a,b),\min(b,c),\min(a,b,c)$ for $(i,k),(j,l),(i,l)$. This follows because all leaves sharing a prefix form one lex interval. The identity
$$
\min(a,b)+\min(b,c)=\min(a+c,b+\min(a,b,c))
$$
follows by expanding the left side as the minimum of $a+b,a+c,2b,b+c$. Thus there is no extra valuation from cancellation in this exchange, including tied cases.

## Recovering the marked stabilizer from minors

For $k>0$ sorted marked addresses $i_1<\cdots<i_k$, put $l_s=v_3(p_{i_s i_{s+1}})$. The number of occupied depth-$j$ nodes is $1+\#\{s:l_s<j\}$: a depth-$j$ prefix starts a new run exactly at those boundaries. Thresholding all $l_s$ at every depth recovers the occupied trie up to child-position names. The missing children are the maximal completely unmarked subtrees, and their full rotation groups are the direct-product factors of the pointwise stabilizer of a monochromatic tower. Thus for $d\ge1$ the exponent is
$$
(3^d-1)/2-1-k(d-1)+\sum_{s=1}^{k-1}l_s .
$$
This gives the actual stabilizer decomposition even when the exponent alone cannot distinguish it. For an arbitrary leaf coloring the missing child positions still indicate where rotations can survive, but the groups of those subtrees then depend on their colors. A chart flip re-expresses the same marks. It is not a physical recoloring or an edge reversal. This distinction matches the exact [marking-budget result](../session311_recursive_symmetry_20261113/MARKING_BUDGETS.md).


The sum of consecutive $l_s$ has an order-free interpretation: it is the maximum spanning-tree weight on the marked addresses with edge weight $w=v_3(p)$. Indeed, at each depth $j\ge1$, a tree has at most $k-a_j$ edges within occupied prefix classes, where $a_j$ is the number of those classes. Its total weight sums these internal-edge counts for $j=1,\dots,d-1$. The path through the lexically sorted marks attains every upper bound, because each prefix class is consecutive. This also gives a greedy rule that does not depend on a chosen digit order: at each step choose a leaf minimizing its maximum $w$ to the current marks. For $0<k<3^d$, if the first $k$ marks have optimal occupancies, put $t=\lfloor\log_3 k\rfloor$. All depth-$t$ classes are occupied and some depth-$(t+1)$ class is empty. The minimizers are exactly the leaves in empty depth-$(t+1)$ classes, and any such choice preserves optimality at every depth. Induction from any first mark proves the rule for arbitrary ties. The exponent drops by exactly $d-1-\max_s w(x,i_s)$ on adding a leaf $x$ to a nonempty set.

## A boundary even for 3-adic orientation charts

Write $\chi(x)$ for the leading 3-adic unit of a nonzero number, identified with $\pm1$. A tournament is the orientation given by $\chi(\det(v_i,v_j))$ for projectively distinct $v_i\in\mathbb Q_3^2$ exactly when some vertex switch of it belongs to the grammar generated from singletons by binary ordinal joins and cyclic ternary joins.

To prove this, choose a generic affine chart with all points finite; the coordinate change can be taken to have determinant one (choose its first row away from finitely many forbidden directions and complete it). Writing each column as $s_i(1,x_i)$ shows that scalar factors only switch vertices. Split the $x_i$ by the first noncommon 3-adic digit, after translating and scaling by a power of three. Two occupied residues give an ordinal join and three give a cyclic join; recurse inside each residue. Conversely choose distinct residues at a split, reverse the order of the two if needed, and translate and scale finite child charts into the prescribed residue classes. Powers of three do not change $\chi$; a factor of $-1$ at selected vertices then realizes any switches.

The order-eight core $B$ has an explicit chart: in the vertex order of the atlas take $x=(12,1,3,9,5,0,2,27)$. Its neutral partner $D$ does not. The induced six-vertex control on $D$ vertices $(0,1,2,3,6,7)$ has outgoing rows $(34,44,41,49,39,0)$, with bit zero the first vertex. For every subset $S$ of size 2, 3 or 4, the two-column sign products $b_{ji}b_{ji'}$ with $i,i'\in S$ are not constant across exterior rows $j$. This is a 50-subset exact check, and invariant under switching since such products change only by the column signs. So no switch has a module of those sizes. But a binary or ternary root on six vertices must have a child of those sizes unless it is binary $1+5$. Splitting the five-child next then supplies a module of size 2, 3 or 4 in the full graph. Contradiction. This is distinct from the old real rank-two [four-vertex sign obstruction](../twoanchor_reset2_friezes_20261007.md) (THM-4602): a real sign obstruction need not obstruct $\chi$.

In particular the full tournament with any $D$ leaf cannot have this 3-adic determinant orientation, since the property is hereditary. The cluster exchange chart nevertheless stays exact on the entire address geometry. This pinpoints which part of the recursion can use a cluster description and which part still needs tournament-specific port data.
