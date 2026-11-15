# Stability below the first Busch ratio change

This note follows the [sharp unrestricted near-floor theorem](../chanpak_near_floor_crossing_20261115/GLOBAL.md). It proves an all-order description of where a small excess over the Busch floor can live. It also gives a local amplification test for an imperfect seam between two cycle factors. Neither result constructs a tournament at every prescribed excess.

Let \(H(T)\) be the directed Hamiltonian-path count and \(C(T)\) the directed Hamiltonian-cycle count up to cyclic rotation. For strong tournaments, write \(f(n)\) for [Busch's minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3). Put \(F=f(n)\) and \(K=H(T)-F\). A decomposition below means *any* recursive use of Busch's strong one-vertex splits, stopping at nearly transitive factors; a split satisfies \(|U|=|A|+|B|-1\). Nodes are tree occurrences, so disjoint branches need not be vertex-disjoint in the tournament.

## 1. No residue-changing split below \(27/25\)

**Theorem.** If \(K<2F/25\), every such decomposition has all leaves of orders 3, 4 or 5 and satisfies \(f(|U|)=f(|A|)f(|B|)\) at every internal node. In particular every order-four leaf contributes 5 and the entire tree has at most two order-three leaves or one order-five leaf. More precisely the exponent of 3 in \(F\), taking values \(0,1,2\), is exactly the count of order-three leaves plus twice the count of order-five leaves.

For a split, the ratio \(f(|A|)f(|B|)/f(|U|)\) is one of \(1,27/25\), by \(f(3)=3,\ f(4)=5,\ f(5)=9,\ f(m+3)=5f(m)\). The ratio of \(H(U)\) to \(f(|U|)\) is at least the product of this ratio and both ratios of the children, since every pair of factor Hamiltonian paths can be merged. Iterating down the tree gives the product of all split ratios and leaf ratios. The nearly transitive tournament of order \(m\) has \(2^{m-2}+1\) paths. For \(m\ge6\) its ratio to \(f(m)\) is at least \(17/15>27/25\): check \(m=6,7,8\), and on advancing three, the numerator grows by a factor greater than 5. All other factors in the product are at least 1. This proves the claim. The leaf budgets follow by unique factorization applied to the product of all leaf floors.

**Sharp threshold for this assertion.** A chain with an order-three and an order-five nearly transitive link has \(H=27\) at order 7, with floor 25, and gives a split of ratio \(27/25\). Attach any number of order-four links; both values multiply by 5 at every step. Thus \(K=2F/25\) can already support a residue-changing split at arbitrarily large order. The path product for a chain is [Moon–Yang, §2](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf).

## 2. Weighted defect and single-branch localization

For a node \(U\), write \(f_U=f(|U|)\), \(E_U=H(U)-f_U\) and \(w_U=F/f_U\), the product of all leaf floors outside its subtree. For an internal node define the local path surplus \(L_U=H(U)-H(A)H(B)\). Each \(L_U\) is an even nonnegative integer. Throughout the band of §1, subtracting the floor product gives exactly
\[
 E_U=L_U+f_B E_A+f_A E_B+E_AE_B .
\]
After weighting and telescoping all internal contributions, the zero leaf excesses give the exact identity
\[
 K=\sum_{U\text{ internal}} w_U L_U+
   \sum_{U\text{ internal}} w_U E_AE_B .                      \tag{1}
\]
Thus any nonfloor node has \(w_U\le K/2\), and any positive local surplus deep down must pay the product of its sibling floors. If its depth is \(d\), the sibling cost is at least \(3^d\). If \(e=v_3(F)\), the sharper bound is \(3^{\min(d,e)}5^{d-\min(d,e)}\), since at most \(e\) siblings can have a factor of 3, and all other siblings have floor at least 5. When \(n\equiv1\pmod3\), this is simply \(5^d\).

If \(P,Q\) are incomparable nonfloor nodes, their subtrees have disjoint **leaf occurrences**, so \(w_Pw_Q/F=F/(f_P f_Q)\) is an integer at least 1. Preserve the two path counts \(H(P),H(Q)\), replace the other leaf counts by their floors, and propagate the product lower bound up the tree. This gives
\[
 K\ \ge\ w_P E_P+w_Q E_Q+\frac{w_Pw_Q}{F} E_P E_Q
 \ \ge\ 2(w_P+w_Q)+4\frac{w_Pw_Q}{F}
 \ \ge\ 4\sqrt{F}+4 .
\]
Consequently if also \(K<4\sqrt F+4\), **all nonfloor nodes are on one branch**, and (1) simplifies to \(K=\sum w_U L_U\) on that branch. Every off-branch factor is a floor tournament, hence a special chain by [Moon–Yang's characterization, §4](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf). In particular \(K=2,4\) in the band puts the entire excess at the root: any proper node already has weight at least 3. The two factors of any top split in those cases are exact special chains, without needing to choose the split cleverly.

## 3. Exact cycle accounting on the branch

Any two factor cycles yield a good cycle of the union: put \(v\) first in each factor cycle, merge the remaining directed lists with Busch's two-list merge, and prepend \(v\) and close back to \(v\). The first/last vertices are respectively a first/last vertex of a factor cycle. Projection distinguishes different pairs. So the number \(J_U=C(U)-C(A)C(B)\) of additional cycles, counted relative to one per factor pair, is nonnegative. The [local surplus inequality](../chanpak_near_floor_crossing_20261115/GLOBAL.md) gives \(L_U\ge3J_U\) (rounded up to the next even integer). Combining with (1), \(K\ge 3\sum_U w_U J_U\). On the single branch all off-branch factors have one cycle, so \(C(T)-1=\sum_U J_U\) there, exactly. If \(K<12\), no proper node can have more than one cycle: a two-cycle strong node has excess at least 4 by the sharp bound and parity, and its weight is at least 3.

A useful additional exclusion at the top seam: each Hamiltonian cycle whose cyclic restriction to a factor of order \(m\) is invalid contributes at least \(m-1\) broken rotations (the cuts outside one invalid projected gap). For a strong triangle an invalid projection instead makes every rotation broken. So for \(K=2\) all cycles are good at the top, and there is exactly one. For \(K=4\) with both factors of order at least six, all top cycles are good, with at most two. Merely knowing the factors are large does not rule out nonclosing broken paths or repeated path shuffles.

## 4. Prefix collisions can persist across a whole row of cuts

Here is a quantitative constraint on the top gluing when its available local surplus is small. Use any factor cycles \(v,a_1,\ldots,a_p,v\) and \(v,b_1,\ldots,b_q,v\). Suppose \(j\ge2\), \(b_j\to a_i\), \(a_i\to b_1\), and, if \(i>1\), \(a_{i-1}\to b_1\). Then
\[
 L_U\ \ge\
 1+(q-j)\mathbf1_{j<q,\ a_i\to b_{j+1}}
  +(p-i)\mathbf1_{i<p,\ b_j\to a_{i+1}} .                      \tag{2}
\]
This can always be rounded up to the next even integer. The first two shuffles of prefixes with lengths \((i,j)\) are \(a_1\cdots a_i b_1\cdots b_j\) and \(a_1\cdots a_{i-1}b_1\cdots b_j a_i\). Pair either prefix with the same merge of the two complementary tails before \(v\). This is the basic duplicate-witness mechanism of [Moon–Yang, §3](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf). The quantitative extension is that if both shuffles can append \(b_{j+1}\), they can then append as long a consecutive \(b\)-tail as desired; each yields two paths in a separate fibre \((i,j')\), \(j'>j\). The same applies in the \(a\)-direction if both can append \(a_{i+1}\). All these fibres are distinct, since the pair of rotations records separately the two prefix lengths. This proves (2). Converse gives the suffix version.

For budgets 2 and 4, this forces many interior reversals either near the end in the relevant factor or screened by opposite orientations immediately after the collision. It gives a local directional test, not a global bound on the number of reversed cross arcs; the screening case is not resolved. The [exact replay](verify.js) checks (2) against literal subset path counting on 5,008 cross matrices and 7,030 marked prefix-collision configurations.
