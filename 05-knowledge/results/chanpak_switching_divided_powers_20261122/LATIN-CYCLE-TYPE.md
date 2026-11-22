# Reading a square-sum path from a tournament cycle cover

**PROVED.** Each symmetric Latin trade in the five-square construction has a cycle type that can be recovered from a correlation on labelled tournaments. The marked Hamiltonian cycle alone detects exactly the connected trades. Marking a cycle cover of a specified type also distinguishes all the disconnected cases.

## A cycle-type test with two matching rows

Let \(P,Q\) be edge-disjoint perfect matchings on \(n\) labelled vertices. Their union is a collection of even cycles of lengths at least four. Fix reference edge directions and let \(p_P,p_Q\) be the products of the orientation signs on the two rows of an ordinary tournament. Let \(\lambda\) be a partition of \(n\) into parts at least three, with \(m_j\) parts of size \(j\), and \(z_\lambda=\prod_j j^{m_j}m_j!\).

Choose uniformly a permutation of cycle type \(\lambda\), orient its cycles along that permutation, and orient the other pairs fairly and independently. This gives uniform pairs consisting of a tournament and a marked directed cycle cover of type \(\lambda\). Then

\[
 \mathbb E_\lambda(p_Pp_Q)=
 \begin{cases}
   \epsilon\,2^{\ell(\lambda)}z_\lambda/n!,&
      \text{the cycle lengths of }P\cup Q\text{ are }\lambda,\\
   0,&\text{otherwise}.
 \end{cases}
\]

Here \(\epsilon\) is the product of the reference signs along traversals of the undirected cycles. It is independent of each traversal's direction because every cycle has even length.

**Proof.** An unplanted edge in \(P\cup Q\) has a free fair sign, making the conditional product mean zero. A planted permutation cover has exactly \(n\) distinct undirected edges, as does \(P\cup Q\), so containment forces equality. This happens for exactly \(2^{\ell(\lambda)}\) permutations if the types agree. The total number of permutations in the conjugacy class is \(n!/z_\lambda\), by writing cycles in order, then forgetting their origins and the ordering of cycles of equal length. All fitting orientations give the stated sign. The single-part case is the [Hamiltonian row correlation](../chanpak_matching_parity_faces_20261122/CYCLE-ROWS.md).

## The same permutation in three places

Use the two labelled shores of the cyclic Latin square of even order \(B\). Symbol zero pairs \(x\) on the first shore with \(t(x)=-x\) on the second. After a symmetric trade on symbols \(A,C\), the positions of symbol \(A\) give a second matching \(R\). Reference all inter-shore edges from the first to the second. For the trade parameters
\[
 g=\gcd(B,C-A),\quad M=B/g,\quad D=(C-A)/g,\quad
 d=\gcd(A+2\ell D,M),
\]
where \(\ell\) is the number of exchanged alternating cycles, the first-return formula makes \(t^{-1}R\) a permutation with exactly \(d\) cycles of length \(B/d\). This is the permutation governing the square-sum path. [Trade theorem and five-square construction](../chanpak_square_latin_trades_20261122/PROOF.md)

On two shores, every step of \(t^{-1}R\) becomes two alternating matching edges. Consequently the only nonzero cycle-cover correlation has partition \((2B/d)^d\), denoting \(d\) equal parts. Its value is
\[
             (-1)^B\,\frac{2^d(2B/d)^d\,d!}{(2B)!}.
\]
Identifying the shores gives the original graph, with the two fixed points of \(t\) supplying its path endpoints. The following three events are equivalent for each trade.

1. Its prescribed square-sum edges form one Hamiltonian path.
2. The two symbol matchings in the traded Latin square have one alternating cycle.
3. Their orientation parities have nonzero correlation under a marked Hamiltonian cycle in a tournament on \(2B\) vertices.

In the connected case the correlation is \(2/(2B-1)!\), since \(B\) is even. Under a marked Hamiltonian **path** the two parities are independent for all the trades, by the matching-face theorem. For example the two triples descended arithmetically from 3–4–5 give \(B=705600\). Each of the \(\binom{4200}{2100}\) central choices is detected by the single-cycle statistic. A single exchange from the all-original matching gives \(d=3\); its signal moves to the cover with three cycles of length 470400. The integral statistic on the square graph and the correlation under a tournament measure carry the same predicate. The correlation is an ensemble average and supplies no tournament of a prescribed total path count.

## One Berggren move puts every seed in the construction

Write a primitive triple with Euclid parameters \(m>n>0\), coprime and of opposite parity. Apply the middle Berggren branch \((m,n)\mapsto (m',n')=(2m+n,m)\). For every seed the resulting primitive triple \((u,v,s)=(m'^2-n'^2,2m'n',m'^2+n'^2)\) satisfies the chamber hypothesis of the five-square theorem. In particular at least one third of the addresses at each positive depth supply graphs with \(\exp(cB^{5/8})\) prescribed Hamiltonian paths for a common \(c>0\). This fraction counts tree addresses; different addresses may give the same integer order. [THM-4057, Euclid-parameter Berggren branches, equation CW21](../../../01-canon/theorems/THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)

**Proof.** Put \(x=n'/m'=1/(2+n/m)\). Then \(1/3<x<1/2\), giving \(3/4<y=v/u=2x/(1-x^2)<4/3\). For the two Euclid pairs \((us, v(u+v))\) and \((vs, u(u+v))\), use the respective ratios
\[
 f(y)=y(1+y)/\sqrt{1+y^2},\qquad f(1/y).
\]
The function \(f\) is increasing on positive arguments, and both ratios lie between \(f(3/4)=21/20\) and \(f(4/3)=28/15<1+\sqrt2\). Thus both differences of squares are positive in the indicated larger-to-smaller order and strictly less than twice the product. This puts both primitive triples in the chamber required by the trade theorem. Throughout this interval \(u,v,s,u+v\) are comparable with absolute constants. The resulting \(B=(2uv(u+v)s)^2\) is comparable to \(s^8\), while the number of independent exchanges \(2uv(u+v)s^2\) is comparable to \(s^5\). The central binomial layer yields the claimed uniform exponent. Exactly one third of the words of any positive length end in this branch.

The Berggren tree thus provides inputs for which a single operation places *every* source in a useful range. Its three legal successors come with unique parents. A ternary construction for Chan–Pak would also need inverse moves defined at every required **integer target**, with legal tournament data and strictly decreasing height. The theorem above supplies its inverse test on trades, where it is a gcd, and makes explicit why projecting to the path count alone loses that test.

**FINITE-EXACT.** Run `python3 check_cycle_types.py` from this directory. At order eight it enumerates all permutations, retains admissible planted covers and averages over the remaining free arc signs. The connected and disconnected matching pairs have correlations 1/2520 and zero on the Hamiltonian cover, and zero and 1/315 on the (4,4) cover. A separate trade with B=12, A=1, C=5 gives three eight-cycles on doubled shores. [Reproduction script](check_cycle_types.py)
