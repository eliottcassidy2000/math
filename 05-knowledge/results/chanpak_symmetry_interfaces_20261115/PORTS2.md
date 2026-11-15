# Matched binary outlets with full slot state, and hidden recursive symmetry

**Exact results.** The two physical matched binary return hosts act on the entire reservoir slot polynomial. Their source path counts agree, but their second-slot readings provably differ. Each host preserves the reservoir automorphism group. A six-vertex slot-twin pair, one rigid and the other with an order-three rotation, gives exponential full-profile ambiguity at every order at least eight, and exact symmetry recurrences on a ternary substitution tower. Chan–Pak remains open.

## 1. The single-entry functional reader

For a tournament \(U\), \(P_U(k)\) counts ordered path covers in \(k\) slots including empty slots, for nonnegative integer \(k\). Use \(P_\varnothing=1\). The [Grinberg–Stanley odd-cycle packing formula](https://arxiv.org/abs/2307.05569), with each cycle and each leftover singleton allocated to one slot, gives
\[
P_U(k)=\sum_{\mathcal C}2^{|\mathcal C|} k^{|U|-\sum_{C\in\mathcal C}|C|+|\mathcal C|},
\]
where \(\mathcal C\) consists of disjoint nontrivial odd directed cycles. In particular \(P_U(1)=H(U)\), \(P_U(2)=Q(U)\), and the polynomial is monic of degree \(|U|\).

Take a host \(B\) and replace a marked vertex \(v\) by a nonempty \(W\), preserving every exterior incidence. Suppose \(v\) has a single in-neighbor or a single out-neighbor. Let \(F=P_{B-v}\), \(O=\sum_{C\ni v,\ |C|\text{ odd}}P_{B-V(C)}\), and \(E=\sum_{C\ni v,\ |C|\text{ even}}P_{B-V(C)}\), summing directed cycles up to rotation. Then
\[
P_{B[v\leftarrow W]}(k)=(F(k)-2kE(k))P_W(k)
+k(O(k)+E(k))P_W(k+1)+k(E(k)-O(k))P_W(k-1). \tag{1}
\]
At most one cycle in a packing crosses the boundary of \(W\), since that cycle must use the unique outside entrance or exit. Its \(W\) vertices form one run. Contracting that run gives an old cycle whose length has the same parity as the required run size (the new cycle must be odd). The odd subset mass, with \(H(W[J]) P_{W-J}(k)\) as weight, is \(\tfrac12(P_W(k+1)-P_W(k-1))\); the even nonempty mass is \(\tfrac12(P_W(k+1)+P_W(k-1))-P_W(k)\). For the signed sum proof, expand odd-cycle packings in both sides: each odd block can be put in the distinguished slot with factor \(-1\) or one of \(k\) other slots. Thus \(P_W(k-1)\) is the signed sum. The crossing cycle contributes \(2k\); packings with no crossing contribute \(FP_W\). At \(k=1\), (1) reduces to the inherited ordinary one-outlet transfer.

The functional operator is injective on polynomials. On a nonzero input \(D\) of degree \(d\), \(FD\) has degree \(|B|-1+d\). The odd-cycle term is \(kO\) times the first difference of \(D\), of degree at most \(d-1\); the even-cycle term is \(kE\) times the second difference, of degree at most \(d-2\). Since the degrees of \(O,E\) are at most \(|B|-3,|B|-4\), respectively, these terms cannot cancel the leading term. Subtract and descend to recover each coefficient.

## 2. Same scalar, different continuation, same automorphism group

Let \(p=P_W\). The minus host \(F_r(W)\) is the triangle substitution \(C_3[A,W\to TT_r,y]\), \(A,y\) singletons, with \(A\to W\to y\to A\). The source is \(F_4(W)\). Adding a singleton beating only \(y\) gives \(F_5(W)\) up to isomorphism: in the near-transitive block order place the newborn just before \(y\). The one backward edge between the singletons gives
\[
P_{F_r(W)}(k)=k^{r+2}p(k)+k(k+1)^r p(k+1)-k(k-1)^r p(k-1).
\]

The plus host \(G(W)\) has quotient rows \((10,28,57,48,33,3)\), replacing vertex \(2\), with unique in-neighbor \(1\), by \(W\). In \(G^+(W)\) add a singleton beating only old vertex \(5\). The coefficients of \(p(k),p(k+1),p(k-1)\) in (1), from literal quotient cycles and complements, are:

| Host | \(a(k)\) | \(b(k)\) | \(c(k)\) |
|---|---|---|---|
| \(G\) | \(k^5+2k^3-2k\) | \(2k^4+4k^3+7k^2+3k\) | \(-2k^4+4k^3-7k^2+3k\) |
| \(G^+\) | \(k^6+4k^4-2k^2\) | \(2k^5+5k^4+11k^3+10k^2+4k\) | \(-2k^5+5k^4-11k^3+10k^2-4k\) |

Write \(h=p(1),q=p(2),s=p(3)\). The two source counts agree at \(H=h+16q\), and \(F_5,G^+\) give \(2H-h,2H+h\). When \(W=TT_m\), both exact binary signs occur from \(1+16\cdot2^m\). The minus move preserves the family under one transitive extension; the plus move is not such a return. Indeed even the second slot distinguishes both pairs:
\[
P_{F_4(W)}(2)-P_{G(W)}(2)=64s+20q+20h,\qquad
P_{F_5(W)}(2)-P_{G^+(W)}(2)=206s+8q+38h.
\]
The differences are positive. A continuation reading \(Q\) cannot freely treat the two scalar-equivalent predecessors as interchangeable.

For \(|W|\ge2\), all four constructions preserve \(\operatorname{Aut}(W)\). In the minus hosts the triangle quotient is indecomposable: a proper module cannot touch two parts, because it would then touch the third and include every vertex. Thus the nonsingleton part \(W\to TT_r\) is characteristic, the cyclic quotient stabilizer is trivial, and the ordered singleton suffix adds no automorphisms. In the plus hosts the six- and seven-vertex quotients are indecomposable and rigid, checked by all proper subsets and permutations. A module meeting \(W\) and the exterior projects to a module of the quotient: any excluded label is witnessed by an actual excluded vertex. If the projection is the whole quotient, an omitted vertex of \(W\) would see both an in- and an out-neighbor among the included exterior vertices, impossible. An exterior-only module likewise projects. This makes \(W\) the unique nonsingleton maximal proper module and the trivial quotient automorphism fixes every exterior label. Every automorphism of \(W\) extends. These constructions are strong because the quotients are strong and all modules are nonempty.

## 3. Exponential ambiguity at every order at least eight

The strong order-six reservoirs \(U=(14,52,56,50,33,1)\) and \(V=(6,44,56,49,35,1)\) have common polynomial \(p_0(k)=k^6+14k^4+18k^2\). At order six the only possible exponent drops in odd-cycle packings are two and four, so the common \(H=33,Q=360\) determines the polynomial. A literal permutation check gives \(\operatorname{Aut}(U)=\langle(1\,2\,3)\rangle\simeq C_3\) and \(\operatorname{Aut}(V)=1\).

For \(n\ge8\), let \(r=\lfloor(n-2)/6\rfloor,\ell=n-2-6r\). Take any binary word of length \(r\) over \(U,V\), let \(W\) be the ordinal sum of the letters followed by \(TT_\ell\), then form \(F_0(W)\). Every output is strong and has order \(n\). The common reservoir polynomial is \(p_0^r k^\ell\), and common output count is \(33^r+360^r2^\ell\). All \(2^r\) outputs are pairwise nonisomorphic: the large part is characteristic, and the canonical order of strong components records the word. Exactly \(\binom rj\) of these isomorphism classes have group \(C_3^j\). Full slot equality follows by the functional formula.

The all-order existence boundary can also be made exact. At order seven the literal pair
\[
(78,52,56,50,33,65,30),\qquad (70,108,56,49,99,65,12)
\]
is strong, has common profile \(k^7+24k^5+86k^3+24k\), and group orders \(3,1\), respectively. The order-six pair, this order-seven pair, and the wrapper together give different groups at a common full profile for **every order \(n\ge6\)**. Exhaustion over all labeled strong tournaments through order five finds no equal-profile group-order discrepancy there; all groups in this small range have order \(1,3\), or \(5\), so equal order also identifies the group. The five-or-fewer verification is a finite boundary check, not an ingredient in the construction.

The same word bank in the two matched binary hosts gives a common predecessor \(33^r+16\cdot360^r2^\ell\) and offsets \(\pm33^r\) around twice that predecessor. For each sign there are \(2^r\) inequivalent strong witnesses, with identical whole slot polynomial within that sign and the stated group distribution. Taking \(W\) transitive gives unit offset but no internal symmetry freedom.

## 4. Sufficiency under arbitrary substitution, and a ternary symmetry recursion

The full slot polynomial is sufficient under **any lexicographic substitution**. In each of \(k\) ordered path slots, collapse every run occupying one module to its label. Consecutive labels follow quotient arcs and are distinct. Fix a \(k\)-tuple of such quotient words, letting \(j_i\) count all occurrences of label \(i\). Its contribution is \(\prod_i L_{W_i}(j_i)\), where the ordered cover count for nonempty runs is
\[
L_W(j)=\sum_{a=0}^j(-1)^{j-a}\binom ja P_W(a).
\]
This is inclusion–exclusion for empty slots. The skeleton is recovered from any final path cover, so there is no overcount. Only \(j_i\le |W_i|\) matter. For the directed triangle, the skeleton counts are coefficients of
\[
\left(\frac{1+x+y+z+xy+yz+zx+2xyz}{1-xyz}\right)^k,
\]
marking occurrences of the labels \(x,y,z\); the single-slot series follows by choosing the starting label and word length.

For height \(t\), substitute the directed triangle at every internal node of a balanced ternary tree, allowing \(U,V\) independently at every leaf. All resulting tournaments are strong of order \(6\cdot3^t\) and share the full slot polynomial by the sufficient-state rule. A node type is its triple of child isomorphism classes modulo cyclic rotation, since the maximal modules are exactly its children. The automorphism group is the direct product of the child groups unless all three child types are equal, in which case it is \(\operatorname{Aut}(Z)^3\rtimes C_3\) with cyclic permutation for their common type \(Z\).

Writing \(J_t(x)=\sum x^{\log_3|\operatorname{Aut}|}\) over the resulting isomorphism classes yields
\[
J_0(x)=1+x,\qquad
J_t(x)=\big(J_{t-1}(x)^3-J_{t-1}(x^3)\big)/3+x J_{t-1}(x^3).
\]
The nonconstant triples have orbit size three, and the constant ones gain one group factor. The first class totals are \(2,4,24,4624\) at heights \(0,1,2,3\).

- At \(t\ge1\) the minimum exponent is \(3^{t-1}\), attained by exactly two classes with group \(C_3^{3^{t-1}}\). At height one the two minima are \(VVV\) (one quotient rotation) and \(UVV\) (one leaf rotation). At every next height the only minimum triples are the two nonconstant cyclic types on the previous two minima. Constant triples gain an extra factor.
- The maximum exponent \(M_t=(3^{t+1}-1)/2\) is uniquely all \(U\) leaves. The next exponent is \(M_t-(t+1)\), uniquely one \(V\) leaf. For the induction, a nonconstant triple loses the root rotation and at least the previous second-place gap; a constant nonmaximal triple loses at least three times the previous gap. Equality in the tighter bound occurs precisely with one second-place and two maximum children. This gives one cyclic type.

Thus the whole slot profile is an exact composable count state, but it is neither an isomorphism state nor an automorphism state. Symmetry choices do not by themselves add scalar counts to the binary gate.

**Replay:** python3 replay_ports2.py checks the coefficient table, independent small full slot evaluations, complete quotient modules and permutations, the order-six palette, and the ternary group recurrence against direct cyclic canonicalization.

## 5. A physical internal reader can distinguish hidden point choices

For a marked \(x\in W\), define \(S_x(k)\) and \(T_x(k)\) by summing \(P_{W-V(L)}(k)\) over directed paths \(L\) starting or ending at \(x\), respectively. Let \(D_x(k)\) be the signed version of \(S_x(k)\) with sign \((-1)^{|L|}\), and \(f_x=S_x(0)\) the full starting count. If \(W\) is strong and \(x\in W\), build \(Y_r(x)\) by putting \(TT_r\) after \(W\) in ordinal order and then adding a vertex \(z\) beating only \(x\). Its complete output polynomial is
\[
P_{Y_r(x)}(k)=k^{r+1}P_W(k)+k(k+1)^r S_x(k)+k(k-1)^rD_x(k).
\tag{2}
\]
The first term keeps \(z\) singleton. Every odd cycle through \(z\) is \(z,L\), then an optional increasing path on a subset of \(TT_r\), closing at \(z\). Its prefix after \(z\) must have even order, so the total mass is half the unsigned-plus-signed sum, with the boundary of one vertex automatically excluded. The cycle weight \(2k\) gives (2). At \(r=0\) the signed term is retained (use the polynomial convention \(0^0=1\)); at \(r\ge1,k=1\) it vanishes. The output is strong whenever \(W\) is strong of order at least three.

Use the two order-seven homogeneous-profile reservoirs from §3, mark vertex \(4\) in the first and \(6\) in the second. The two sets of values are
\[
(S_x(1),D_x(1),S_x(2),D_x(2))=(165,-27,993,-285),\quad(165,-25,997,-273).
\]
At \(r=0\) the two physical outcomes are \(273,275\), adjacent odds, despite identical complete predecessor profiles. For every \(r\ge1\) they coalesce on the exact unbounded count \(135+165\cdot2^r\), while their second-slot counts differ by \(8\cdot3^r+24\). The output groups remain \(C_3\) and \(1\). For \(r=0\) the new \(z\) is uniquely of outdegree one; for \(r\ge1\) only \(z\) and the last transitive vertex have outdegree one, and only the latter points to an outdegree-one vertex. Fixing \(z\) fixes the marked port; the stabilizer of that port in the first reservoir is the full \(C_3\), and in the other it is trivial. Thus both agreement and disagreement of scalar outputs can be forced by the positioning of a physical reader, with distinct symmetry types.

There is also a two-growth-scale exact reader inside \(F_r(W)\). Put \(T=F_r(W)\) and retain \(x\in W\). Then the ordinary increment of adding \(z\) beating only \(x\) is
\[
2E_x(T)=2^{r+1}S_x(1)+2\cdot3^r S_x(2)-f_x-T_x(1),
\tag{3}
\]
where \(E_x(T)\) is the even-start weight with residual \(H\). For the all-start count in \(T\), split according to how the path encounters \(y,A\). If it uses neither, the complement has remaining \(A,y\) and a two-slot readout of the remaining \(W\) and padding. This contributes \(2^r S_x(1)+3^r S_x(2)-f_x\): subtract \(f_x\) for the case the entire middle is selected, where \(A,y\) form only an edge. If it reaches \(y\) but not \(A\), it gives \(2^r S_x(1)\). If it uses \(y,A\), each padding vertex independently joins before \(y\), after \(A\), or to the complement, and each remaining \(W\) slot splits into an optional path after \(A\) and the complement, giving \(3^r S_x(2)\). A full path ending \(x\) in \(T\) must use \(y,A\) consecutively, and all padding precedes \(y\). This gives the weighted end count \(T_x(1)\) from \(W\). Subtract it from the total-start weights using rooted signed reciprocity to get (3). Equal two start profiles with a differing final sum \(f_x+T_x(1)\) would lock a constant shift at all widths; this is a sufficient recipe, not a claim that every such signature has a physical witness.

For \(r\ge2\) the new \(z\) in (3) is the only outdegree-one vertex: \(y\) and each padding vertex now beat \(z\) and another vertex, while \(W\) vertices have outdegree at least \(r+1\). Consequently the output group is exactly \(\operatorname{Aut}(W)_x\). In the ternary bank, a marked leaf orbit can lose many ancestral cyclic rotations without changing the unrooted input slot state. Its scalar output is then controlled by (3), not by the unrooted profile alone.
