# Neutral sign pairing at the full pointed interface

The starting point is [Session 313's exact two-count offset](../session313_zeckendorf_target_selection_20261113/ADDITIVE_PORTS.md). The new point is that the two cuts differ only in one mixed response **at every slot count**. Pairing two switches with opposite signs transports a chosen incidence without changing the full cover polynomial or the named apex interface. A probe at another named vertex reads exactly where the switch landed.

## 1. A full-polynomial local identity

Let $R$ be the cyclic order-five tournament with outgoing masks $(24,17,3,6,12)$ (zero-based); equivalently $i$ beats $i+3,i+4\bmod5$. A future apex beats the cut $U$. For an integer $t\ge0$, define:
- $S_U(t)$: ordered covers of $R$ with a distinguished nonempty path starting in $U$, and $t$ further ordered slots allowing empty paths;
- $E_U(t)$: same, with the distinguished path ending in the complement of $U$;
- $B_U(t)$: two distinct ordered nonempty distinguished paths, one ending in the complement of $U$ and the other starting in $U$, and $t$ further unrestricted slots.

At masks 3 and 5, these responses are
$$
\begin{aligned}
P(t)&=t^4+2t^3+8t^2+7t+3,\\
S_3(t)=S_5(t)&=2P(t),\qquad E_3(t)=E_5(t)=3P(t),\\
B_3(t)&=2(t+1)(3t^2+6t+13),\\
B_5(t)&=B_3(t)+2(t+1).
\end{aligned}
$$
Also $Q_R(z)=z(z^4+10z^2+4)$, where $Q$ counts all ordered slots including empties. The equalities can be checked by elementary subset decomposition: assign vertices to the distinguished paths, count their starts/ends by subset DP, and put the remainder in the unrestricted slots. The degree is at most 4 when one distinguished path is nonempty and at most 3 when there are two, because choosing the indices of the nonempty slots gives a binomial-polynomial expansion. So evaluations at $0,1,2,3,4$ certify the identities. The standalone replay supplies those exact evaluations.

## 2. Mirror pairs give exponentially many exact same-count strong tournaments

Take an ordinal sequence of $n=2d$ blocks, first $d$ copies of $R$, then their converses in reverse order. Add a common apex $x$. In left block $i$ choose out-cut $U_i=3$ or $5$; the mirrored right block uses the complementary out-cut in the same local vertex labels. The map reversing the block order, preserving each local vertex label and fixing $x$, is an anti-involution. The tournament is strong: the apex reaches the first strong block, the last strong block reaches the apex, and cross arcs between blocks run forward.

Count ordered covers with $z\ge1$ slots. Choose the one containing the apex. If its immediate predecessor and successor lie in different blocks, each of those blocks uses only a single constrained endpoint, with the other slots ordinary. Their responses agree for the two cuts at **every** number of slots. If one or both immediate neighbors do not exist, the same conclusion holds. The only possible change is when both neighbors lie in the switched block. Then that block contributes $B_U(z-1)$. In each other block the apex-containing path can visit on only one side of the apex, by the ordinal order, so its $z$ slots contribute $Q_R(z)$. Summing over the distinguished slot gives exactly
$$
\Delta Q_T(z)=z\,[B_5(z-1)-B_3(z-1)]Q_R(z)^{n-1}
             =2z^2Q_R(z)^{n-1}
$$
for one switch, including the dual switch in a right block. Thus a mirror pair adds $4z^2Q_R(z)^{2d-1}$. With $k$ of the $d$ pairs switched, the exact count is $H_{d,0}+4k15^{2d-1}$, and the full polynomial depends only on $k$, not on which pairs were chosen. Switching two opposite bits is an exact neutral operation.

These are distinct isomorphism classes, and are rigid. The apex is the unique vertex whose deletion leaves a nonstrong tournament. Deleting it leaves the ordinal chain of strong blocks. Deleting a vertex in any interior block still leaves paths through the first block, last block, and apex to connect everything. For the first block, starting from either surviving out-cut, directed reachability after deleting local vertex $v$ gives all other four vertices (the reached masks for $v=0,1,2,3,4$ are $30,29,27,23,15$ for both cuts). The last block follows by converse. All vertices can then reach the last block and the apex. Any isomorphism must fix $x$, hence preserve the ordered SCC positions after deleting $x$. The only automorphisms of the cyclic block are its rotations: a fixed vertex fixes its two ordered out-neighbors and its two ordered in-neighbors. The cuts $\{0,1\}$ and $\{0,2\}$ are not rotation-equivalent, and neither admits a nonidentity stabilizer. Therefore each $k$ gives exactly $\binom dk$ distinct members within this bank at one count and one polynomial.

The equality survives a **named selective interface**. With the slot containing $x$ specified, distinguish four cases: $x$ is the only vertex, the initial but not final vertex, the final but not initial vertex, or interior. The first three counts are invariant under all switches because no block supplies a doubly constrained mixed split. The interior count is invariant for fixed weight by subtracting those three from $Q_T(z)/z$. Requiring the other slots to be nonempty is inclusion–exclusion. Thus when an exterior sees the bank uniformly except possibly at $x$, each exterior itinerary has the same number of fillings: the order of the bank runs is prescribed, and the legality of exterior incidences depends only on which run has $x$ and its endpoint status there. This applies to one full path or to ordered path covers. The result does not give independence for arbitrary additional named internal ports.

## 3. A probe that recovers a chosen bit, with both signs

Fix any left block position $k$, with $d\ge2$. Add a vertex $p$ at that ordinal position, beating the common apex, beating all five local vertices except local vertex 0 (which beats $p$), and inheriting the ordinary ordinal edges to every other block. This is a literal six-vertex replacement in that slot. Its internal $h,q$ values are $33,360$ for either cut. The five scalar responses $(s,S,e,E,b)=(S_U(0),S_U(1),E_U(0),E_U(1),B_U(0))$ are
$$
(8,93,24,231,69)\quad(U=3),\qquad
(9,94,23,226,69)\quad(U=5).
$$

Compare two unprobed words of the same weight with bit $k$ changed from 3 to 5. All other changes have net $-1$ among the other $2d-1$ single sites, including the right-hand mirrors. In the probed graph those compensating switches change only $b$, so contribute $-2\cdot33\cdot15^{n-2}$ in total. The ordinary one-apex path formula, partitioning by the last block before the apex and first after it, gives the exact change
$$
\frac{\Delta H_{\rm probe}}{15^{n-2}}
 =c_k+\sum_{i=k+1}^{n}E_i\,8^{i-k-1}
 -210\sum_{r=0}^{k-2}8^r,\qquad
c_k=15[k=1]+9[k>1]-s_{k+1}-66 .
$$
Here $s_i=6,E_i=63$ on the left, and $s_i=9,E_i=42$ on the right. To see each term, the empty-left case contributes $\Delta s$ only at $k=1$; consecutive one-visit boundaries contribute $e_{k-1}\Delta s$ and $\Delta e\,s_{k+1}$; the double-visit cases contribute $\Delta S\,E_i$ to the right and $S_j\Delta E$ to the left. Every intermediate double slot contributes 120 rather than 15, accounting for the ratio 8. The probe's internal $h,q$ do not depend on the bit, so all unaffected terms cancel.

This is always positive. For $k=1$, the first positive summand is 63, exceeding $-c_1=57$. For $k\ge2$, $n-k-1\ge k-1$, so the two sums have difference at least
$$
42\left(8^{k-1}-4\sum_{r=0}^{k-2}8^r\right)
=42(3\cdot8^{k-1}+4)/7\ \ge168,
$$
and $c_k\ge-66$. A probe at local vertex 1 instead gives the same $h,q$ and response change $(-2,-13,-3,-22,-15)$. Every affected term is then nonpositive and the compensation is negative, so this probe has strictly negative response. These are opposite signs, not a claim of equal magnitudes.

At $d=2$, the two words $(3,5)$, $(5,3)$ have $H=33{,}787{,}125$ and $Q(2)=25{,}979{,}389{,}200$. The vertex-0 probe in the first block gives respectively $74{,}776{,}365$, $75{,}458{,}115$; their difference $681{,}750$ agrees with the formula. These values were independently checked by literal 21- and 22-vertex subset DP. This realizes selective sign control at an exact pre-existing count with a reusable apex interface, but does not claim arbitrary targets as the common count.

## 4. Why the obvious strong $7\cdot3^j$ growth cannot just add a triangle

For any partition into nonempty blocks $B_i$ of an actual tournament $T$, the positive odd-cycle packing expansion gives
$$
H(T)\ge\prod_i H(T[B_i]).
$$
Indeed packings lying within blocks contribute exactly the product. Equality holds exactly when every SCC of $T$ lies in one block. If there is a mixed directed triangle, its singleton packing adds 2. Conversely the connected components of the directed-triangle overlap hypergraph, including unused single vertices, are exactly the SCCs. A union of overlapping triangles is strong. Any outside vertex seeing such a strong union in both directions would form a triangle with an arc from its out-neighbors to its in-neighbors inside the union, joining the component. So these components are modules and their quotient is transitive (a quotient triangle would lift). Thus a crossing SCC forces a mixed triangle, and if no SCC crosses, no mixed directed cycle is possible.

In particular reorienting only the cross arcs of the exact ordinal construction $C_3\to T\to C_3$ cannot make it strong at count $9H(T)$. This rules out all such strongifications, not just single cross-arc flips or those respecting converse. It leaves changing the old interiors or a different factor certificate open. Direct finite scans also found no extension by one mirrored pair (keeping the old anti-involution and induced old graph) from Paley $7$, $H=189$, to $567$ or $1701$, or from the supplied $9$-vertex self-converse $H=567$ witness to $1701$ or $5103$. Those failures are scoped to the chosen witnesses and extension protocol, not all strong self-converse tournaments at those counts.

## 5. Which vertex-transitive modules can carry such a switch?

For any vertex-transitive core the single-end responses $S_U(t)$ and
$E_U(t)$ depend only on the size of $U$. Let $K_{b,a}(t)$ count two
disjoint distinguished runs ending at $b$ and starting at $a$,
with $t$ other slots (so $K_{a,a}=0$). All column sums are the same,
say $c(t)$, by vertex transitivity. Then
$$
 B_U(t)=|U|c(t)-
        \sum_{\{a,b\}\subseteq U}\big(K_{b,a}(t)+K_{a,b}(t)\big).
$$
Thus a same-size cut change is precisely a change in weighted pair
energy in the induced subset. If the core action is transitive on
unordered pairs as well, every pair has the same weight and the
mechanism cannot distinguish any two cuts of the same size. In
particular the Paley 7 action by translations and square multipliers
is transitive on unordered pairs: the sign change switches squares
and nonsquares. The cyclic five-vertex core instead has two pair
orbits, visible to the two masks. This gives a test for new modules
rather than a coincidental match of two small outputs.
