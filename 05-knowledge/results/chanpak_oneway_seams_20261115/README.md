# One-way Busch seams cannot produce the first two stable excesses

**Proved.** Suppose a strong tournament \(T\) has \(H(T)=f(n)+K\), with \(K\in\{2,4\}\) and \(K<2f(n)/25\). In **every** strong one-vertex Busch split of \(T\), the cross-arcs between the two factors off their common vertex go in both directions. Thus the one-way gluing of two exact floor factors is not a candidate for either absolute target at this order. The mixed-cross-arc case is open.

Notation and the reduction to two exact floor factors are from the [near-floor stability theorem](../chanpak_busch_stability_20261115/README.md). By [Moon–Yang, §§2 and 4](https://citeseerx.ist.psu.edu/document?doi=78c085212810b466efaf3f903baeb284a695f719&repid=rep1&type=pdf), each factor is a special chain of nearly transitive links. In its chain order, all arcs go forward except one bottom-to-top reversal in each link; successive links share the bottom of the earlier link and the top of the later one. Link orders belong to \(\{3,4,5\}\). Here path ends and starts count with direction, never modulo converse.

## 1. Exact surplus for a one-way seam

Let the common vertex be \(v\), and assume \(A-\{v\}\) beats \(B-\{v\}\). In each factor take an ordered two-path cover \((P,Q)\) of the vertices other than \(v\), allowing empty paths. A prefix has status \(0\) if empty, \(+\) if its last vertex beats \(v\), and \(-\) otherwise. A suffix has status \(0\) if empty, \(+\) if \(v\) beats its first vertex, and \(-\) otherwise. Let \(a\) count \(A\)-covers of type \((-,0)\) or \((-,+)\), and \(a_+\) only the second type. Let \(b\) count \(B\)-covers of type \((0,-)\) or \((+,-)\), and \(b_+\) only the second. Write \(E_A\) for the number of \(A\)-paths ending at \(v\) and \(S_B\) for the number of \(B\)-paths starting there. Then the local path surplus is
\[
 L=H(T)-H(A)H(B)
   =a\bigl(H(B)-S_B\bigr)+b\bigl(H(A)-E_A\bigr)+a_+b_+ .
\]
Indeed, before \(v\) and after \(v\), the portions from \(A,B\) must occur in that order. Each whole path is uniquely \(P_AP_BvQ_AQ_B\). If only the restriction to \(A\) is broken, the failure must be its prefix, and a nonempty valid \(B\)-prefix screens it. The \(A\)-suffix must be empty or valid. This counts the first product. The second is the dual case. If both break, \(P_B\) and \(Q_A\) must both be nonempty and valid at \(v\), yielding the final product. Unbroken pairs give exactly \(H(A)H(B)\) paths.

In a special chain, every start or end at a specified vertex has at most \(4/9\) of the total path count. To verify this, its decomposition into links is exact multiplicative, and every pair of paths at one boundary has a forced unique join, by the preceding order rule and the chain path product. Thus a uniform chain Hamiltonian path projects to a uniform Hamiltonian path in each link. If the full path starts or ends at \(v\), so does its projection to any link containing \(v\). Direct counts in links 3, 4, 5 give maximum endpoint shares \(1/3,2/5,4/9\), respectively. Hence a nonzero \(L\) is at least \((5/9)\min(H(A),H(B))\). If \(B\) is a triangle, \(b=0\) and \(L=2a\), as follows directly from its two other vertices. By converse one can always arrange that a triangular factor is \(B\).

## 2. All small broken-prefix profiles of an exact floor factor

For a special chain \(R\) rooted at \(v\), write \(O=N^+(v)\). Its \(a=a_0+a_+\) profile counts a nonempty prefix ending in \(O\) and a suffix which is empty (\(a_0\)) or starts in \(O\) (\(a_+\)). Apart from the overall bottom root with \(a=0\), every \(a\le2\) occurs in these exceptions, listed in full chain order:
\[
\begin{array}{c|l}
\text{link orders}&a\text{ at all vertices}\\ \hline
(3)&0,0,0\\
(4)&2,4,1,0\\
(5)&6,16,7,1,0\\
(3,3)&12,14,2,2,0 .
\end{array}
\]
The four rows are literal two-path cover counts and are included in the replay. The exclusion of other small profiles works at every order, as follows.

* **Overall bottom.** Its only outneighbor is the top of the last link. This vertex cannot be both prefix end and suffix start. After removing the bottom, the ordinary chain order has a terminal transitive segment after that top, so no Hamiltonian path ends there. Thus \(a=0\).
* **Root outside the last link.** All vertices of the last link lie in \(O\). The earlier vertices minus \(v\) have a Hamiltonian path in their natural order, except for swapping the endpoints of a triangle whose middle vertex was removed (if both endpoints are present). The last vertex of that path beats every vertex in the last link except possibly its top; the sole possible failure is the reversal from that top into the previous link. Hence append any Hamiltonian path of the last link except possibly the one starting at its top, giving at least two covers of type \((-,0)\). The top-started path in a nearly transitive link is unique. Also append only a single nontop vertex of the last link to the earlier path, and put all remaining vertices of that link into the suffix in a Hamiltonian order. Both endpoints lie in \(O\), giving one additional \((-,+)\) cover.
* **Root strictly inside the last link, with an earlier chain \(P\).** If the link has order at least 4, deleting \(v\) leaves a shorter nearly transitive link of order at least 3. Its bottom is in \(O\), and the unique path of that shortened link ending at its bottom extends with *every* Hamiltonian path of \(P\), by the exact chain merge at the link's top. Thus \(a_0\ge H(P)\ge3\). If the last link is a triangle, after deleting its middle \(v\), its bottom can end a Hamiltonian path exactly when the preceding path in \(P\) does not end at the common vertex \(t\): the bottom loses to all earlier vertices except \(t\). Here \(a_0=H(P)-E_P(t)\). This is 2 when \(P\) is a single triangle, and at least 3 otherwise by the endpoint share bound and \(H(P)\ge5\).
* **Root at the top of the last link, shared with a preceding chain.** Let \(w\) be the top of the preceding chain's last link, so \(v\to w\). If the final link has order 4 or 5, the local \(a_+\) at its top is at least 2 from the single-link counts. Prepend a Hamiltonian path of all earlier vertices to each local prefix. There is one additional cover: place all earlier vertices except \(w\), then the first intermediate of the last link, in the prefix; use \(w\) followed by a Hamiltonian path of the remaining last-link vertices in the suffix. This one has a different suffix start. If the last link is a triangle, with remaining vertices \(z\to d\) and \(v\to z\), choose a suffix starting at \(w\), followed by any chosen subset of the preceding link's intermediate vertices in increasing order, then \(d\). Place a Hamiltonian path of the unused earlier vertices, then \(z\), in the prefix. If the preceding link has order at least 4 there are already at least four choices. If it is also a triangle, two choices arise this way. With at least one more earlier link there is a third choice: start the suffix at \(w\), use its backward arc to the top of that earlier adjacent link, then go directly to \(d\). The remaining case is exactly the two-triangle row.

All paths of unused vertices in these constructions exist by Rédei's path existence theorem, including the empty path when needed. The arc conditions used to attach them follow from the chain order.

## 3. Eliminate the one-way branch at the prescribed excesses

In the stability band the two root floors multiply to \(f(n)\), so the target excess equals \(L\).

If \(L=2\), the endpoint bound forces a triangular factor. Then \(a=1\), so the other factor is a single link of order 4 or 5. Either way the total order is at most 7, where \(2<f(n)\,2/25\) does not hold.

If \(L=4\), at least one factor has order at most 4. If it is a triangle, \(a=2\), forcing the other factor to have order at most 5. If neither is a triangle, up to converse take \(A\) of order 4. If \(H(B)\ge9\), positive \(a\) would alone cost at least 5, hence \(a=0\). Its bottom has \(E_A=1\), so \(L=4b\). Therefore \(b=1\), which by the dual of §2 forces \(B\) to have order 4 or 5. If \(H(B)<9\), it already has order 4. All cases have \(n\le8\), where \(4<f(n)\,2/25\) does not hold. This proves the claimed obstruction.

The [exact replay](verify.js) checks 273 marked special-chain profiles on 30 chains, including the four exceptional rows, and evaluates 1,303 marked one-way joins by literal subset DP. The all-order exclusion relies on the structural argument, not that range. It constrains the top seam for these two exact targets; it is not a nonattainment theorem without treating mixed cross-arcs.
