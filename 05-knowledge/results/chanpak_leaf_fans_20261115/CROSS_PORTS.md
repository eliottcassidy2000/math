# Cross-port response and obstruction to two fixed fans

Write $H=H(T)$, with empty count one. Add a leaf $x$ by making it beat only its designated root $q$ among the old vertices, and write $L_qT$ for the result. The defect is $\delta_q=H(L_qT)-2H$.

## Exact cross-port law

Let $E_v,O_v$ be the sums, over nonempty paths starting at $v$ of even and odd vertex order respectively, of the Hamiltonian count of the remaining vertices. Let $J_{pq}$ be the same sum over odd-order paths from $p$ to $q$. For $p\ne q$, $D_{rs}$ sums $H(T-V(P)-V(Q))$ over ordered disjoint old paths $P,Q$ starting at $q,p$, of respective parities $r,s$.

Then $H(L_qT)=H+2E_q$, and
$$
\boxed{\delta_p(L_qT)-\delta_p(T) =
-2E_q+2(O_p-J_{pq})+2(D_{eo}+D_{oe})+4D_{ee}.}
$$
In particular, attaching one leaf at each of $q,p$ in opposite orders changes the resulting count by
$2[(O_p-J_{pq})-(O_q-J_{qp})]$ (order $q,p$ minus $p,q$).

**Proof.** Use the odd-cycle packing expression for $H$: disjoint odd cycles of length at least three have weight two each, uncovered vertices weight one. This follows from the [Grinberg–Stanley tournament expansion](https://arxiv.org/abs/2307.05569). Here cycles are counted up to rotation. The term from packs using one new vertex with old vertices is $2E_q$ or $2E_p$, since the outgoing arc is forced and the old path must have even order. For two new vertices $x,y$, added in that order, $x\to y$. If they are consecutive on the same cycle, the cycle uses $x\to y\to p$, followed by an old odd-order path ending anywhere except $q$, so its contribution is $2(O_p-J_{pq})$. If separated on the same cycle, their old paths start at $q,p$, are disjoint, and have opposite parities; disjointness also guarantees the old incoming arcs to $x,y$ are present. If in different cycles, the two old paths both have even order. In each case the factor for unused old vertices is their own Hamiltonian path count. Subtracting $2H(L_qT)$ yields the formula.

## Third visits destroy fixed opposite signs

For distinct $p,q$ let $T_k=L_q^k T$. For $k\ge1$ there exist nonnegative integers $a,b,c$, independent of $k$, such that
$$
H(L_pT_k)=a+b(2^k-2)+c(3^k-3\cdot2^k+3).
$$
Moreover, $c>0$ **if and only if $q$ beats an old vertex outside $\{p,q\}$**.

Indeed the $k$ leaves form a transitive module $F$. With the latest added vertex denoted $y$, the only outneighbors of $F$ outside $F$ are $q,y$. Compress every maximal visit to $F$ to a placeholder. There are at most three visits. A compressed path with $j$ placeholders has exactly $j!S(k,j)$ lifts: partition the labelled fan vertices into $j$ ordered nonempty subsets, and use their unique transitive orders within the visits. A three-visit skeleton exits to $q,y$ once each. The $y$-excursion must start $y\to p$; if the $q$-excursion is to be disjoint it needs an outneighbor $r\ne p$ in the old tournament. Conversely, given $q\to r$, a skeleton is $F,q,r,F,y,p,F$, preceded by a Hamilton path of all remaining old vertices. Its final vertex, if any, beats $F$.

In contrast $H(T_k)$ is affine in $2^k$ for $k\ge1$ by the very same compression with only one outneighbor for $F$. Therefore if $c>0$, then $\delta_p(T_k)\to +\infty$. There is a finite test too: along $k=1,2,3$,
$$
\delta_p(T_3)-3\delta_p(T_2)+2\delta_p(T_1)=6c.
$$
Thus all three defects cannot equal a fixed sign.

**Corollary.** No two distinct permanently designated vertices maintain defects $+1,-1$ under arbitrary words of attachments at those vertices. More strongly, stability of both signs through the three consecutive repeats at either root is impossible. Otherwise the preceding test in both directions forces each old root's outneighbors to lie in the other singleton. One root is then a sink $r$ and the other $s$ beats only $r$. The sink has defect $-H$, so $H=1$ and the tournament is transitive, since any directed triangle would give an additional positive term in the odd-cycle formula. Attaching $k$ times at $r$ keeps the tournament transitive and makes $s$ have $k+1$ outneighbors. A path starting at $s$ in the new transitive tournament is uniquely specified by its subsequent vertices; exactly $2^k$ choices give even-order paths. Thus $\delta_s=2^{k+1}-1$, already $3$ for $k=1$. This does not rule out reselecting an armed root after an attachment.

The counting basis is the same ordered-surjection count as in the prior outlet-width argument; the root-specific nonzero coefficient test and finite obstruction use the exact forced arc $y\to p$.

## Exact transport of repeat-stable ports

Call $p$ repeat-stable when every Hamiltonian path of $T-p$ ends in an in-neighbor of $p$; equivalently, its defect is unchanged by its own leaf append, by the same-port law. For a leaf $x$ at $q$, the entire repeat-stable menu of $L_qT$, for old order at least two, is:
- $q$ always.
- A distinct old $p$ if and only if the old outneighborhood of $q$ is exactly $\{p\}$.
- $x$ if and only if $q$ lies outside the final strong component of $T$.

For the first claim, deleting $q$ makes $x$ a new sink: every remaining Hamilton path ends at $x$, which beats $q$. For the second, if $q$ is not a sink in $T-p$, that tournament has a Hamiltonian path ending at some vertex other than $q$: possible ends are all vertices of the terminal strong component, by the ordered component decomposition and a Hamiltonian cycle on a nonsingleton strong component. Append $x$; its endpoint is beaten by $p$, ruling out stability. If $q$ is a sink of $T-p$, every Hamiltonian path of $(L_qT)-p$ ends $x,q$, so stability holds exactly when $q\to p$. The newborn case follows immediately since its only old outneighbor is $q$.

In fact a repeat-stable vertex must belong to the final component: a vertex outside it beats every member of it, and after deleting that vertex the final component is unchanged. Consequently when every move chooses a repeat-stable port, no newborn is eligible. One can change to another eligible port only when the currently chosen $q$ has a unique outneighbor $p$. Stability of $q$ implies $p$ was not a sink (otherwise all Hamiltonian paths of $T-q$ would end at $p$). The new leaf is beaten by $p$, so the outdegree of $p$ is now at least two and cannot decrease. After changing to it, no other port is available again. **Every walk restricted to repeat-stable ports has at most one change of root.** This does not bound raw defect $+1$ or $-1$ ports that lack individual repeat-stability.

### Exact classification of repeat-stable unit switches

Suppose a switch from a repeat-stable port $q$ to another repeat-stable port $p$ is possible after $k\ge1$ repeated appends at $q$. Write $A=T-q$. The menu forces $N_T^+(q)=\{p\}$, so $T=L_p A$, and old stability of $q$ means $p$ lies outside the terminal component of $A$. Consequently no Hamilton path of $A$ ends at $p$: the end count $d_p(A)$ is zero. By the newborn rule, $\delta_q(T)=2d_p(A)-H(A)=-H(A)$. It is a unit only when $A$ is transitive, and then it is $-1$.

Let $T_k=L_q^kT$, with fan $F$ of size $k$. An odd cycle through a nonempty fan subset exits it to $q$, then follows a nonempty old path starting at $p$ before returning into the fan. If the old path has even order, the fan subset must have positive even order (there are $2^{k-1}-1$ choices); otherwise there are $2^{k-1}$ odd subsets. The terms using no fan give $H(A)+2E_p(A)$. Hence $H(T_k)=H(A)+2^k(E_p(A)+O_p(A))$ for all $k\ge1$, even if $q$ was not repeat-stable, as long as $N^+_T(q)=\{p\}$. A path starting at $p$ first follows an old path in $A$, then any subset of $F$ in transitive order, then optionally $q$. Paths may never leave $F$ except to $q$, and may never leave $q$ since its only outneighbor $p$ was used already. For each old path, the $k$ fan positions and optional $q$ give $2^k$ even-order choices, except exclude the impossible path $p,q$. In all cases the remaining old vertices, remaining fan and $q$ if present occur in that order in the complement, so the complement count is just that of the remaining subset of $A$. Hence
$$
E_p(T_k)=2^k(E_p(A)+O_p(A))-H(A-p),
\qquad
\delta_p(T_k)=2^k(E_p(A)+O_p(A))-H(A)-2H(A-p).
$$
If the old stable port carried a unit, $A$ is transitive and the old outdegree $d$ of $p$ in $A$ is positive. Thus $E_p(A)+O_p(A)=2^d$ and the new defect is $2^{k+d}-3$. **The only repeat-stable unit switch is from $-1$ to $+1$, with $k=1$ and $p$ penultimate in the transitive $A$.** Conversely that case works for every transitive $A$ with at least two vertices. Its count goes $3\to5\to11$ through the first switch, preceded by count 1 in $A$; extra leading vertices do not change those counts. The unit-only stable grammar cannot branch into a full binary tree.

Equivalently every fixed-seed stable-unit trajectory consists of a single dyadic ray $H_j=2^j(H_0+\delta)-\delta$, with at most the exceptional transition from the ray starting at $3$ to the $+1$ ray starting at $5$. The latter has totals $6\cdot2^j-1$. This rigidity concerns a proof architecture, not the overall Hamiltonian spectrum.
