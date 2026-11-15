# Sharp outlet costs, and a same-source pair that cancels them

A single one-outlet operation cannot increment an unbounded Hamiltonian count by a sublinear amount. But changing the outlet **within one strong component** can cancel two large responses down to a literal unit in the even shift. The two statements are compatible, and both hold for physical tournaments, not only formal transfer values.

## Notation and inherited identities

Write $H(T)$ for the number of directed Hamiltonian paths, including $H(\varnothing)=1$, and $Q(T)=\sum_{J\subseteq V(T)}H(T[J])H(T-J)$ for the two-slot count with empties. For a root $u$, let $s_u,O_u,E_u$ count weighted paths starting at $u$, unrestricted in size, odd in size, and even in size respectively. Each directed path $P$ is weighted by $H(T-V(P))$. Thus $s_u=O_u+E_u$. Let $t_u$ be the ordinary full Hamiltonian count ending at $u$.

We use the repo's [squarefree rooted reciprocity](../chanpak_ternary_emptying_20261115/ROOTED_SLOT_RECIPROCITY.md) in the form $O_u-E_u=t_u$. For completeness, in squarefree variables let $P$ be 1 plus the path series and $R_u^\pm$ be the rooted start/end series. The identity $R_u^+=-P(x)R_u^-(-x)$ follows by expanding a path ending at $u$ on its penultimate vertex and using induction on the deletion of $u$. After multiplication by $x_u$, terms using $u$ again vanish, and the penultimate in-neighbors complement the out-neighbors that can follow $u$. Applying this to the converse gives $R_u^-=-P(x)R_u^+(-x)$; the whole-vertex coefficient is exactly the stated signed identity. Summing the start/end identity over roots gives $P(x)P(-x)=1$. Thus in every nonempty tournament, the signed two-slot sum, with sign $(-1)^{|J|}$ on the first slot, vanishes. Each size parity contributes $Q(T)/2$.

We also use the [Grinberg–Stanley odd-cycle formula](https://arxiv.org/abs/2307.05569): a packing of $k$ disjoint nontrivial directed odd cycles contributes $2^k$ to $H$, with the empty packing allowed. Add a vertex $z$ beating one old root $u$ and losing to every other vertex. All odd cycles newly containing $z$ traverse an even path in $T$ starting at $u$; the path has at least two vertices and its last vertex necessarily beats $z$. The other cycles lie in the residual $T-V(P)$. The exact outlet count is consequently $H(T)+2E_u(T)$.

## 1. A sharp linear lower bound on positive outlet increments

[THM-4876(7)](../../../01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md) already proves the sole-exit lower bound of $(H+1)/2$ on the positive increment and observed the sharper $(H+3)/2$ through order nine. The result here gives that sharper bound at all orders for $H\ge5$, with an exact degree- and endpoint-dependent inequality. The general asymptotically sharp multiexit bound $(H+1)/2$ was since proved in the [Session214 forced shuffle](../session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md). This sharper additive statement is specifically for a single exit of the appended vertex.

Let $d=|N^+(u)|>0$, $R=T-u$, $W=N^+(u)$, and $D=\sum_{x\in W}t_x(R)$. Then
$$
H(T)\le (2+2/d)E_u(T)-(2+1/d)D.
$$
To see this, $E_u(T)=\sum_{x\in W}O_x(R)$ by deleting the initial $u$. Split every full path on $R$ just before $x$ to obtain an $x$-started second slot and a first slot on the complement. This is injective, so $H(R)\le s_x(R)=2O_x(R)-t_x(R)$. Sum over $W$. Separately, count packings in $T$ by whether a cycle meets $u$. Those avoiding $u$ contribute $H(R)$. Contributions through $u$ use an even path starting at some $x\in W$ in $R$, with the extra condition that it closes back to $u$. Ignoring closure gives at most $2\sum E_x(R)=2(E_u(T)-D)$. Addition proves the estimate.

In particular $H(T)\le4E_u(T)$, and the exact integer strengthening is
$$
E_u(T)\ge2\quad\Longrightarrow\quad H(T)\le4E_u(T)-3.
$$
If $d\ge2$, the first bound already gives $H(T)\le3E_u(T)$; for $E_u(T)=2$ use oddness of $H(T)$. For $d=1$ with sole out-neighbor $x$, the claim follows directly if $t_x(R)\ge1$. If $t_x(R)=0$, then $x$ belongs to a nonterminal strong component of $R$. In fact the vertices which can end a Hamiltonian path are exactly those of the terminal strong component, by its spanning cycle or by singleton triviality. Write the ordinal decomposition around the component $B\ni x$ as $R=A\to B\to C$ with $C\ne\varnothing$. Rooted paths at $x$ ignore the prefix $A$ and split the suffix $C$ between their path and the complement, so
$$
s_x(R)=H(A)s_x(B)Q(C)\ge 2H(R),
$$
using the split-at-root injection and the two empty-slot choices in $Q(C)$. As $t_x(R)=0$, $s_x(R)=2E_u(T)$; again $H(T)\le3E_u(T)$. For $E_u(T)=1$, the first bound and oddness give $H(T)\in\{1,3\}$, both attained. When $E_u(T)=0$, $u$ is a sink and the outlet adds nothing.

For exact unbounded sharpness, let $y,s_1,\ldots,s_r$ be a transitive spine $(r\ge1)$, let $x$ beat only $y$ in this spine, and add $u$ beating only $x$. In $R=T-u$, odd cycles through $x$ include $x,y$ and an odd subset of $s_1,\ldots,s_r$, so $H(R)=1+2^r$. Paths starting $x$ either are singleton, or go through $y$ and any subset of the tail. All complements are transitive and $s_x(R)=1+2^r$. Exactly one full path ends $x$, obtained by following the whole spine and then $x$. Therefore $E_u(T)=O_x(R)=1+2^{r-1}$. Every even path from $x$ in $R$ does close through $u$; thus $H(T)=H(R)+2E_x(R)=1+2^{r+1}=4E_u(T)-3$.

So an equal-source repair using no old-root switch and one added outlet cannot have a small relative step: a positive increment $2E_u(T)$ is at least half the base count, with asymptotically sharp factor. All statements have a converse version with arcs reversed.

## 2. Switching across distinct strong components also costs a linear step

If $u,v$ lie in different strong components of $T$, with the component of $u$ earlier in the ordinal decomposition, then
$$
E_u(T)-E_v(T)\ge H(T).
$$
Write $T=P\to A\to M\to B\to C$ with $u\in A,v\in B$ and $A,B$ the indicated strong components. The blocks $P,M,C$ may be empty, and $Q(\varnothing)=1$. Counts $H,Q$ multiply over ordinal sums. Splitting each suffix into the root path's optional contribution and the residual contribution gives
$$
E_u(T)=\tfrac12 H(P)s_u(A)Q(M)Q(B)Q(C).
$$
Indeed the suffix after $A$ is nonempty, and its two-slot parity cancellation balances the even and odd root terms independently of the prefix.

If $C$ is empty, then $E_v(T)=H(P)H(A)H(M)E_v(B)$. By two-slot parity in $B$, $E_v(B)\le Q(B)/2-H(B)$: the selected even slot must be nonempty, and the empty first slot has weight $H(B)$. Also $s_u(A)\ge H(A)$ and $Q(M)\ge H(M)$. Subtract to get the bound. If $C\ne\varnothing$, then $E_v(T)=H(P)H(A)H(M)s_v(B)Q(C)/2$, while $s_v(B)\le Q(B)/2$. For the latter inequality, a two-cover counted by $s_v$ has $v$ in the first slot, and exactly half the two-covers have this property by swapping slots. The difference is at least $H(P)H(A)H(M)Q(B)Q(C)/4\ge H(T)$, since $Q(X)\ge2H(X)$ for nonempty $X$. Equality is realized in the two-vertex transitive tournament, switching from its source to its sink.

## 3. A changing port inside one component gives adjacent outcomes

Let $W$ be any nonempty tournament and put $h=H(W)$, $q=Q(W)$. Take three singleton vertices $A,x,y$ with $A\to W\to x\to y$ in total-block order, except for the single reversed arc $y\to A$. Call this tournament $S$, and let the common source be $T=S\to s$ where $s$ is one new sink. The reversal makes $S$ strongly connected: traverse $y,A$, a full path of $W$, then $x$, and close to $y$. Both $x,y$ are therefore ports in that same strong component of $T$.

A directed cycle crossing the $W$ block and the exterior, or entirely in the exterior, has to use $y\to A$; otherwise all cross-block arcs follow the total order. Such a cycle takes $y,A$ and a path on an odd subset of $W\to x$. After counting its complementary packings and doubling for the cycle, the contribution is $Q(W\to x)=2q$ by two-slot parity. Packings avoiding the backward arc contribute $h$. Thus $H(S)=h+2q$. The identical argument after deleting $x$ gives $H(S-x)=h+q$.

Appending the sink $s$ gives $E_i(T)=s_i(S)$: each pre-existing rooted path admits one even option between appending $s$ or leaving it to the complement. At $y$, the singleton contributes $h$; every other path begins $y,A$ and optionally uses $W$ and $x$, contributing $2q$. At $x$, the singleton contributes $h+q$, the path $x,y$ contributes $h$, and paths beginning $x,y,A$ with optional $W$ contribute $q$. Thus
$$
H(T)=h+2q,\qquad E_y(T)=h+2q,\qquad E_x(T)=2h+2q .
$$
Add the same vertex $z$ beating only $y$ or beating only $x$, respectively. The two output counts are exactly
$$
3(h+2q) \qquad\hbox{and}\qquad 3(h+2q)+2h.
$$
Both outputs are strong because $S\to s\to z\to S$ and $S$ is strong. In particular for $W=TT_m$, $h=1,q=2^m$: the same physical source of count $1+2^{m+1}$ gives the two adjacent odd outcomes
$$
6\cdot2^m+3,\qquad6\cdot2^m+5\qquad (m\ge1).
$$
The two increments from the common intermediate source are large; only their difference is two. This avoids the single-port lower bound and the cross-component switching barrier. It is an unbounded absolute family with a unit gap, not an iteration giving all targets. With general $W$ the gap is $2H(W)$, so the shape does not give a unit pair on arbitrary reservoirs.

[Literal replay](verify_ports.py) checks the inequalities on all labeled tournaments through order five, the sharp family and same-source transitive family through several widths, and the general block formula on all labeled reservoirs through order four, by subset path DP.

## 4. Whole port palette of a nearly transitive source

Take $S_d$ in natural order $A,r_1,\ldots,r_d,y$, with the sole backwards arc $y\to A$, where $d\ge1$. Its path count is $H=1+2^d$. The all-path root weights are
$$
s_y=1+2^d,\qquad s_A=2^{d+1}-1,\qquad
s_{r_j}=2^{d-j+1}+2^j3^{d-j}-[j=1].
$$
At $y$, either stop or follow $A$ and any subset of the interiors. At $A$, any increasing subset can be used except the singleton subset $\{y\}$, because $y$ beats $A$. For $r_j$, put $a=d-j$ and $b=j-1$. If a starting path omits $y$, it takes an arbitrary subset $X$ of the $a$ later interiors; its complement still contains $A,y$ and has $b+a-|X|$ interiors. The weight is $1+2^{b+a-|X|}$ except when there are no remaining interiors, where the reversed two-vertex edge has count $1$. Summing gives $2^a+2^b3^a-[b=0]$. If the path uses $y$ but not the reversal, there are $2^a$ choices with weight one. If it uses the reversal, each earlier interior is available either after $A$ or to the complement, and each later interior is available before $y$, after $A$, or to the complement. Thus this case contributes $2^b3^a$, again with weight one per choice. The total is the stated palette.

For $k\ge1$ append a transitive sink block $TT_k$ in ordinal order. Its optional subset split contributes $2^{k-1}$ choices to each parity, so the new even weights at the old ports are $2^{k-1}s_i(S_d)$. The different outlets at the old ports now give the full count palette
$$
H+2^k s_i(S_d).
$$
Every output is strong: $S_d$ beats every sink-block vertex, all these beat the new outlet, and the outlet beats a vertex in the strong $S_d$. The two bottom ports give a unit gap for $k=1,d\ge2$, as in §3; at $d=1$ all three ports tie, as they should for a directed triangle.

## 5. Physical binary returns with both signs on the same unbounded count arm

### One-outdegree barrier to the plus sign

Any root $u$ with unique out-neighbor $x$ in $T$ has the exact signed response
$$
2E_u(T)-H(T)=2t_x(T-u)-H(T-u).
$$
Indeed the failed splits counted by $s_u(T)-H(T)$ have a prefix ending at $x$. Their second slot must be just $u$, since any longer $u$-started path uses $x$ again. Their number is $t_x(T-u)$. Also $t_u(T)=H(T-u)-t_x(T-u)$, since any full path ending at old $u$ has last predecessor other than $x$. Subtract these using signed reciprocity.

If $x$ is not a sink of $T-u$, then $t_x(T-u)\le(H(T-u)-1)/2$. This bound is inherited from THM-4876(7): inject paths ending at $x$ into those not ending there by moving $x$ just before its first out-neighbor in the old linear order. The preceding vertex, if it exists, beats $x$; the inverse deletes $x$ and appends it. Parity gives the strict half-bound. Thus an outlet at $u$ reaches at most $2H(T)-1$. If $x$ is a sink of $T-u$, then every other vertex beats both $u$ and $x$, and $H(T)=H(T-u)=E_u(T)$. The outlet instead gives $3H(T)$. In particular no degree-one old root realizes the exact plus gate $2H(T)+1$ for $H(T)>1$.

### Exact classification of the degree-one minus lane

The same signed response gives an exact inverse: a degree-one old root $u$ can give $2H(T)-1$ if and only if $R=T-u$ and its marked vertex $x$ satisfy $H(R)=2t_x(R)+1$. The known [all-order endpoint equality classifier](../../../session213/cycle-arithmetic/EQUALITY.md) gives two possibilities. In the zero case $t_x(R)=0$, $R$ is transitive and $x$ is not its sink. If $r$ vertices follow $x$ in that order, $T$ has $H(T)=1+2^r$: all nontrivial odd cycles use $u,x$ and an odd subset of those $r$ later vertices. In the positive case $R=TT_j\to S_d$ and $x$ is the penultimate vertex of $S_d$. The prefix, of count 1, has no effect on the root weight. At $d=1$, $S_d$ is the directed triangle and $H(T)=5$, already in the first type. At $d\ge2$, the palette and the endpoint equality give
$$
s_x(R)=2+2^d,\qquad t_x(R)=2^{d-1},\qquad H(T)=3+3\cdot2^{d-1}.
$$
Conversely each stated source is physical and the endpoint identity gives the minus gate. Thus the complete set of *predecessor counts* supporting this gate with a degree-one old root is
$$
\{1+2^r:r\ge1\}\ \cup\ \{3+3\cdot2^r:r\ge1\}.
$$
A sink old root is excluded by the degree-one hypothesis. This is an exact thin obstruction for arbitrary target selection in that lane, not a restriction on the full spectrum or on higher-degree gates.

### Two different source identities

For any nonempty $W$ with $h=H(W),q=Q(W)$, there exist two physical sources of *the same* count $H=h+16q$, with outputs $2H-h$ and $2H+h$ from one respective outlet.

For the minus gate, use $S$ of §3 with middle block $W\to TT_3$, omitting the extra sink. In §3 notation the new block count is $h$ and its two-slot count $8q$, so $H(S)=h+16q$. The bottom vertex $y$ has $s_y=H(S)$ and $t_y=h$ (ending at $y$ forces the remaining ordinal path before it). Hence $2E_y=H(S)-h$ and the outlet gives the minus gate.

For the plus gate use the six-vertex base with outgoing bit masks $(10,28,57,48,33,3)$, substitute $W$ at vertex $2$, and add a new vertex beating only old vertex $5$. In the base and augmented quotients the only incoming arc of vertex $2$ is from vertex $1$. Thus a full Hamilton path visits the substituted block at most twice: apart from possibly starting there, every entry requires that unique outside vertex. Count legal quotient words with one and two nonadjacent occurrences of mark 2 and every other vertex once. These counts are respectively $(33,16)$ in the source and $(67,32)$ after the outlet. A one-run word contributes $h$ and a two-run word contributes $Q(W)-2h$, giving source $h+16q$ and output $3h+32q=2H+h$.

The four word counts can be verified with the literal replay or by the following small certificate using ordinary path counts. Before inflation the source has 33 paths and the augmented quotient 67. After replacing vertex 2 by a transitive pair, these rise to 65 and 131. One marked run on that pair has a unique internal order, while two marked runs have exactly two assignments. Therefore the two-run counts are $(65-33)/2=16$ and $(131-67)/2=32$. Both sources used in these transfers are strong. For the plus case the six-vertex source is strong by reachability in the stated bits, and it stays strong under inflating any nonempty block because each old vertex is on a directed cycle using other quotient vertices. Adding the outlet preserves strongness because it beats a vertex of the component and at least one other vertex of that component beats it.

For $W=TT_m$, $m\ge1$, the two constructions therefore certify both binary children $2H-1,2H+1$ for every $H=1+16\cdot2^m$, with $m$ recovered exactly from the target predecessor. The minus child remains in this arm; the plus child does not, so the result is not a repeatable full binary grammar. The sources in the two signs need not be the same tournament.
