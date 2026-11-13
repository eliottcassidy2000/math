# Exact adjacent corridors with an unbounded changing reservoir

Write $h=H$, $q=Q(2)$ for two ordered slots allowing empties. In a block with a cut $U$ (the apex's out-neighbors), $s,e$ are full paths starting in $U$ and ending outside $U$, respectively; $S,E$ are their optional-other-slot analogues; $c$ is the count after adding that apex. The starts and ends here belong to the marked path, which is required nonempty.

## Serial propagation of the exposed end

For two ordinal blocks $P\to Z$ with one common apex, the joined terminal parameters are
$$
 h'=h_P h_Z,\qquad s'=s_P h_Z,\qquad E'=E_P h_Z+q_P E_Z,
$$
$$
 c'=h_Z c_P+h_P(c_Z-s_Z)+e_P(s_Z-h_Z)+E_ZS_P.
$$
Indeed the first and last local terms are $(s_P+b_P)h_Z$ and $(e_Z+b_Z)h_P$, where $b$ counts two nonempty sides within a block. A forward seam contributes $e_Ps_Z$, and a backward seam $E_ZS_P$. These partition all paths including endpoint cases, giving the formula after $c=s+e+b$. The expression for $E'$ splits by the block in which the marked path ends. This is the old common-apex law; the new use is to preserve everything it sees except a pure exposed-end change.

Two literal prefixes give such a change:

| outgoing masks on 7 vertices | out-cut mask | $h,q,s,S,E,c$ | $e$ |
|---|---:|---|---:|
| $(16,1,3,7,14,31,63)$ | 91 | $9,336,9,227,82,73$ | 2 |
| $(0,33,3,7,15,29,63)$ | 101 | $9,336,9,227,82,73$ | 0 |

Bits label vertices from zero. The replay computes these by subset DP. Thus for **any** nonempty suffix $Z$ with a physical cut, the two reservoirs have identical $(h',s',E')$ and their new anchors differ by exactly
$$
c'_{e=0}-c'_{e=2}=2(h_Z-s_Z).
$$
The difference is supplied by exposed-end interference, not by independent local offset digits. If $Z$ is strong and its in-cut nonempty, every vertex of the full reservoir $P\to Z$ can reach the apex. A nonempty binary prefix beginning in 1 therefore makes both full output tournaments strong. The [terminal law from Session 317](../session317_exact_anchor_atlas_20261113/TERMINAL_ATLAS.md) then changes both counts by the same amount for **any** choice of that binary prefix, not just for a prescribed leading bit.

## Exact digit control of the gap

Use the strong suffix $Z=w\to a\to TT_m\to w$ with $m\ge1$ and ordered chain $b_1,\ldots,b_m$. Then $h_Z=2^m+1$. The Hamilton starts at $w$ and $a$ each have count one, and the start at $b_j$ has exactly $2^{m-j}$ paths. In the latter case the path must reach $w$, then $a$, and resume in the chain, possibly with no resume. Each later chain label may independently go before or after $w$; earlier chain labels go after $a$.

Hence choose as the in-cut any nonempty subset $J$ of the $b_j$. The missing start mass is exactly $r=\sum_{j\in J}2^{m-j}$, and every $1\le r<2^m$ occurs. The prefixes above yield strong outputs with exact gap $2r$, all with a **matched changing reservoir count** $h'=9(2^m+1)$. For any fixed $r$, increase $m$ without bound and recode its binary digits. This differs from merely producing every even gap: the whole terminal response $(h',s',E')$ is matched, and survives any common binary-prefix construction.

A particularly transparent adjacent corridor uses the in-cut $\{a\}$ instead of the chain. Its missing start mass is one, its end response equals $h_Z$, and $c_Z=2h_Z-1$. To verify the latter two, in any tournament $R$ add a vertex $v$ with sole in-neighbor $u$: its optional-slot end response is $H(R)+E_R(u)$ and the new full count is this minus the old start count at $u$. Apply this first to $a$ with $Z-a=TT_m\to w$, then to the common apex with sole in-neighbor $a$. For $t=2^m$, the resulting two strong path counts after any binary prefix $B\ge1$ are
$$
(t+1)[309+418B+9D(B)]-e,
\qquad e\in\{0,2\},\qquad D(B)=B-K(B)-(B\bmod2).
$$
Here $K$ is the terminal-atlas binary defect. Its recurrence gives $D(B)\ge0$ by induction. At $B=1$ the $e=2$ chart gives $727(2^m+1)-2$, divisible by seven when $m\equiv2\pmod3$, an infinite explicit strong $7$-divisible ray with a strong neighbor at distance two. No primality of the cofactor is claimed.

## Exact mixed-power weights and a density obstruction

For the chain cuts $J$ above, the optional-slot end weight of $b_j$ is
$$
 E_j=2\cdot3^{j-1}2^{m-j}+2^j-[j=m],\qquad E_Z=\sum_{j\in J}E_j.
$$
To count it, let the first slot be the path ending at $b_j$. When both singleton labels $w,a$ are in it, its word is $[B],w,a,B$ with the leading run optional. The $j-1$ earlier chain labels have three destinations, the later labels two. This gives $3^{j-1}2^{m-j}$. When only $a$ lies in the marked slot, there are $2^{j-1}$ choices. Only $w$ cannot occur there. When neither lies there, the second slot has a cyclic join on its remaining chain labels and $w,a$, with count $1+2^k$ at remaining chain size $k>0$ and 1 at $k=0$. Summing over the free earlier labels gives $2^{j-1}+3^{j-1}2^{m-j}-[j=m]$. This proves the formula and exposes literal mixed powers, not an approximate rational-base analogy.

These exact gap controls do not imply absolute coverage. In the whole chain-cut family with both literal prefixes and every binary $B\ge1$, each terminal offers at most $X/E'$ prefix indices with count at most $X$ (use $K(B)\le B$ and $c'\ge s'$). If $j=\max J$, then $E_Z\ge 2^m(3/2)^{j-1}$, while there are $2^{j-1}$ sets with that maximum. Since $E'\ge336 E_Z$,
$$
\#\{H\le X\}\le\frac{2X}{336}\sum_{m\ge1}2^{-m}\sum_{j=1}^m(4/3)^{j-1}=X/56.
$$
The double sum equals 3, by reversing the order. Allowing $B=0$ adds only $O(X^\alpha)$ for any $\alpha>\log_3 2$: the terminal count is at least $227 E_Z$ from the join formula, and $\sum_{m,j}2^{j-1}(2^m(3/2)^{j-1})^{-\alpha}$ converges for such $\alpha$. The inequality is a genuine absolute-density bound on this chain-cut grammar, not on all possible suffixes accepted by the pure-end pair. The limitation locates the missing step: change the end-weight system or find a much wider matched endpoint menu, not merely more exact even gaps.

## Independent cancellation example

A second mechanism uses a single physical eight-vertex prefix $(4,65,50,7,11,27,61,127)$ with cuts 201 and 177. Both have $(h,q,s,E)=(91,2392,91,1058)$, while $(e,S,c)$ is $(36,1679,921)$ or $(39,1682,919)$. Here none of the three changing coordinates is inert. Start with $R$ equal to $w\to a\to TT_m\to w$ as above, and add $v$ with sole in-neighbor $w$. For the suffix cut omitting $v$ alone, the same sole-in-neighbor identity gives $(h_Z,s_Z,E_Z,c_Z)=(3t-1,2t-1,3t,5t-1)$ with $t=2^m$. Indeed $E_R(w)=2t-1$ by splitting the chain before the predecessor of $w$; the start at $w$ is one; and the start at $v$ is $H(R)-1=t$. The serial law gives $(h',s',E')=(273t-91,273t-91,10350t-1058)$ and anchors $8037t-921$, $8037t-919$. The large drift cancels exactly at every $m$. This is a second adjacent corridor with changing count, but no overlap across all prescribed targets is claimed.
