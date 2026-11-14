# A two-block source and a changing-reader triple

For a tournament $T$, write $h_T=H(T)$ and let the ordered two-path-cover count be
$Q_T=\sum_{X\subseteq V(T)}H(T[X])H(T-X)$, with $H(\varnothing)=1$.
For nonempty $A,B$, build $S=K(A,B)$ as the cyclic substitution
$r\to A\to B\to r$. At root $r$, let $E_r,O_r$ be the weighted even/odd sums
over paths starting at $r$, weighted by the Hamilton count of the remaining
tournament. Put $u=h_Ah_B$ and $d_A=Q_A-h_A$. Then
$$
(H(S),E_r,O_r)=\left(2u+d_A(Q_B-h_B),\ \tfrac12d_AQ_B,\ u+\tfrac12d_AQ_B\right).
$$

**Proof.** The signed two-cover sum
$\sum_X(-1)^{|X|}H(T[X])H(T-X)$ vanishes for nonempty $T$.
In the inherited odd-cycle packing expansion (each cycle contributes a
factor 2), each cycle or leftover singleton independently chooses one of
the two sides and has odd cardinality. If $r$ is uncovered in the packing
of $S$, its contribution is $u$. Every packing with $r$ on a cycle uses a
directed path on a nonempty subset of each of $A,B$; it cannot visit
either block twice without revisiting $r$. The odd-cycle indicator
multiplied by the cycle weight 2 is $1+(-1)^{|P_A|+|P_B|}$.
The unsigned contribution is $d_A(Q_B-h_B)$ and the signed contribution
is $(-h_A)(-h_B)=u$, since the empty subsets were removed.
A rooted path other than $r$ visits a nonempty subset of $A$ and an
arbitrary subset of $B$. Its complement is an ordinal join of what
remains of the two blocks, so its total weight summed over paths is
$d_AQ_B$. The signed part vanishes by the signed sum for $B$.
The singleton rooted path adds $u$ to the odd coordinate.

## Physical affine-product reader

For a tournament $R$ with entry $w$, let $a=H(R)$, let $q_e,q_o$ be
its weighted even/odd path sums for paths **ending** at $w$, and put
$s=a+q_e+q_o$, $c=2q_e$. Take the ordinal join $S\to R$ and reverse
just the arc $r\to w$. The resulting count is
$$
F_{a,s,c}(A,B)
 =(2a+c)h_Ah_B+(Q_A-h_A)(sQ_B-a h_B). \tag{1}
$$
Indeed, packings not using the seam contribute $aH(S)$. An odd cycle
using the lone backward arc comprises a path starting at $r$ in $S$
and a path ending at $w$ in $R$. They have opposite vertex parities,
so their forward connection is never the deleted arc $r\to w$.
The complement counts factor and the cycle has weight 2. This proves
the inherited seam row $(a,2q_o,2q_e)$ directly. The source formula then
gives (1). This is a formula for two jointly realized source coordinates,
not a relaxation to arbitrary integer points.

## Three consecutive counts at two independent heights

The reader row masks use zero-based vertices.

| Reader | Entry | a | q_e | q_o | s |
|---|---:|---:|---:|---:|---:|
| (26,60,57,48,32,1) | 2 | 15 | 7 | 12 | 34 |
| (14,20,24,18,1) | 1 | 15 | 8 | 11 | 34 |
| (26,60,57,48,32,1) | 0 | 15 | 9 | 10 | 34 |

Literal subset DP verifies the table. These reader choices keep the
entire product term of (1) fixed and change only the intercept by
$2h_Ah_B$. In particular, $A=TT_k,B=TT_m$, with $k,m\ge1$ independent,
gives the infinite exact triple
$$
(2^k-1)(34\cdot2^m-15)+\{44,46,48\}.
$$
The first triple is 97,99,101. The outer rows use different entries of
the same actual reader; the middle row genuinely changes the reader.

The readers $C_3$ (masks (2,4,1), entry 0) and (4,13,8,1)
(entry 1, the universal source over $C_3$) give $(a,s,c)=(3,6,2)$
and $(3,6,0)$. They yield $(2^k-1)(6\cdot2^m-3)+\{6,8\}$
at both unbounded heights, or gap $2h_Ah_B$ for unrestricted $A,B$.

The earlier slope-72 67/69 reader pair instead has effective row
difference $(2,-4,0)$. Thus reading $K(A,B)$ gives difference
$2(3h_A-Q_A)h_B$. From $A=TT_1$ to $TT_2$ the polarity changes from
$+2h_B$ to $-2h_B$; at $TT_2,TT_m$ the respective counts are
$216\cdot2^m+5$ and $216\cdot2^m+3$. These are the only possible
first blocks with $|3h_A-Q_A|=1$. Indeed, cutting each Hamilton path
at any of its $|A|+1$ gaps injects to the ordered covers, so
$Q_A\ge(|A|+1)h_A$. For $|A|\ge3$ this implies
$Q_A-3h_A\ge h_A$, and if $h_A=1$ then $A$ is transitive by the
odd-cycle formula, giving $Q_A\ge8$. The changing readers bypass
this all-order width-two restriction of the old readers.

## Exact inverse and prime-target controls

A finite reader list with two transitive source slots gives only
$O((1+\log X)^2)$ outputs up to $X$, since the two heights are
individually $O(\log X)$. For the triple, a target $N$ requires
$z\in\{44,46,48\}$ and $k\ge1$ with $2^k-1\mid N-z$ and
$((N-z)/(2^k-1)+15)/34$ an integral power of two of exponent at
least 1. This is sufficient as well as necessary inside the family.
With unrestricted $B$, equation (1) instead requires the **same
actual tournament** to realize both coordinates $(h_B,Q_B)$.

In the transitive triple, branch 44 is divisible by 11 exactly if
$10\mid k$ or $m\equiv2\pmod {10}$, since 2 has order 10 modulo 11
and $34\cdot2^m-15\equiv2^m-4$. Branch 46 is divisible by 23
exactly if $11\mid k$ or $m\equiv4\pmod {11}$, since 2 has order 11
and $15/34\equiv16\pmod {23}$. Branch 48 is divisible by 3 exactly
if $k$ is even. Thus the **whole triple is composite** whenever
$110\mid k$, for every $m$. A local exact step is mathematically
separate from prime placement.

With unrestricted blocks, every output in (1) is divisible separately
by $\gcd(h_A,Q_A)$ and by $\gcd(h_B,Q_B)$. For transitive blocks and
a finite reader list, infinitely many primes and infinitely many $7p$
targets are missed unconditionally: a set with $O(\log^2 X)$ count
up to $X$ has convergent reciprocal sum by dyadic grouping, whereas
both prime and $7p$ reciprocal sums diverge. This is restricted to
that family, not Chan–Pak.

## Source-swap determinant and an exact nontransitive compensation

Swapping $A,B$ keeps $H(S)$ and $O_r-E_r$ fixed, but moves each of
$E_r,O_r$ by $(Q_Ah_B-Q_Bh_A)/2$ (first orientation minus second).
In the spirit of the older [Giuga block-lattice notes](../giuga_block_lattice_20261014f.md),
the determinant sees the two joint coordinates and the orientation,
not just a scalar fingerprint; no transfer of a realizability theorem
from that lattice is asserted.

If both blocks have order at least 2 then $4\mid Q_A,Q_B$. Indeed,
an odd-cycle packing term of the cover polynomial at 2 has size
$2^{b+j}$, where $b$ counts packing blocks and $j$ counts cycles.
If $j=0$ then $b\ge2$; otherwise $b+j\ge2$. Thus a same-reader swap
can move the output only by a multiple of four. Using different
readers in one $(a,s)$-fiber removes the parity obstruction:
$$
F_{a,s,c}(A,B)-F_{a,s,c'}(B,A)
 =(s-a)(Q_Ah_B-Q_Bh_A)+(c-c')h_Ah_B. \tag{2}
$$

For a literal witness take $A=(16,17,3,23,4)$ and
$B=(56,57,35,4,44,8)$, with $(h_A,Q_A)=(5,64)$ and
$(h_B,Q_B)=(19,256)$. The determinant is $-64$, and $h_Ah_B=95$.
Use the $C_3$ reader on $K(A,B)$ and the singleton-source-over-$C_3$
reader on $K(B,A)$. Then (2) is $3(-64)+2\cdot95=-2$, with actual
counts **88021 and 88023**, independently recounted. Both source
blocks are nontransitive. Using the same reader in both orientations
would instead give gap $-192$. The first output is strong; in the
second the reader's $C_3$ is a terminal component. Both counts are
composite; this witness demonstrates compensation, not prime placement.

Since $c,c'\in[0,2(s-a)]$, an exact gap $\pm2$ in (2) requires
$$
\left|Q_A/h_A-Q_B/h_B\right|\le 2+\frac{2}{(s-a)h_Ah_B}.
$$
Since $Q_B/h_B\ge|B|+1$, fixing $A$ bounds the order of $B$
uniformly **even with arbitrary growth and changes in the reader pair**,
provided it shares $(a,s)$ within each comparison. This does not
bound reader order or absolute output counts. In the concrete lane
above, the equation for gap $-2$ reads $Q_B=(202h_B+2)/15$. Integrality
forces $h_B=19+30t$ and $Q_B=256+404t$, $t\ge0$, so
$Q_B/h_B\le256/19$ and $|B|\le12$. An unbounded source-order argument
along this lane must change both source blocks or leave the
equal-$(a,s)$ setup.

## Checks

Run replay_source_switch.py (also works in optimized mode). Literal subset DP
checks 121 pairs of nonempty blocks and their swaps, five readers,
60 further glued counts and the two nontransitive compensation counts.
The exact scalar inverse is checked through 3999.
