# Exact integer windows from independent cut switches

This is a conditional construction, not a solution of the physical palette problem.
The contiguous-box specialization overlaps [CRT_AMPLIFIER.md](CRT_AMPLIFIER.md);
the residue-extremal form and physical count-three/five results below are separate.
It uses the [independent-switch law](../chanpak_selective_target_20261114/SWITCH.md):
build a tournament by putting modules $B_i$ in ordinal order and extending
by one apex, with a separate mark $U_i\subseteq B_i$ in every module. Let
$h_i=H(B_i)$, $Q_i=Q_2(B_i)$. For a nonempty directed path in a module, use
its complement Hamilton count as weight. Let $E_s,O_s$ be the even/odd
weights of paths starting in $U_i$, and $E_e,O_e$ the corresponding end
weights for paths ending outside $U_i$. The internal cut coefficient $c_i$
counts weighted even paths with both constraints. Holding
$(h_i,Q_i,E_s,O_s,E_e,O_e)$ fixed while changing $c_i$ keeps every
crossing cycle contribution fixed. In the odd-cycle packing expansion the
only changed terms are internal apex cycles, multiplied by the Hamilton
counts of all other modules. Hence for $P=\prod_i h_i$, the exact
count satisfies

$$
H=H_0+2\sum_i d_i P/h_i,\qquad d_i=c_i-c_i^0.
$$

Here $H_0$ is the absolute Hamilton count of the all-baseline construction,
including all crossings, and is *not* a free parameter.

## The carried-window theorem

Assume $k\ge2$, $h_i>1$ pairwise coprime, and put $w_i=P/h_i$.
Suppose module $i$ has a physical response-matched palette of internal
digits $0,\ldots,m_i h_i-1$, with $m_i\ge1$. Put

$$
K=\sum_i(m_i-1),\qquad F=(k-1)P-\sum_i w_i.
$$

Then the composite achieves $H_0+2t$ for **every** integer

$$
F+1\le t\le(K+1)P-1.
$$

The immediately neighboring half-gaps $F,(K+1)P$ are absent in this
palette product. With $K\ge k-1$ this interval has at least $P$
positions and includes an entire absolute period (length $2P$ in
ordinary count). In particular it suffices to double all but one local
palette; if every palette is doubled then one obtains $2P+\sum_i w_i-1$
consecutive odd counts.

Proof and decoder. For prescribed $t$, compute

$$
u_i\equiv t w_i^{-1}\pmod{h_i},\qquad 0\le u_i<h_i .
$$

Other weights vanish modulo $h_i$, so these digits are forced modulo
$h_i$. Their weighted sum $s=\sum_i w_i u_i$ satisfies
$s\equiv t\pmod P$ and ranges from $0$ to $kP-\sum_i w_i$.
The extra choices are $d_i=u_i+h_i v_i$, with every
$0\le v_i\le m_i-1$ available, and they contribute exactly
$P\sum_i v_i$. Therefore $t$ is attainable iff
$0\le(t-s)/P\le K$. If $t>s-P$ and $t<s+(K+1)P$, both forced
inequalities hold. The advertised range guarantees these strict bounds
for **all** $s$. The lower adjacent failure is the maximal canonical
$s$ minus $P$; the upper failure is $s=0$ plus $(K+1)P$.
Distribute the exact target carry greedily up to the individual capacities,
then choose the corresponding literal tournament modules. This proof
reads like the conductor calculation in a numerical semigroup generated
by the cofactors $w_i$, except that the upper boundary is equally
essential: finite physical palettes give bounded carries.

**More general digits.** Arbitrary physical palettes may be far smaller
than long contiguous intervals. If each residue fiber in every module consists of
a consecutive range of lifts $u_i+h_i[L_i(u_i),U_i(u_i)]$, then for a
residue $r\bmod P$ the attainable quotients are exactly

$$
q(r)+\left[\sum_i L_i(u_i),\sum_i U_i(u_i)\right],\qquad
q(r)=\frac{\sum_i w_i u_i-r}{P}.
$$

The absolute bounds separate over the blocks: set

$$
A=\sum_i w_i\max_{u_i}\bigl[u_i+h_i L_i(u_i)\bigr],\qquad
B=\sum_i w_i\min_{u_i}\bigl[u_i+h_i U_i(u_i)\bigr].
$$

Then $[A-P+1,B+P-1]$ is wholly attainable and both adjacent half-gaps
are absent. Indeed the largest lower missing value in a residue is its
least attainable half-gap minus $P$, and the smallest upper missing
value is its greatest attainable half-gap plus $P$. Maximizing and
minimizing these sums separates because the canonical local residues
vary independently by CRT. This need not be the longest unrelated run
elsewhere. An aligned absolute period $[zP,(z+1)P-1]$ exists exactly
when $\lfloor A/P\rfloor\le z\le\lfloor B/P\rfloor$. A missing residue in any local palette precludes a full
period, and a common residue system alone does not guarantee one. To
overlap windows from distinct core tuples one must compare their *actual*
$H_0$'s and endpoints. This theorem isolates the carry condition without
purporting to place the baseline.

**Obstruction for a finite catalogue.** With noncoprime counts, the gcd
of switch weights is $P/\operatorname{lcm}(h_1,\ldots,h_k)$, whether or
not digits are consecutive. Thus a repeated finite catalogue cannot create
unit-spaced windows by taking more independent switches in one composite.
Its modulus grows with $P$; adding any number of transitive tape
vertices does not alter it. This is an obstruction to the exact additive
mechanism, not to all ways of changing cores.

## A physical obstruction at count three, even with arbitrary padding

If $H(B)=3$, then $B=TT_a\to C_3\to TT_b$. Indeed Hamilton counts
multiply over the ordinal SCC decomposition, a cyclic SCC has count at
least 3 by the positive odd-packing formula, and [Busch's strong bound](https://doi.org/10.37236/1141)
forces a strong count-three SCC to have order at most three. Let $x,y$
be the mark words, as binary numbers of lengths $a,b$, and $j\in
\{0,1,2,3\}$ the number of marked triangle vertices. With
$S=E_s+O_s,T=E_e+O_e$, the full two-slot count and start/end reads are

$$
\begin{aligned}
Q&=12\cdot2^{a+b},\\
S/3&=2^b(4x+j)+y,\\
T/3&=2^a(4\operatorname{rc}_b(y)+3-j)+\operatorname{rc}_a(x).
\end{aligned}
$$

Here $\operatorname{rc}_\ell$ means bit reversal and complement at
fixed width $\ell$, including any leading zeros. The formulas follow
from factorization of the two-path count over ordinal blocks, and
by rotation each triangle vertex contributes start weight 3 and end
weight 3. A selected vertex contributes its one slot within $C_3$
times the Hamilton count of every earlier block and the two-slot count
of every later block; end weight factors oppositely.

Interpret $S/3$ as a binary word $w=x\,\Vert\,\operatorname{bin}_2(j)
\,\Vert\,y$ of length $a+b+2$. The observable discrepancy is

$$
T/3-\operatorname{rc}_{a+b+2}(S/3)
=2^a(\operatorname{rev}_2(j)-j).
$$

If nonzero, its sign determines $j$ and its absolute value determines
$a$, thus also $b,x,y$. A fixed value of $j$ determines a marked
triangle up to rotation. If the discrepancy is zero, the entire triangle
is monochrome. Replacing it by two transitive same-color vertices does
not change the normalized transfer: for the monochrome triangle
$(h,Q,n)=(3,12,3)$, with total start and signed start $(9,3)$ if
all marked, zero otherwise, and total end and signed end $(9,3)$ if
all unmarked, zero otherwise. For two transitive same-color vertices
these are all divided by 3. By the
[cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md),
the whole marked base gives $H(B+\mathrm{apex})=3H(TT_w+\mathrm{apex})$,
regardless of the hidden placement of that triangle. This also follows
by splitting the apex cycle according to which block supplies its
endpoints.

In both cases the unpadded apex count $h+2c$, hence $c$, is
**determined** by $(Q,S,T)$ alone. Consequently every profile-matched
physical palette at $h=3$ has a single value of $c$, even across
unbounded transitive padding. The unrestricted transitive case $h=1$
has the same uniqueness because $S$ directly records the complete
mark word and $Q=2^{|B|}$ records its width. Neither the smallest cyclic
multiplier nor multiplier one is a source of independently matched
carry digits. This only concerns independent ordinal-apex switches,
not other constructions with these counts.

Replay: `python replay_crt_windows.py` verifies exact conductor endpoints,
the target decoder, and small triangle-padding controls by literal
subset-endpoint counts. The threshold and count-three theorem are
proved above at every height.

## A literal unit switch at the next multiplier, five

There are two distinct ways to position the sole strong four-vertex
count-five core among three transitive vertices. Use these two bases and
marks (outgoing masks use bit $i$ for vertex $i$):

$$
\begin{array}{c|c|c|c}
\text{rows}&U\text{ mask}&H(B+\text{apex})&(h,Q,S,T)\\
(118,124,120,113,96,64,0)&50&69&(5,256,102,149)\\
(126,124,120,48,96,64,8)&38&67&(5,256,102,149).
\end{array}
$$

The total start/end responses match, but their signed parts need not.
There is a general neutral repair: prepend one **unmarked** singleton,
append one **marked** singleton. They are a source and sink, respectively,
even with the apex included, so do not change $H(B+\mathrm{apex})$.
In the base the left singleton doubles the old total end weight and adds
$h$ from its singleton path; it leaves the total start weight unchanged.
The right singleton doubles the old total start weight and adds $h$
likewise, leaving the end weight unchanged. Start and end signed parts
both equal $h$: the signed count of marked starts is the number of full
paths ending in a mark, which now all end at the marked sink; the
converse statement holds for unmarked ends because the source is
unmarked. Equivalently, the internal signed terms cancel via the
parity-changing optional vertex.

Thus the repaired pair has equal complete response
$(h,Q,E_s,O_s,E_e,O_e)=(5,1024,102,107,149,154)$ and local
$c$ choices $31,32$. In every ordinal context this switch gives
an exact $2H_{\rm left}H_{\rm right}$ shift, not merely an equality
modulo a proposed target period. **This is one physical carry digit,
not a residue-complete palette modulo five.** It contrasts with the
all-order impossibility at multiplier three.
