# Reciprocal readers and exact changes of strong source

The construction below changes both the source and the reader without moving the absolute count. A second source pair moves the count by exactly two without moving the reader. In both cases, a target can be decoded exactly on the stated physical lane; neither lane claims global coverage.

## 1. The literal source/reader law

For a rooted tournament $(A,u)$, write $h=H(A)$,
$E=\#\{\text{full Hamilton paths ending at }u\}$, and
\[
 P=\sum_{U\ni u}\#\{\text{Hamilton paths in }A[U]\text{ starting at }u\}
                    H(A-U),\qquad H(\varnothing)=1.
\]
Let the reader $B(w)$ have a transitive spine and one apex $v$, where bit $w_i=1$ means $v\to i$. Form a tournament on $A$ and the **converse** of $B(w)$ with all arcs from $A$ to the converse reader except $v\to u$.

The general one-back-edge splice, previously also proved by odd-cycle packing in the [bilateral-glue note](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/BILATERAL_GLUE.md), has a direct positive-path proof. Paths avoiding the backward arc occupy all of $A$, then all of the reader; subtract exactly the paths with the forbidden seam $u\to v$, namely $E$ times the reader's rooted ending count before converse. Paths using $v\to u$ consist of an optional initial run in $A$, a nonempty run in the converse reader ending at $v$, a nonempty run in $A$ starting at $u$, and an optional final run in the converse reader. Choices split independently across the two sides; neither forward seam can be the forbidden edge since its endpoints are used at the backward seam. Thus the total is the product of ordinary counts, plus the product of $P$ values, minus the product of $E$ values.

For nonempty $w$, the reader's $P$ is $1+p$, $p=[w]_2$: after the apex, choosing the first spine contact allows each later vertex to go in its rooted run or its complement run. Its $E$ is $1-e$ for $e=\operatorname{last}(w)$. Therefore
\[
 H_{\rm splice}=hH(w)+P(1+p)-E(1-e).
\]
Write $h(a)$ for the one-apex count of the binary expansion of the positive integer $a$, and $\overline w$ for **reverse-complement**, not bitwise complement. If the source has integral positive ratio $P/h=a$, replace the reader word by an auxiliary completed word $W=w\,\overline{\operatorname{bin}(a)}$. Exact inversion weights across the concatenation seam give
\[
 H(W)=H(w)+h(a)-1+ap+(a\bmod 2)e.
\]
The product $ap$ accounts for all cross inversions, and $(a\bmod 2)e$ adds the repeated weight for adjacency at the seam. Hence
\[
 H_{\rm splice}=
 hH(W)+h\bigl(1-h(a)+a-(a\bmod2)e\bigr)-E(1-e).
\]
This is a physical affine selector, not subtraction of counts from several different tournaments. It also specifies a limitation to this strategy: changing among a fixed finite bank of these integral-ratio sources gives only affine transforms of slices of the one-apex spectrum, not an independently tunable integer offset. It does not exclude covering the full spectrum through a sufficiently diverse union.

Converse preserves $H(W)$ but can change both the word and the source required to remove a suffix. For two sources of equal $h$ with ratios $a,b$, choose $W$ beginning with $\operatorname{bin}(b)$ and ending with $\overline{\operatorname{bin}(a)}$. The source of ratio $a$ reads $W$ minus that suffix; the source of ratio $b$ reads $\overline W$ minus its $\overline{\operatorname{bin}(b)}$ suffix. The completed one-apex baseline is preserved exactly. Match the two displayed affine offsets to obtain an exact return. The special case $a=b=1$ and $W$ beginning **11**, ending **00** gives the common count $h(H(W)+1)-E$ whenever the sources' $E$ values match. If the $E$ values differ, the change is exactly the negative $E$ difference. That contrast costs no adjustment at all in the completed word.

## 2. Actual strong sources, not formal profiles

With vertex labels $0,\ldots,6$ and marked vertex 6, take outgoing row masks
\[
 A=(42,124,121,112,97,64,1),\quad
 B=(40,109,121,112,99,64,1).
\]
Both are strong, with $(h,P,E)=(31,31,9)$, but their sorted outdegree lists are $(1,1,3,3,3,5,5)$ and $(1,1,2,3,4,5,5)$. The two side assignments above therefore give a genuine change of source, and of reader unless $W=\overline W$. Every resulting graph is strong: the marked source reaches its whole side, and every vertex in the converse reader reaches its marked apex (the original reader begins with 1). The back edge and the forward edges between the two sides complete the connections.

The two returns below have literal subset-DP checks on **both** outputs, and their full-graph degree multisets distinguish the outputs even up to converse.

| prescribed absolute target | $W$ | $\overline W$ | order of each output |
|---|---|---|---:|
| $2719$, prime | 1100100 | 1101100 | 14 |
| $92743=7\cdot 13249$, prime cofactor | 110000010000 | 111101111100 | 19 |

Neither output is inferred merely from a congruence or difference.

## 3. Selectable endpoint and an exact two-unit repair

Use the strong rooted sources
\[
 C=(112,97,99,39,110,64,8),\qquad
 D=(248,97,227,230,238,128,160,2)
\]
with the last vertex marked. These are from the ready construction with a source core followed by a transitive tail and a final marked vertex beating only the core root. Literal profiles are $(59,59,9)$ and $(59,59,11)$. Keep the same reader $w=W$ minus its final zero. Now every framed $W$ gives the absolute pair
\[
 \bigl(59(H(W)+1)-9,\;59(H(W)+1)-11\bigr).
\]
The repair is exactly $-2$ for the same reader (or for $\overline W$ with its last zero removed). For prescribed $11083$, the baseline is fixed, $H(W)=187$, found with $W=$ **11001000**. Both literal outputs are strong and have counts $11083$ (prime) and $11081=7\cdot 1583$ (prime cofactor), checked by subset DP. Changing source gives the target offset while preserving the same physical reader. The distinct orders also prevent isomorphism.

Structural source-search rule: if $v$ beats exactly one vertex $u$ and the remaining tournament is $Q$, then
\[
 H(T)=H(Q)+P(Q,u)-E(Q,u),\quad
 E(T,v)=H(Q)-E(Q,u),\quad P(T,v)=H(Q)+P(Q,u)=H(T)+E(Q,u).
\]
The first identity is the single-back-edge splice with one side a singleton; the second restricts full paths of $Q$ to endpoints other than $u$. For the third, the rooted run is either the singleton $v$, or it starts $v,u$ and then uses a subset of $Q$ starting at $u$; its complementary path count is just the untouched factor. Thus choosing the source endpoint data is equivalent to looking for predecessor rooted $Q$ with a chosen sum $P(Q,u)+[H(Q)-E(Q,u)]$ and controlled term $H(Q)-E(Q,u)$, not guessing the count of the extension. In particular if $Q=C_0\to D_0$ with $D_0$ nonempty and $u$ in $C_0$, then $E(Q,u)=0$ and $P(T,v)=H(T)$, $E(T,v)=H(C_0)H(D_0)$. The independent factorization $P(Q,u)=P(C_0,u) F_2(D_0)$ uses the two labelled path slots (empty slots allowed); this is the search rule that produced the two count-59 sources.

## 4. Target-first decoding and an unbounded limitation

**Scope of endpoint motion.** For a strong non-singleton source with \(P=H\), the [unsigned equality classification](https://github.com/eliottcassidy2000/math/blob/2d70eb9cd31437bbfa8cf418c6285f79fc8afaaa/05-knowledge/results/chanpak_11_13_26_ready_barrier_20261114/READY_CLASSIFICATION.md) forces outdegree one at the root. The preceding deletion identities force \(E(Q,u)=0\); strength forces \(u\) into the first SCC \(C_0\) of \(Q\), so the later part \(D_0\) is nonempty. Every Hamilton path of \(D_0\) supplies \(|D_0|+1\) distinct ordered two-slot covers (all cuts including the two empty-slot placements), and different paths supply disjoint sets of these covers because concatenation recovers the path. Since \(P(C_0,u)\ge H(C_0)\) by splitting each Hamilton path just before \(u\), we get
\[
 H=E+P(C_0,u)F_2(D_0)\ge (|D_0|+2)E.
\]
In particular \(E\le H/3\), with equality precisely at a singleton \(D_0\) and \(P(C_0,u)=H(C_0)\). This is sharp, starting at the directed triangle and iterating the singleton-tail construction. Thus offsets from a fixed equality source cannot range over an entire congruence system modulo \(H\); one must change the scale as well as the offset.

For given odd $N$ on the count-31 lane, solve $M=(N+9)/31-1$, rejecting a noninteger. On either count-59 branch use $M=(N+E)/59-1$ with $E=9$ or $11$ as prescribed, then decode the same framed baseline. Use the [half inverse](https://github.com/eliottcassidy2000/math/blob/c598c9ba1c079e42d82c61347185b995b03c84cc/05-knowledge/results/chanpak_long_cycle_and_transfer_20261114/HALF_INVERSE.md) on $H(W)=M$, restricting both half integers $x,y$ to start with binary **11**. This is exactly equivalent to the required first two and last two bits of $W$. The monotonic search on each balanced cut stays valid on these narrower intervals, so the decoder constructs the physical outputs or rules out this particular route at **all orders**, in $O(\sqrt M\log^2 M)$ arithmetic operations for positive $M$.

The doubly framed restriction makes a visible infinite exclusion. At each framed width $L$, the half-product identity gives $M\ge xy+1\ge9\cdot 2^{L-4}+1$. The bound $h(x)\le x+1-(x\bmod2)$, and the same bound for $y$, yield $M\le (x+1)(y+1)\le 2^L$. Consecutive framed widths have a multiplicative gap: no framed count in $(2^L,(9/8)2^L+1)$, $L\ge4$. It transfers to genuine gaps inside the correct $N\bmod 31$ residue. To reach all targets one must move beyond fixed framed readers or build a complementary independently changing source bank, not just swap the two existing witnesses more often.

For the fixed-reader repair only, the initial **11** is unnecessary: the physical reader must start with 1, and the completed word must end in **00**. In the balanced inverse, drop the initial **11** condition on \(x\), keeping it on \(y\); width three has the single word **100**. In particular the prescribed prime \(1171=59(19+1)-9\) supplies the physical word **10100** with baseline 19. The same reader gives \(1169=7\cdot167\) from the other source. The baseline 19 is absent from the doubly framed family at every width, so refusing the unnecessary frame would lose this actual pair. The multiplicative holes above are proved for the *doubly framed* decoder, not the widened decoder.
