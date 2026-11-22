# Exact +2 on independent physical branches: a change of unit line

**PROVED:** two physical branching flanks can be grown independently, and an eight-vertex middle switch changes their ordinary Hamiltonian-path count by exactly two for every pair of histories. The absolute image has a complete inverse by positive divisor tests and the separate history inverses. Three charts give infinitely many actual counts in every odd class modulo 256; their union still omits classes modulo 512. A specific twin-prime pair has exact certificates. Chan–Pak remains open.

## 1. A second unit line with two independent moves

Use the [rooted migration](../chanpak_bank_fusion_20261117/migration/PROOF.md) starting from the transitive pair rooted at its last vertex. Its state $(h,R)$ is the ordinary Hamilton count and the number of choices of a marked-terminal path and a complementary path. Exactly one full path ends at the terminal, and none starts there. A leading-zero address retains these properties, with transfer
$
M(k,u)=\begin{pmatrix}F(u)&u+1\\ 2^k-\mathrm{rev}_k(u)&2^k\end{pmatrix},
\quad
F(0)=1,\quad F(2v+1)=F(v),\quad F(2v)=F(v)+v+[v\text{ odd}]\ (v>0).
$
Here $k$ denotes address width, including the leading zero, and matrices act on columns. The inherited alphabet $\mathcal G_1=\{A,B,C\}$ preserves $3h-R=1$; its $h$-maps are $82h-21,5976h-1579,96264h-24407$.

There is a different physical binary alphabet $\mathcal G_2=\{P,Q\}$ preserving $5h-2R=1$. Its literal moves and products are
$
\begin{aligned}
P&=(8,101):&\quad M_P&=\begin{pmatrix}37&102\\90&256\end{pmatrix},
&h'&=292h-51,\\
Q&=(9,169)\ \text{then}\ (3,2):&
M_Q&=\begin{pmatrix}861&2046\\2150&5116\end{pmatrix},
&h'&=5976h-1023.
\end{aligned}
$
These follow by direct substitution, and $(5,-2)M=(5,-2)$ in each case. All addresses have a leading zero; there are no restrictions on word adjacency or repetitions. Starting from $h=1$, both maps preserve $h\equiv1\pmod4$. On that residue $P$ gives $h'\equiv1\pmod{16}$, while $Q$ gives $h'\equiv9\pmod{16}$. The inverse checks $h=1$ first, otherwise selects the last letter by these residues, undoes it by exact affine division and recurses on a smaller positive count. Reject invalid residues and divisions. Thus every $\mathcal G_2$ word has a unique physical history at its exact count. The prior $\mathcal G_1$ inverse uses $v_2(27h-7)$ for its terminal $A$-run and the mod-3 residue for the preceding letter.

For a third, unary chart, $\mathcal G_3=\{D\}$ with $D=(6,26)$ has matrix $\bigl(\begin{smallmatrix}19&27\\42&64\end{smallmatrix}\bigr)$, preserves $7h-3R=1$, and has $h'=82h-9$. Exact iteration is invertible by division after testing the seed. All three families preserve a unit affine line but its slope is not universal.

## 2. Three literal +2 charts

Take the previous four-channel graft: ordinal blocks $I\to B\to O$, reverse only $x\to i$ and $o\to y$, with $i,o$ the flank terminals and $x,y$ distinct middle marks. Converse the left flank. Its signed start sum and the right signed end sum vanish because the original rooted flanks have no full path starting at the terminal. For middle signature $(g,r,l,m)$, the ordinary count is
$
N=g hk+rR_I k+l hR_O+mR_I R_O.
$
Here $g=H(B)$, $r$ counts a selected path ending at $x$ and its complement, and $l$ does the same for a path starting at $y$. The joint coordinate splits $m=M_1+M_2$: $M_1$ counts a selected directed $y$-to-$x$ path and an ordered two-slot cover of its complement, and $M_2$ counts disjoint paths ending at $x$ and starting at $y$ plus a Hamilton path on the rest. The two terms distinguish the case where the two reversed arcs lie on the same odd cycle from the case where they lie on distinct odd cycles. In both cases the signed terms in the oddness filters vanish on the rooted flanks, leaving exactly these products. This argument also gives the single-arc terms. Every listed middle has $y$ reaching $x$, so the reverse arcs join the flank paths into a strong assembled tournament.

Give the left flank family $\mathcal G_q$ and the right $\mathcal G_t$, for $q,t\in\{1,2,3\}$. Their unit equations are $(2q+1)h-qR_I=1$ and $(2t+1)k-tR_O=1$. Thus in integer increment coordinates
$
h=1+qu,\quad R_I=2+(2q+1)u,\qquad
k=1+tv,\quad R_O=2+(2t+1)v.
$
This yields $N=Auv+Bu+Cv+d$. For a source switch to add a constant $s$ on the whole independent grammar, the change in $(g,r,l,m)$ must be
$
s\big((2q+1)(2t+1),-q(2t+1),-(2q+1)t,qt\big).
$
Indeed this tensor reads $s((2q+1)h-qR_I)((2t+1)k-tR_O)=s$. Conversely, an affine bilinear formula constant on a product of two unbounded sets has zero coefficients of $uv,u,v$, giving exactly this one-dimensional direction. The slope-$3$-only target had needlessly fixed the direction.

All masks below are outgoing and use the standardized middle marks $x=0,y=1$. The column $d$ simultaneously identifies the lower or upper source of the chart.

| $q,t$ | $d$ | Outgoing masks | $g,r,l,M_1,M_2$ |
|---|---:|---|---|
| 2,1 | 6207 | 220,117,24,18,192,29,172,46 | 377,468,757,476,369 |
| 2,1 | 6209 | 158,180,40,66,76,25,39,124 | 407,456,747,303,546 |
| 2,2 | 7507 | 148,93,16,5,104,15,173,62 | 251,498,760,619,566 |
| 2,2 | 7509 | 86,60,40,65,44,73,150,63 | 301,478,740,461,732 |
| 1,3 | 6239 | 52,13,8,161,78,86,143,55 | 317,562,567,522,394 |
| 1,3 | 6241 | 30,148,32,22,100,139,47,93 | 359,548,549,334,588 |

The first two charts are independently branching on both flanks. The third chart has a unary flank and provides a useful distinct arithmetic phase. The coefficients are:

| $q,t$ | $A$ | $B$ | $C$ | $d_{\rm low}$ |
|---|---:|---:|---:|---:|
| 2,1 | 20311 | 14572 | 8654 | 6207 |
| 2,2 | 43209 | 17882 | 18144 | 7507 |
| 1,3 | 29214 | 8633 | 21116 | 6239 |

The tabulated signatures are finite exact counts, from the end-path subset recurrence with singleton initializations, a $y$-starting variant and the complementary subset sums in the four-channel proof. The companion replay also checks literal assembled tournaments with small transitive flanks. The archive guided discovery; no completeness of that archive is needed for the infinite constructions.

## 3. Absolute, not just modular, inverse

For each of the six literal sources and each proposed ordinary target $N$, apply the exact identity
$
(Au+C)(Av+B)=A(N-d)+BC.
$
On physical states the two factors are positive. For every ordered divisor pair of the right side, recover $u,v$ by exact division, rejecting negatives and wrong residues modulo $A$. Then run both independent history inverses at $1+qu,1+tv$. Every accepted pair produces the target and every physical pair is recovered. This is a finite absolute membership test, not an efficient factoring claim or uniqueness assertion across distinct factorizations. The two branching charts have disjoint source variants modulo 8: their respective lower/upper classes are $7/1$ and $3/5$. Thus the target itself separates those four sources before any factoring; the additional unary chart can overlap and is tested separately.

Both independent branching charts give infinitely many distinct physical exact neighbor pairs; fixing one history and lengthening the other suffices. This is still a sparse absolute atlas. All three single-flank alphabets have word-count exponents strictly below one, since the sum of inverses of their expansion factors is less than one. The two-flank bound $O(X^\sigma\log X)$, with $\sigma$ the larger of the two one-flank exponents, follows from a dyadic product sum as for the +6 construction. No assertion is made that most or all target integers arise.

## 4. A dyadic atlas through 256 and a sharp next test

The arithmetic on the increment variables starts at 0 and is as follows:
$
\begin{array}{c|l}
\mathcal G_1 & 82 u+60,\ 5976 u+4396,\ 96264 u+71856\\
\mathcal G_2 & 292 u+120,\ 5976 u+2476\\
\mathcal G_3 & 82 u+24 .
\end{array}
$
For modulus $M$, compute each reachable set $S$ by closing $\{0\}$ under its displayed maps. To distinguish a mere hit from arbitrarily long histories, start from $S$ and iteratively replace the current set by the union of its forward images, until stabilization $S_\infty$. The set after $j$ iterations consists exactly of phases with histories of length at least $j$, by splitting off the final $j$ letters. A phase in the limiting set therefore has arbitrarily long histories, with unbounded counts because the physical maps grow strictly. Conversely outside it histories have bounded length. A target residue is infinite in the finite atlas exactly when at least one contributing input phase is in $S_\infty$, since all count formulas grow strictly in either variable.

The companion replay performs precisely this finite closure. At $M=256$ the output residue set with arbitrarily long histories is all 128 odd residues. At $M=512$, out of the 256 odd classes, the only ones never hit are $411,413$; the only additional classes without infinitely many distinct outputs are $479,481$. In fact the latter two classes contain **exactly one count each** in this atlas: $621848031,621848033$, from the first chart at empty left history and right history $C$. To certify this, enumerate only histories that remain outside $S_\infty$ at every prefix. The stable set is forward-invariant, so all histories finishing outside it are included. At modulus 512 this leaves 11, 6 and 6 actual input addresses respectively in $\mathcal G_1,\mathcal G_2,\mathcal G_3$. Their finite product tests give just that pair. These are obstructions for the literal atlas, not for Chan–Pak itself. In particular, querying congruences only to the successful precision 256 cannot distinguish this sparse atlas from full coverage.

## 5. A physical twin-prime witness and an exact provenance filter

In the first chart, modulo 3 one has $u\in\{0,1\}$ and $v\in\{0,1\}$. Here $u=0$ means the left history never uses $Q$; $v=1$ means the last right letter other than $A$ is $B$. The lower count reduces to $uv+u+2v\bmod3$. It follows that **both** members avoid 3 exactly when $u=0,v=1\bmod3$. In the symmetric $(2,2)$ chart the lower count reduces to $1+2u\bmod3$, so one of the two outputs is always divisible by 3. Merely having a +2 switch does not imply a viable prime-pair phase.

The actual histories $P$ and $AACAAB$ in the first chart give $h=241,\ k=19266225173740301$, and the twin primes
$
47\,124\,685\,853\,114\,258\,307\,047,\qquad
47\,124\,685\,853\,114\,258\,307\,049.
$
The replay certifies both by complete predecessor factorizations and modular-order (Lucas) tests; the smaller auxiliary primes are verified by trial division. This is one finite twin-prime result, not infinitely many primes in these charts.

The substantive transfer from [Pasten–Vidaux §§4–6](https://arxiv.org/pdf/2611.24748) is the separation of *local phase* from its exact integer witness: the phase selects eligible histories, while the determinant factorization and rooted inverse recover actual witnesses. It also shows concretely why a low-modulus approximation can miss the obstruction at the very next scale. The [linked RH/GRH remark](https://alpo.ge/remark.pdf) motivates tracking what repeated lifts preserve: changing the invariant affine line changes the rank-one tensor that an actual middle switch must realize. Neither paper supplies a theorem used here, and the announcement is not treated as an RH/GRH proof dependency.

**Next:** increase the set of invariant slopes and find middle signatures in their corresponding tensor directions; a useful next criterion is whether the new chart fills a missed $512$-phase *with arbitrarily long histories*, not just at its seed. An absolute covering mechanism remains open.

Exact finite replay: `python3 check_unit_switches.py`.
