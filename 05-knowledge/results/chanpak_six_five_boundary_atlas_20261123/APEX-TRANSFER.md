# An apex cut reduces the count to a positive word plus a nearest-neighbor signed word

**PROVED scoped / FINITE-EXACT.** The [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) is inherited; the two-ended observation statement concerns only ordinal apex contexts.

Write $H(U)$ for the number of Hamilton paths, with $H(\varnothing)=1$. Use the odd-cycle-packing expansion $H(U)=\sum_{\mathcal C}2^{|\mathcal C|}$ over sets of mutually vertex-disjoint directed odd cycles of length at least three (cyclic rotations are identified). Fix an apex $z$, and let $A_1\to\cdots\to A_r$ be the SCCs of its deletion. The argument also works with any nonempty ordinal blocks in place of maximal SCCs. Set $D=\prod_i H(A_i)$.

For a single block $A$, sum over its nonempty directed paths $P$, assigning weight $t^{|V(P)|}H(A-V(P))/H(A)$, for $t\in\{1,-1\}$. Let $a(t)$ be this sum with the restriction $z\to\operatorname{first}(P)$, $b(t)$ with $\operatorname{last}(P)\to z$, $c(t)$ with both restrictions, and $d(t)$ be one plus the unrestricted sum. Use the transfer
$$
M_A(t)=\begin{pmatrix}1&a(t)&c(t)\\0&d(t)&b(t)\\0&0&1\end{pmatrix}.
$$
Define $M(t)=\prod_i M_{A_i}(t)$ in ordinal order. Then
$$
H(T)=D\bigl(1+M(1)_{13}+M(-1)_{13}\bigr).
$$
Indeed, packings avoiding $z$ give $D$. A cycle through $z$ is one path obtained by traversing a strictly increasing sequence of blocks, with the endpoint restrictions and allowing intermediate blocks to be skipped; the other packing cycles are in the deletion remainder and contribute the stated products. The sum at $t=1,-1$ assigns coefficient two exactly to paths with even vertex order, the ones yielding eligible odd cycles. The ordered transfer reads exactly these choices.

## Exact locality of the signed channel

One has
$$
H(A)d(-1)=\sum_{S\subseteq V(A)}(-1)^{|S|}H(A[S])H(A-S)=0
$$
for **every nonempty** block, regardless of parity. For a fixed union odd-cycle packing of $A$, each cycle (odd vertex order) and each uncovered singleton can be independently assigned to either argument of $H$; changing its assignment reverses the sign. At least one such item exists. Thus
$$
M(-1)_{13}=\sum_i c_i(-1)+\sum_{i=1}^{r-1}a_i(-1)b_{i+1}(-1).
$$
In particular signed interactions between nonadjacent blocks vanish exactly, without positivity or asymptotic approximations. This extends to coloring each block among ordered slots whose signs are $\pm1$: any fixed packing contributes its usual weight times $(m_+-m_-)^{n-\sum_{C\in\mathcal C}(|C|-1)}$. The exponent is the number of cycles plus uncovered singletons; each such item has odd order. Thus signed slot covers are the ordinary slot polynomial evaluated at the difference of the two signed slot counts.

## Seven literal H probes recover the entire two-ended context observation

Consider a word $W$ of such blocks, its normalized transfers $M_W(t)$ with entries $a_\pm,b_\pm,c_\pm,d_\pm$, and its deletion mass $D$. This applies without assumptions of strong connectivity on the probes. Let O and I be new singleton blocks with apex contacts $z\to O$ and $I\to z$, respectively; other cross-block arcs follow word order. The following seven absolute Hamilton counts recover every count in an arbitrary left-right ordinal context of $W$:
$$
W,\quad WI,\quad WII,\quad OW,\quad OOW,\quad OWI,\quad WOI.
$$
Write $F(X)$ for the actual count of the apex tournament built on word $X$. The differences $F(WI)-F(W)$, $F(WII)-F(W)$ are $D(a_+-a_-)$, $D(3a_+-a_-)$, so recover both scaled $a$'s. The prepending probes analogously recover $Db_\pm$. The mixed difference for $OWI$, subtracting the two one-sided observations and adding the empty-context observation, is $D(d_++d_-)$. Since $d_-=0$ for nonempty $W$, this recovers $Dd_+$. Finally $F(WOI)-F(W)=2D+2Da_+$, recovering $D$, and $F(W)-D=D(c_++c_-)$. These are all the required data: in upper triangular multiplication, when $W$ is in the center, its two top-right entries only enter through their sum. Translation into probes needs no negative, fractional, or synthetic tournament weights.

This is a word in an upper triangular monoid, not a free group: the negative channel has a zero middle diagonal and cannot be inverted. The two channels distinguish distant positive interaction from genuinely local signed interaction.

## Exact 7 / 6 / 5 linear ranks

For each nonempty word, $a_-$ comes from its last block only, and $b_-$ from its first block only, by multiplying the signed transfer matrices with zero middle diagonals. Thus if one end block and its apex contacts are fixed, its negative-channel contribution is a known multiple of $D$ and six observables suffice. With both end blocks and contacts fixed, five suffice: $(D,Da_+,Db_+,Dd_+,H(T))$. For guards $O$ initially and $I$ finally, the signed exterior entries vanish outright. These guards make the actual apex tournament strong: the apex reaches the first singleton and then every later block, and the last singleton receives from every earlier block and reaches the apex.

These are exact **linear** ranks, not claims about nonlinear encodings. Regard the ordinal left/right contexts as columns of a matrix and candidate central words as rows. The seven probes above are columns, in the order shown. For all singleton words, the first seven rows below have determinant $128$:

| central word | W | WI | WII | OW | OOW | OWI | WOI |
|---|---:|---:|---:|---:|---:|---:|---:|
| O | 1 | 3 | 5 | 1 | 1 | 5 | 5 |
| I | 1 | 1 | 1 | 3 | 5 | 5 | 3 |
| OO | 1 | 5 | 11 | 1 | 1 | 9 | 9 |
| OI | 3 | 5 | 9 | 5 | 9 | 11 | 9 |
| IO | 1 | 3 | 5 | 3 | 5 | 9 | 5 |
| II | 1 | 1 | 1 | 5 | 11 | 9 | 3 |
| OOI | 5 | 11 | 23 | 9 | 17 | 23 | 19 |
| OIO | 3 | 9 | 19 | 5 | 9 | 19 | 15 |
| OII | 5 | 9 | 17 | 11 | 23 | 23 | 15 |
| OOII | 11 | 23 | 47 | 23 | 47 | 51 | 37 |
| OIOI | 9 | 19 | 39 | 19 | 39 | 45 | 31 |

With first guard $O$, drop column OOW, take rows (O, OO, OI, OOI, OIO, OII); the determinant is $64$. With first guard $O$ and final guard $I$, also drop WII, take rows (OI, OOI, OII, OOII, OIOI); the determinant is $64$. The counts can be checked by ordinary permutation enumeration or the singleton matrices (one of $a,b$ is $t$, the other is zero, $c=0$, and $d=1+t$).
