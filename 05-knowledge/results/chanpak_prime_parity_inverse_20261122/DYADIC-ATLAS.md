# Repairing the 411/413 phase by a new physical move

**PROVED.** One new legal migration letter repairs both missing $411,413\bmod512$ phases of the previous +2 atlas. It also turns its $479,481\bmod512$ phases into infinite families, without loss of the exact history inverse. Adding two literal middle charts gives infinitely many counts in **every odd class modulo $2^{18}$**; at modulus $2^{19}$ this enlarged atlas has exactly two omitted odd phases. In every odd class modulo $2^{13}$, both flanks can grow independently and branch freely. None of these modular results implies full absolute coverage. There is also a literal one-arc switch with the **absolute values 411 and 413**.

## 1. An extra letter that preserves the unit line but changes the phase

Use the existing rooted state $h=H$, marked-end complementary-cover count $R$, exactly one full path ending at the root, and none starting there. A leading-zero migration of width $k$ and address $j$ uses
$
M(k,j)=
\begin{pmatrix}F(j)&j+1\\2^k-\mathrm{rev}_k(j)&2^k\end{pmatrix},
$
where $F$ is the one-apex word count. This is the physically rooted transfer of the [migration construction](../chanpak_bank_fusion_20261117/migration/PROOF.md). The original letters $A,B,C$ preserve $3h-R=1$, and act by $82h-21,5976h-1579,96264h-24407$.

Append a new letter $E$, with chronological address sequence $(6,8),(9,72),(10,340)$. The middle move is not a return on its own. Together the three matrices give
$
M_E=\begin{pmatrix}13371651&14153823\\40114950&42461470\end{pmatrix},
\qquad (3,-1)M_E=(3,-1).
$
All three addresses have leading zero: every $ABCE$-word, including arbitrary repetitions and adjacency, is physical. On the unit line the new map is
$
h\longmapsto55833120h-14153823.
$
It is strictly increasing relative to its input on $h\ge1$, and preserves $h=1\bmod4$.

**No inverse loss.** In the shifted coordinate $z=27h-7$, $A$ multiplies by 82 and increases $v_2(z)$ by one. $B,C$ reset $v_2(z)$ to 3 and 2 respectively, and set $h\bmod3$ to 2 and 1; the seed has valuation 2 and residue 1. $E$ acts by $z\mapsto55833120z+8678612$. The multiplier has valuation 5 and the added term has valuation 2, so $E$ always resets to valuation 2; its $h$-residue is 0 modulo 3. Thus the trailing $A$-run has length $v_2(z)-3$ if $h\equiv2\bmod3$, and $v_2(z)-2$ otherwise. Undo that run (validating exact division), stop if the seed is reached, else decode the preceding letter from $h\bmod3=2,1,0$ as $B,C,E$. Reject wrong residues, nonpositive predecessors and failed exact divisions. Strict growth means termination, and forward multiplication certifies every accepted word.

This is valuation propagation without identifying a coarse residue with an exact witness. The same signed unit line enables old middle switches unchanged, but its residue orbit genuinely changes.

## 2. Explicit, absolute witnesses for the formerly missing phases

The old atlas's third chart couples the $ABCE$ flank to the unary $D$ flank; the two absolute outputs are
$
N_s=29214uv+8633u+21116v+6239+2s,\qquad s\in\{0,1\},
$
where $u=h-1$ on the $3h-R=1$ line, and $v=(k-1)/3$ on the $7k-3R_O=1$ line. In these coordinates the positive affine maps are
$
A:u\mapsto82u+60,\quad C:u\mapsto96264u+71856,\quad
E:u\mapsto55833120u+41679296,\quad D:v\mapsto82v+24.
$
Take $D$ on the unary flank ($v=24$). On the other flank, take any $ABCE$-word $W$, then append $CEA$ in chronological order. The multiplier of this suffix has valuation $3+5+1=9$, and its value at the seed is $328982880701372\equiv444\bmod512$. Therefore the actual numbers are $709769u+513023$ and $709769u+513025$, always respectively $411,413\bmod512$. For empty $W$ they are
$
233501850252532616091,\qquad233501850252532616093.
$
Taking arbitrary longer prefixes gives infinitely many distinct values in both phases, by strict growth. For the old formerly finite-only classes, use the old first chart with empty left history and right history $W AEC$. This suffix again has multiplier valuation 9, and gives $v=176\bmod512$, hence the two counts are $479,481\bmod512$. Both guarded families are completely physical, not progression heuristics.

For the exact small targets, an 8-vertex source with outgoing masks
$
(10,192,35,134,15,27,61,117)
$
has $H=411$; reverse only the arc on vertices $5,6$ (zero-based) to obtain
$
(10,192,35,134,15,91,29,117)
$
with $H=413$. The exact same subset recurrence verifies analogous single-arc switches $479/481$ and $551/553$ in the replay. These are explicit small absolute realizations, separate from the growing atlas.

## 3. Two extra literal middle pairs

For each middle $B$, use the same ordinal graft $I\to B\to O$, with only reverse arcs $x\to i,o\to y$, the left rooted flank conversed, and $x=0,y=1$. The four-channel count is
$
N=g hk+rR_I k+l hR_O+mR_I R_O
$
for $g=H(B)$, marked end/start reads $r,l$, and $m=M_1+M_2$. As in the original graft, $M_1$ counts a selected $y$-to-$x$ path plus an ordered two-path cover on its complement, and $M_2$ counts disjoint end/start paths at $x,y$ and a path on the rest. The two cases correspond to the reverse arcs on one or two mixed odd cycles. The oddness filters lose no information because the signed marked sums of the flanks are zero. Each listed middle has $y$ reaching $x$, hence the complete assembly is strong via the two root paths and return arcs.

The left line is now $5h-2R_I=1$ (the independently branching $P,Q$ alphabet); the right is $3k-R_O=1$ (the enlarged $A,B,C,E$ alphabet). Use $u=(h-1)/2$ and $v=k-1$. The signed tensor requirement for a constant $+2$ is exactly $\Delta(g,r,l,m)=2(15,-6,-5,2)$, since it reads $2(5h-2R_I)(3k-R_O)$. Two literal middle pairs giving that tensor are:

| $d$ | Outgoing masks | $g,r,l,M_1,M_2$ |
|---:|---|---|
| 5855 | 208,25,3,21,164,15,62,110 | 425,623,570,384,377 |
| 5857 | 90,48,35,134,44,137,62,87 | 455,611,560,265,500 |
| 7787 | 60,77,72,192,14,30,177,55 | 459,650,810,661,441 |
| 7789 | 142,28,80,52,97,71,139,62 | 489,638,800,418,688 |

Their count formulas are $18800uv+13855u+7947v+5855/+5857$ and $25558uv+18428u+10801v+7787/+7789$, respectively. The table uses exact marked subset counts; its replay checks literal assembled grafts with several small transitive flank pairs. For any chart $N=A_0uv+B_0u+C_0v+d$, retain the exact integer membership inverse
$
(A_0u+C_0)(A_0v+B_0)=A_0(N-d)+B_0C_0.
$
Factor the positive right-hand side, recover nonnegative $u,v$ by exact division, and decode each history in the selected alphabet. When sources overlap modulo a chosen modulus, test each of the finitely many source variants. No modular approximation is accepted as an absolute certificate.

There are also explicit **two freely growing prefixes** for the original missing classes, even modulo 1024: in the new $5855/5857$ chart use arbitrary left prefixes followed by $QPQP$, and arbitrary right prefixes followed by $BACB$. Their multipliers have valuations $3+2+3+2=10$ and $3+1+3+3=10$; their seed phases are $u=936,v=172\bmod1024$, giving $411,413\bmod1024$. The counts with both prefixes empty are $4920127039032761792995273678235$ and that value plus two. The second new chart similarly sends suffixes $PPQQ,EBB$ to $551,553\bmod1024$.

## 4. Exact finite-state parity theorem, with arbitrarily long histories

Use the three old charts plus the two new ones, feeding $ABCE$ wherever the $3h-R=1$ line appears. The charts have positive coefficients, and each source's upper count is always exactly two above its lower count. In the increment coordinates the affine maps, all seeded at zero, are:
$
\begin{array}{c|l}
3h-R=1 & 82u+60,\ 5976u+4396,\ 96264u+71856,\ 55833120u+41679296\\
5h-2R=1 & 292u+120,\ 5976u+2476\\
7h-3R=1 & 82u+24.
\end{array}
$
For each modulus $M=2^r$, close $\{0\}$ under these maps to obtain each finite orbit $S$. A second exact set $S_\infty$ stabilizes by repeatedly taking the union of images of the current set, starting from $S$. After $j$ iterations this is exactly the set of residues with a history of length at least $j$, by removing the final $j$ letters. Thus $S_\infty$ consists exactly of phases with histories of unbounded length. At a dyadic modulus it has a stronger suffix interpretation: a history of length at least $r$ has total multiplier divisible by $M$, so *every* preceding prefix can be admitted without changing its phase.

For each chart evaluate both outputs on its product of $S$'s. A residue has infinitely many counts exactly when one input belongs to $S_\infty$; if neither belongs, the finite collection of contributing residues has only bounded-length words. Evaluating on $S_\infty\times S_\infty$ characterizes the stronger two-suffix capacity. This proof is exact for all history lengths; the finite part is the explicit orbit iteration in the standalone replay.

| Atlas and criterion | Exact finite output |
|---|---|
| Three old charts plus $E$, at least one freely growing flank, modulo 512 | All 256 odd residues |
| Same atlas at modulo 1024 | Exactly $411,413,551,553$ omitted |
| Five-chart atlas, at least one freely growing flank, modulo $2^{18}$ | All 131,072 odd residues |
| Five-chart atlas, ordinary or freely growing, modulo $2^{19}$ | Exactly $141915,141917$ omitted |
| Using only the four genuinely branching charts, both flanks freely growing, modulo $2^{13}$ | All 4096 odd residues |
| Same stronger criterion modulo $2^{14}$ | Exactly $13915,13917$ omitted from the two-suffix phase set |

The finite calculations use that multiplying the congruence modulus is not the same as increasing the size of a tournament. They establish sharp frontiers for these specific finite banks and leave open whether further lines or middle choices can give all dyadic levels. In particular arbitrarily long histories in every class through $2^{18}$ do not solve the absolute inverse for any specified integer.

Exact finite checks: `python3 check_dyadic_atlas.py`.
