# Fixed-seed padding cannot sustain a nonzero neutral frame switch

This closes the proposed double-zero seed in [FRAME_REPAIR §4](FRAME_REPAIR.md): **it cannot exist**. More generally a fixed exterior tournament and fixed contacts admit only finitely many unequal transitive padding pairs with a nonzero, exact two-slot-neutral switch. Identical nontransitive pads cannot evade the finiteness. There are nevertheless genuine one-time padded switches, including an unequal-pad prime pair 1800 apart.

Notation: $f_T(q)$ counts ordered $q$-slot path covers allowing empty slots, $H=f(1)$ and $B=f(2)$. Add a vertex $v$ to an exterior $U$; fix a subset $I$ of $U$ that $v$ beats. Switch just $v$, and let $\Delta$ denote after minus before. Let $\mathcal E_q$ and $\mathcal S_q$ be the weighted end and start totals
$$
\mathcal E_q=\sum_{\varnothing\ne S\subseteq U}
       \sum_{x\in I\cap S}\operatorname{end}_{U[S]}(x)\,f_{U\setminus S}(q),
$$
with $\mathcal S_q$ using starts. Write $D_q=\mathcal E_q-\mathcal S_q$. The earlier notation $E_q,O_q$ splits $D_q$ according as the internal support $S$ has even or odd size.

## 1. Adjacent reads give the lost parity coordinate

**Identity:** $E_q-O_q=D_{q-1}$ for every integer $q\ge1$. In fact the sign-weighted end total is $-\mathcal S_{q-1}$, and the sign-weighted start total is $-\mathcal E_{q-1}$. Here the sign attached to $S$ is $(-1)^{|S|}$.

Work in the commuting squarefree ring $\mathbb Z[x_i]/(x_i^2)$. Let $A$ be the ordinary 0/1 adjacency matrix, $X=\operatorname{diag}(x_i)$, $u=X\mathbf1$, and
$$
s=(I-XA)^{-1}u,\qquad e=(I-XA^{\mathsf T})^{-1}u,\qquad
F=1+\mathbf1^{\mathsf T}s.
$$
The components of $s,e$ count rooted starts and ends. Expansion counts paths, since repetitions of a vertex vanish. On each support $S$, the coefficient in $F^q$ counts $f_{U[S]}(q)$ by division of the vertices into $q$ ordered slots.

The tournament condition $A+A^{\mathsf T}=J-I$ and $X^2=0$ give
$$
I+XA^{\mathsf T}=(I-X)(I-XA+u\mathbf1^{\mathsf T}).
$$
The rank-one update, applied to $u$, gives $-e(-X)=s/F$. Reversing the tournament gives $-s(-X)=e/F$, since reversal preserves $F$. Thus $e(-X)F=-s$ and $s(-X)F=-e$. Multiply by $F^{q-1}$, sum the marked components, and take the coefficient of the whole support. This proves both identities including at $q=1$ with $F^0=1$.

In particular $E_2=O_2=0$ would imply $D_1=0$. The proposed seed cannot have a nonzero offset at any positive symmetric padding width.

## 2. Exact transfer for two possibly unequal padding blocks

Take nonempty pads $L,R$ on the two sides of $U$ in ordinal order $L\to U\to R$. Make $v$ lose to every vertex in $L$ and $R$, and switch $v$ at the end. For compactness put $L_j=f_L(j)$ and $R_j=f_R(j)$. The exact laws are
$$
\begin{aligned}
\Delta H &=R_1 L_2\mathcal E_1-L_1R_2\mathcal S_1,\\
\Delta B/2&=R_2 L_3\mathcal E_2-L_2R_3\mathcal S_2+\Delta H.
\end{aligned}
$$
To prove the full transfer, the two-sided padding does not allow any path using both left and right ends to contribute to the signed change: both endpoints have the same contact bit. For a path ending in $U$, choose its subset of $L$ freely and leave $R$ unused in the path; paths starting in $U$ use the opposite rule. Before parity projection, the two marked endpoint series are therefore multiplied by $R_q L_{q+1}$ and $L_q R_{q+1}$, respectively, after counting the path on chosen padding vertices in one extra slot. Signed support replaces $q+1$ by $q-1$ on that pad, while switching the marked end/read by §1. Since an odd cycle through $v$ uses a positive even support away from $v$, the difference is
$$
\Delta f(q)=q[ R_q(L_{q+1}\mathcal E_q-L_{q-1}\mathcal S_{q-1})
              -L_q(R_{q+1}\mathcal S_q-R_{q-1}\mathcal E_{q-1})].
$$
At $q=1$ use $L_0=R_0=0$; at $q=2$ the two lower-index terms give precisely $\Delta H$ inside the bracket.

**Symmetric transitive case.** For $L=R=\mathrm{TT}_a$, $a\ge1$,
$$
\Delta H=2^a D_1,\qquad
\Delta B=2^{a+1}(3^aD_2+D_1).
$$
A nonzero offset can be neutral at at most one positive width. If at width zero the original switch is neutral, then $E_2=0$ implies $D_2=-D_1$. At positive width it has exact defect $\Delta B=2(1-3^a)\Delta H$; no fixed seed can remain neutral and move.

**Unequal transitive case.** For pad lengths $\ell,r\ge1$, set $d=\ell-r$, $A_d=3^d\mathcal E_2-2^d\mathcal S_2$ and $C_d=2^d\mathcal E_1-\mathcal S_1$. Then
$$
\Delta H=2^r C_d,\qquad \Delta B=2^{r+1}(3^r A_d+C_d).
$$
If $I$ is nonempty, all four marked totals are strictly positive. For neutrality with nonzero offset the signs of $A_d$ and $C_d$ must oppose. Hence $d$ lies strictly between $\log_{3/2}(\mathcal S_2/\mathcal E_2)$ and $\log_2(\mathcal S_1/\mathcal E_1)$, regardless of their order. This interval is finite. Each possible integer $d$ admits at most one $r$ since $3^r=-C_d/A_d$ must hold exactly. The empty contact set is inert. This gives both an all-order obstruction to unbounded fixed-seed asymmetric lengths and a finite exact selector; comparisons can use rational powers to avoid rounding.

**Identical arbitrary pads.** With $L=R=P$, a nonzero offset requires $D_1\ne0$; neutrality is equivalent to $f_P(3)/H(P)=-D_1/D_2$ (and $D_2\ne0$). For order $m$, $f_P(3)/H(P)\ge\binom{m+2}2$: cut every Hamiltonian ordering according to each weak composition of $m$ into three parts. Concatenating its slots recovers the original ordering and the three lengths, proving injectivity. Thus changing the shapes of both identical pads cannot keep a fixed seed neutral at arbitrarily large orders.

The same order bound applies if fixed, possibly different pads $L,R$ each acquire a common extra ordinal factor $P$. In normalized coordinates $t_L=L_2/L_1,\ u_L=L_3/L_1$ and likewise for $R$, the neutrality expression is $A+C=0$, where $A=t_R u_L\mathcal E_2-t_L u_R\mathcal S_2$ and $C=t_L\mathcal E_1-t_R\mathcal S_1$. Multiplying each pad profile by the profile of $P$ multiplies every $t$ by $f_P(2)/H(P)$ and every $u$ by $f_P(3)/H(P)$. At nonzero shift the latter ratio must therefore equal the single fixed number $-C/A$. No assertion is made for independently changing nontransitive pads.

## 3. Physical one-time moves and their prime read

The background pass supplied the identity and three symmetric-width examples. The unequal pads and the two directed-triangle pads below were then constructed and replayed from literal arcs. $I$ is recorded as a bit mask over the exterior in [padding_replay.py](padding_replay.py).

| Pad pair | Exterior order | $H$ before and after | Shared $B$ |
|---|---:|---:|---:|
| Transitive $1,1$ | 6 | 57, 75 | 3696 |
| Transitive $2,2$ | 7 | 709, 529 | 107584 |
| Transitive $3,3$ | 8 | 15029, 15461 | 6443776 |
| Transitive $2,3$ | 6 | 905, 725 | 138368 |
| Transitive $3,2$ | 8 | 7013, 8813 | 2032160 |
| Directed triangle, directed triangle | 8 | 69975, 69579 | 11628288 |

The triangle pair is selected by the exact ratio 11, not by a power of 3: $D_1=-11$, $D_2=1$ and $f_{C_3}(3)/H(C_3)=33/3$. It realizes a same-lattice step of 396. The second unequal transitive pair has $(\mathcal E_1,\mathcal S_1,\mathcal E_2,\mathcal S_2)=(1073,1696,8012,12043)$. It gives a same-lattice step of 1800. With the established slope-6 positive reader and 16 ordinal singletons this gives the prime quotients $799077833573$ and $799077835373$ and their exact sevenfold partners. The replay checks both primes by complete trial division and checks the switches by subset DP independently of the transfer formulas.

**What survives:** the earlier ordinal growing-width compiler uses a selected inner frame, and is unaffected. **What must change:** to get indefinitely many nonzero neutral moves within this enlarged lattice architecture, one must change the marked exterior with scale or independently control the padding profiles (or use a different interface). Simply widening the pads, including with common nontransitive additions, cannot solve the prescribed-baseline problem.
