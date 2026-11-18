# A fixed prime support cannot sustain a fixed-core neutral switch

**Proved:** a marked core with a fixed finite set of primes in the normalized pad profiles permits only finitely many nonzero neutral switches, even when the two pads vary independently and their topologies are arbitrary. A separate exact inverse selects a mirror-marked transitive core from two pad profiles. **Physical example:** a pair using new pad primes 347 and 349 chooses that core and gives two verified prime counts 673,920 apart, and both sevenfold partners. This is not arbitrary-prime coverage.

## 1. Four positive terms and an S-unit obstruction

Keep the marked exterior $U,I$ of the [padding identity](PADDING_OBSTRUCTION.md) fixed, with $I$ nonempty, and use its positive four totals $\mathcal E_1,\mathcal S_1,\mathcal E_2,\mathcal S_2$. Put $L_j=f_L(j),R_j=f_R(j)$ for the nonempty pads. The exact two-slot neutrality law is
$$
 X+Z=Y+W,\quad \Delta H=Z-W,\qquad
 (X,Y,Z,W)=(\mathcal E_2 L_3R_2,\,\mathcal S_2 L_2R_3,\,
                  \mathcal E_1 L_2R_1,\,\mathcal S_1 L_1R_2).
$$
Every term is positive. A vanishing proper subsum of $X-Y+Z-W=0$ must split the terms into two cancelling pairs. If $X=Y$ and $Z=W$, the switch is inert. The only other possibility is $X=W$ and $Z=Y$. Those two equalities force
$$
 (L_3/L_1)(R_3/R_1)=
       \mathcal S_1\mathcal E_1/(\mathcal E_2\mathcal S_2)\le1,
$$
which is impossible: both left ratios strictly exceed 1 and $\mathcal E_2\ge\mathcal E_1$, $\mathcal S_2\ge\mathcal S_1$ by slot inclusion on the complement. Thus every nonzero neutral switch yields a **nondegenerate** equation.

Set $t_L=L_2/L_1$, $u_L=L_3/L_1$, and likewise for $R$. Dividing by $W$ gives the equation with fixed coefficients
$$
 \frac{\mathcal E_2}{\mathcal S_1}u_L
 -\frac{\mathcal S_2}{\mathcal S_1}\frac{t_L}{t_R}u_R
 +\frac{\mathcal E_1}{\mathcal S_1}\frac{t_L}{t_R}=1.
$$
If $t_L,t_R,u_L,u_R$ are positive rational numbers supported on a fixed finite prime set $S$, the variable triple $(u_L,(t_L/t_R)u_R,t_L/t_R)$ lies in the product of three copies of the positive $S$-unit group. Its rank is at most $3|S|$. By [Evertse, Schlickewei and Schmidt, Theorem 1.1](https://annals.math.princeton.edu/2002/155-3/p04), this gives finitely many triples. Each determines both $u_L$ and $u_R$. Each pad's order is then bounded because $f_T(3)/H(T)\ge\binom{|T|+2}{2}$ by cutting Hamiltonian paths into three possibly empty slots. Hence there are only finitely many pad pairs up to isomorphism. No effective order bound is claimed here.

Ordinal joins multiply profiles. Thus independently building the pads as arbitrary words from a fixed finite library cannot repair this, since the normalized prime support remains finite. In particular, restricting normalized readings to Ogg/Monster prime support will miss any genuinely unbounded fixed-core family. This is a property of the *pad readings*, not of the output prime or of an individual tournament's group order. It does not say that the number of primes in each individual example must grow.

## 2. A complementary exact height bound and an inverse in the other direction

The neutrality law can be written as one rational score collision:
$$
 \frac{\mathcal E_2L_3-\mathcal S_1L_1}{L_2}
 =\frac{\mathcal S_2R_3-\mathcal E_1R_1}{R_2}.
$$
For a fixed common score, both orders are effectively bounded, even without prime support restrictions. In fact for a tournament of order $n$ and integer $q\ge1$,
$$
 \frac{f_T(q+1)}{f_T(q)}\ge\frac{n+q}{q}.
$$
To prove it, choose one of the $q$ slots and split its path at any position, including before its first or after its last vertex. There are $n+q$ marked splits per input. A resulting $(q+1)$-cover has at most $q$ preimages, one for each possible merge of adjacent slots, whether or not the join edge exists. Consequently the left score of an order-$n$ pad is at least $\mathcal E_2(n+2)/2-\mathcal S_1/2$; the corresponding right bound is $\mathcal S_2(n+2)/2-\mathcal E_1/2$. These give explicit order bounds from a proposed common score. They do not by themselves limit the number of common *scores*.

There is also a finite inverse for the **core** when the two pad profiles are given. Restrict to $U=\mathrm{TT}_m$ in its linear order and to a nonempty contact set $I$ invariant under reversal. Then starts and ends agree under reversal, and
$$
 \mathcal E_1=\sum_{i\in I}2^i,\qquad
 \mathcal E_2=\sum_{i\in I}3^i2^{m-1-i}
$$
for zero-based positions. A marked end at $i$ can bring any subset from among the $i$ predecessors, leaving all successors in the other slots. This proves both formulas. A nonzero neutral switch requires a positive ratio
$$
 K=\frac{L_1/L_2-R_1/R_2}{L_3/L_2-R_3/R_2}
  =\frac{\mathcal E_2}{\mathcal E_1}.
$$
For nonzero $\Delta H$ the numerator is nonzero, so if the denominator vanishes there is no such switch. As a weighted average of $(3/2)^i 2^{m-1-i}$, $K$ lies between $(3/2)^{m-1}$ and $2^{m-1}$. Hence only
$$
 1+\log_2 K\le m\le1+\log_{3/2} K
$$
can work. For each such $m$, solve the exact subset-sum $\sum_{i\in I}(3^i2^{m-1-i}-K2^i)=0$ with the variables grouped into reversal pairs and possible midpoint. This is a finite exact inverse for this class of cores, not a way to ensure that a core exists at a prescribed absolute path count.

## 3. Outside the prime bank: a matched pair and two exact target primes

The background pass supplied a pair of order-eight pad cores:
$$
 L_0=(226,236,73,81,103,76,128,60),\quad
 R_0=(174,68,136,210,167,78,149,34).
$$
Append two ordinal singletons **within each pad**. The resulting profiles $(H,f(2),f(3))$ are $(359,22208,393903)$ and $(355,22336,396171)$. The normalized $f(2)$ numerators include 347 and 349, outside the Monster prime support. Their ratio $K=13/5$ selects only $m=3$ in the above permitted order interval. Of its nonempty mirror masks, only the endpoints give this ratio. They supply $(\mathcal E_1,\mathcal S_1,\mathcal E_2,\mathcal S_2)=(5,5,13,13)$ and the common score $461/2$.

The actual 24-vertex frames have $H=39965675,39291755$ at the shared $B=232641769472$. A single cut at $v$ links them. The replay recomputes pad profiles on subsets and counts outer paths independently by quotient words with up to three visits per pad, multiplying by the ordered nonempty-cover counts for the used visits. This checks the two-slot equality and both Hamiltonian counts without assuming the transfer law.

Use the paired slope-62 readers from the signed-height discovery table, with intercepts 1 and 7 and slopes 62 and 434. The two literal marked readers are checked in the replay. Three ordinal singletons on the input scale its common $B$ by eight. This gives two prime counts and both sevenfold partners:

| Prime count | Sevenfold partner |
|---:|---:|
| 115390357623787 | 807732503366509 |
| 115390356949867 | 807732498649069 |

The separation is 673,920. The replay proves primality by a deterministic Pocklington certificate: known prime factors of each $p-1$ have product greater than $\sqrt p$, and the base-2 congruences force every prime divisor of $p$ to be 1 modulo that product. A composite would then have no divisor at most $\sqrt p$. The isometry assertion is for the incoming frames, not for the graphs after joining on the reader.

**Remaining obligation:** find rational-score collisions across unbounded scores using new prime support, and couple the same choice to an absolute target. Fixed-smooth pad generation is now excluded, not just sparse in the samples. The exact core decoder and the sample outside the bank do not settle unbounded collisions, arbitrary prime placement, or Chan–Pak.

Replay: **python3 prime_support_replay.py** in this directory. The function **build(flip,scale)** materializes the two prime pairs.
