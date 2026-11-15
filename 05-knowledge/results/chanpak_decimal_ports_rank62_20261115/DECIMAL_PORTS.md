# The decimal ladder as a sharp multi-port signed readout

Write $L_k=1+1000^k$ for $k\ge1$, with a separate transitive starting value $1$. The relation $L_{k+1}-1=1000(L_k-1)$ is an exponential **dilation**, not the unipotent step $x\mapsto x+1$. Here is a fixed multi-port model for the dilation, with an exact obstruction to removing its normalization in unmarked ordinal layouts. It complements the [literal two-apex signed counter](TWO_PORT_SIGNED.md).

## Why fixed unmarked linearly growing blocks cannot realize the raw ladder

Fix a tournament quotient, a set of growing vertices $i$, and fixed nonempty modules at all other positions. Inflate growing vertex $i$ by a transitive block $TT_{\alpha_i k+\beta_i}$, where $\alpha_i$ is a positive integer and $\beta_i$ any fixed integer; take $k$ large enough that all blocks are nonempty.

**If the growing vertices contain a directed triangle,** select $m=\min_i(\alpha_i k+\beta_i)$ vertices in each of its three blocks. Orders that visit the blocks cyclically give at least $(m!)^3$ different Hamilton paths on that induced set. Inserting all other vertices injectively gives at least that many in the full graph. This is superexponential in $k$ and cannot equal $L_k$ eventually.

**If the growing vertices are transitive,** collapse each run of a growing module in a Hamilton path to its quotient label and keep every fixed actual vertex as its own symbol. Let $a_{\mathbf r}$ count the schedules with $r_i$ visits to growing block $i$. Each excursion between consecutive visits to $i$ must contain a distinct fixed actual vertex: a closed walk using growing labels only is impossible in their transitive quotient. Independently for the test of each $i$, permuting its interior excursions acts freely on these schedules and preserves every visit count. Consequently
$$
 (r_i-1)!\mid a_{\mathbf r}\quad\hbox{for each }i,\qquad
 H(k)=\sum_{\mathbf r}a_{\mathbf r}\prod_i r_i!S(\alpha_i k+\beta_i,r_i).
$$
The sum is finite with $r_i\le1+\#\{\text{fixed vertices}\}$. The filling factor counts surjections from the ordered vertices of a transitive block onto its ordered nonempty runs. Expand each factor as $\sum_{d=1}^{r_i}(-1)^{r_i-d}\binom{r_i}{d}d^{\alpha_i k+\beta_i}$, valid for positive block length. A term at exponential base $1000$ therefore has $\prod_i d_i^{\alpha_i}=1000$. All $d_i$ use only primes 2 and 5, and some $d_i\ge5$; its schedule coefficient is divisible by $4!$. With negative shifts the denominator can contain only 2 and 5. Thus the coefficient of $1000^k$ is in $3\mathbf Z[1/10]$, and is divisible by 3 if it is an integer. Distinct positive exponential bases are linearly independent on every infinite tail.

**Consequences.** No fixed unmarked substitution layout of this kind realizes $1+1000^k$ eventually, however many linearly growing transitive blocks it uses. Even a finite signed integer combination of such layouts with transitive growing quotients has coefficient divisible by 3 when its coefficient at base $1000$ is an integer. This rules out the entire fixed layout mechanism, not multi-contact, variable quotients, or realizability of individual rungs. It also explains why the $\{2,3\}$-smooth absolute rails do not simply accept 5 as a third repeated tile.

## A sharp signed dilation with actual, fixed deletion ports

Let $P_U(j)$ count ordered covers with $j$ slots, allowing empties. The [integral unmarked interface](../chanpak_symmetry_interfaces_20261115/OBSERVABLES.md) identifies the native coordinates $d_r=(r-1)!c_r$, where $c_r$ counts ordered covers with no empty slots. To make the signed clock concrete, use five strong marked hosts:
1. the directed triangle;
2. the natural order-four tournament with its extreme arc reversed, marking vertex 1 in zero-based order;
3. for $r=3,4,5$, the selector on $z$ and pairs $O_i,I_i$, $1\le i<r$: $z\to O_i$, $I_i\to z$, $O_i\to I_i$, $I_j\to O_i$ if $j\ne i$, with transitive orders within the two groups.

If $B_j(U)$ denotes the count after filling the marked vertex of host $j$ by $U$, the normalized coefficient rows in $(d_1,\dots,d_5)$ are
$$
\begin{pmatrix}
3&1&0&0&0\\
5&2&0&0&0\\
13&10&1&0&0\\
91&101&22&1&0\\
877&1250&411&40&1
\end{pmatrix}.
$$
The rows are exact counts of finite excursion templates; the [replay](verify_decimal_probes.py) constructs them independently from the arcs. Comparing row vectors gives
$$
P_U(2)=B_2(U)-B_1(U),\qquad
24P_U(5)=652B_1(U)-566B_2(U)+149B_3(U)-20B_4(U)+B_5(U).
$$
Using $W_2=TT_{3k-3}$ and $W_5=TT_{3k}$ for $k\ge2$ gives
$$
\bigl[B_2(W_2)-B_1(W_2)\bigr]
\bigl[652B_1(W_5)-566B_2(W_5)+149B_3(W_5)-20B_4(W_5)+B_5(W_5)\bigr]
=3\cdot1000^k.
$$
Each term after expanding the product is one actual Hamilton count: take the ordinal sum of its two filled probes. This realizes the **sharp** integer multiplier 3 in the signed acyclic class above. Putting both radices into a single unshifted block $TT_{3k}$ would cost $9!$: the highest-slot term at 10 has excursion factor $9!$, and the integral triangular probes extract $9!P(10)$. Separating the 2 and 5 packets and shifting the binary one avoids all of that factor except 3.

All ten counts can be read by fixed deletion sets of **one** tournament template. Share $W_2$ between two disjoint copies of the probe exteriors, share $W_5$ between five, and choose arbitrary fixed arcs between unused exterior copies. Orient the entire first scaffold toward the second. Retain one exterior per scaffold when reading a coordinate. The other exterior copies are deleted, so their mutual arcs never affect the readings. For any fixed digit parameter $t$, a third ordinal factor containing the one-spine, two-apex signed counter gives a fixed set of deletion-port readings with signed sum
$$
6t\cdot1000^k.
$$
Increasing $t$ gives a unipotent register after dividing by the current radix stride; increasing $k$ dilates the signed channel. These are physical **port measurements**, not a new tournament whose ordinary count equals their signed sum. In particular a counter whose stride varies with $k$ is not a proof of a stationary absolute-unit radix operation. The outstanding task is a nonlinear contact that converts such signed information into a positive ordinary count without the factorial cost or the strict-extension barrier.

## Two further stationary mechanisms that cannot carry a nonbinary radix

For a fixed unmarked host and a repeated **arbitrary** module $U^{\triangleright k}$ (ordered sum of $k$ copies), one has $P_{U^{\triangleright k}}(j)=P_U(j)^k$. Expanding the host's nonempty-slot formula, its largest base is $P_U(R)$, with coefficient $a_R$, where $R$ is its greatest allowed number of visits. The bases $P_U(j)$ are strictly increasing for nonempty $U$ (one can use the new slot alone for a full path). To realize $1+B^k$, necessarily $a_R=1$, hence $R\le2$ by the excursion factorial. $R=1$ gives only a single exponential. $R=2$ gives $P_U(2)^k+(a_1-2)H(U)^k$. Equality forces $H(U)=1$, so $U$ is transitive and $B=P_U(2)=2^{|U|}$. Conversely a triangle host filled by $TT_{dk}$ supplies radix $2^d$. This excludes base 1000 even if the repeated ordinal core is not specified in advance.

The same iff holds for any **fixed one-apex ordinal repetition** $P R^k S$, permitting arbitrary fixed painted blocks inside the period. On the state $(h,N,S,d)$ of the inherited cut transfer, the matrix of one nonempty period $R$ has characteristic polynomial $x(x-a)^2(x-q)$, $a=H(R),q=P_R(2)\ge2a$. This follows from the triangular update in the order $h,d,S,N$: the four diagonal entries are $a,0,q,a$. Hence for $k\ge1$ every fixed final scalar is $Cq^k+(D+Ek)a^k$, with the zero mode killed. A literal output $1+B^k$ requires $a=1,q=B$, again transitive $R$ and a power-of-two radix. The converse uses repeated transitive all-marked tiles, a terminal unmarked guard, and a single apex. The two-deletion signed counter survives this barrier because its coefficient of the growing mode cancels **between** distinct physical readings; it does not cancel within either positive scalar.
