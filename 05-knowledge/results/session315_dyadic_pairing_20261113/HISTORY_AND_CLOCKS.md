# Four adjacent orders, two prime clocks, and three different meanings of prime support

## Square-sum transfer at 23–26 retains paths with ports, not just totals

Use $Q_n$ on $1,\ldots,n$, joining distinct labels with square sum. The prior [square-sum census](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md), [endpoint-aware order-24 obstruction](../decoder_prime_square_20260925.md), and [Monster/Ogg and square-sum note](../monster_levels_square_sum_7_21_20261022.md) already establish the existence pattern and analyze the ten paths at 25. The following comparison records the exact interface that changes between consecutive orders.

Throughout $Q_{23},\ldots,Q_{26}$, vertex 18 has sole neighbor 7 (its self-sum 36 is excluded, and $49-18=31$). Orient each Hamilton path from 18. Write $N$ for the new vertex, and $U_N$ for its earlier neighbors. A path ending at $N$ is precisely an old path from 18 to some $a\in U_N$ with $N$ appended. If $N$ is interior, deleting it gives an **ordered** pair of disjoint spanning paths of the old graph with endpoints $18\to a$ and $b\to w$, where $a,b\in U_N$ are distinct and $w$ is the final endpoint. Conversely glue those two paths through $a,N,b$. This is an exact recurrence for counts with chosen ports; it never assumes that $a,b$ themselves were adjacent. In particular the count of paths in the previous graph is not a sufficient state.

| $N$ | $U_N$ | append at end | ordered interior port counts | total |
|---:|---|---:|---|---:|
| 23 | $2,13$ | 0 | $2,13:2;\ 13,2:1$ | 3 |
| 24 | $1,12$ | 0 | both zero | 0 |
| 25 | $11,24$ | 0 | $11,24:8;\ 24,11:2$ | 10 |
| 26 | $10,23$ | 2 | $10,23:6;\ 23,10:4$ | 12 |

The entries can be checked by a marked bitset DFS beginning at 18; the recurrence is independent of the finite check. The zero at 24 does not kill the next row: ten order-24 two-cover interfaces survive. At 26 the two mechanisms coexist, unlike 23 or 25. This is a concrete example of why the tournament metagraph should retain endpoint-marked states or witnesses, not only size, Hamilton count, or a contraction with unlabeled boundary. It is not an isomorphism from the square-sum graph to the metagraph.

## The first twenty-five primes give an exact 72-state sieve with independent side registers

Consider $99\cdot2^k+23$, $k\ge0$. Its first-25-prime sieve means avoidance of every prime up to 97, not primality. The primes dividing it in this set, and the unique prohibited exponent classes for each, are as follows:

| prime | period | prohibited residue |
|---:|---:|---:|
|5|4|3|
|13|12|1|
|17|8|1|
|19|18|9|
|29|28|27|
|37|36|23|
|53|52|25|
|59|58|35|
|61|60|0|
|67|66|56|
|83|82|15|

For $p=2,3,11,23$ divisibility is impossible directly. For the other primes, compute the cyclic list of powers of 2 in $\mathbb F_p^\times$ and test whether it contains $-23/99$. This proves the table and that all remaining primes are inert for this particular expression. Rules for 29 and 37 are redundant: both imply $k\equiv3\bmod4$, already killed by 5.

The remaining periods factor into a base modulus 72 and five distinct side primes $13,29,5,11,41$. Reject $r=k\bmod72$ if it meets any of the rules for 5,13,17,19. At each remaining $r$, an outer rule is enabled iff $r$ agrees with its prohibited residue modulo, respectively, $4,2,12,6,2$, for 53,59,61,67,83. An enabled rule removes exactly one value of its associated side-prime register. Thus the number of survivors over a CRT period is obtained by multiplying five independent factors $p$ or $p-1$ at each of the 72 base states and summing. This gives period $L=61,209,720$, survivors $32,353,278$, and density $5,392,213/10,201,620$; the conditional densities by parity are $103/110$ (even) and $5600/46371$ (odd).

This period is minimal. For each squarefree side prime $p$, a live base state enabling its rule is respectively $r=5,5,0,2,5$. The shift $L/p$ changes just that side register. Choose its prohibited value and all others allowed; the shift changes a rejected exponent into an accepted one. For the 2-factor, $L/2$ sends base $17$ to $53\bmod72$, changing rejection by 17 to an allowed base state. For the 3-factor, $L/3$ sends base $45$ to $21$ or $69\bmod72$, removing rejection by 19. In both cases the side registers can be chosen to avoid all rules at the live base state and remain unchanged under the shift. Any least period divides $L$, and if proper it divides some $L/p$, so the tests suffice. The raw lcm before removing redundant rules is $428,468,040$.

The census and minimum period are a small arithmetic example of eliminating inactive contexts and retaining the exact shared interface. The 25 sieve primes are **not** the Ogg prime set; even 23, prominent in moonshine and the constant term here, never divides this expression.

## Real tournament length clocks at 23 and 25

Use the [apex overlap law](../session313_zeckendorf_target_selection_20261113/ADDITIVE_PORTS.md) and [length bypass](../session313_zeckendorf_target_selection_20261113/LENGTH_BYPASS.md). For an ordinal chain of $n$ identical blocks with a common apex and uniform cut, let $h=H(T)$, $q=Q_T(2)$, and $s^*,e^*$ be its two-slot start/end cut responses. If $p\mid h,\ p\nmid q$, all non-overlap contributions vanish modulo $p^k$ for $n\ge k+2$. The overlap formula gives
$$
H_n\equiv s^*e^*q^{n-2}(1-h/q)^{-2}\pmod{p^k}.
$$
Indeed the terms are $s^*e^*\sum_{t=0}^{n-2}(t+1)q^{n-2-t}h^t$; higher terms in the geometric derivative vanish modulo $p^k$.

* For **23**, take the strong tournament with outgoing masks $(48,1,3,7,14,30)$. It has $h=23,q=296$; apex out-cut mask 1 gives $(s,s^*,e,e^*,b)=(1,25,22,244,24)$. Here $q\equiv20\bmod23$ has order 22 and $q^{22}\equiv139=1+6\cdot23\bmod529$. Binomial lifting gives full unit period at every precision, and $25\cdot244$ is a unit. One fixed cut suffices; appending ordinal copies gives nonunits. In contrast $2$ has order 11 modulo 23, so merely multiplying a square response by neutral sink factors of 2 cannot move it into the nonsquare coset.
* For **25**, take rows $(50,8,3,5,14,30)$, with $h=25,q=304$. Add one sink within each block that beats the apex. The [sink rule](../session313_zeckendorf_target_selection_20261113/TERNARY_COUNT_CHART.md) leaves $h=25$ and gives $q'=608$, primitive at every $5^k$ since $608\bmod5=3$ and $608^4\equiv21\bmod25$. Cuts 1 and 2 have padded $(s^*,e^*)=(84,556),(56,503)$, whose products have opposite square classes mod 5. Each alone covers all unit residues by varying length. Ordinally joining the 5-path order-four core supplies arbitrary nonunit valuations.

All uniform constructions just specified for units can be taken strong: the apex reaches the first interior block, the last block can reach the apex, the unpadded blocks are strong, and every added sink is reachable from its old block and beats the apex. The periods $22\cdot23^{a-1}$ and $4\cdot5^{b-1}$ share only a parity condition. For simultaneous prescribed units at the **same** length $n$, use the 23-cut as specified and select the 25-cut by the square class that makes the two required parities agree. CRT gives the common length, arbitrarily large. This is a relation between physical response systems, not a connection to the square-sum paths at those sizes.

The 23-clock is still onto all its unit residue targets if its length $n$ is required to pass the entire first-25 sieve for $99\cdot2^n+23$. The common part of the length clock period and $L$ is 22. Every residue mod 22 has a live lift for the sieve: take base 4 mod 72 for an even class and 5 mod 72 for an odd class. Rule 67 is disabled at both bases, so the prescribed mod-11 register never obstructs. All other side registers can be chosen to pass. Then CRT with the remaining power of 23 gives arbitrarily large admissible lengths. This proves coexistence of two exact modular observers, not primality of the exponential expression or absolute path-count coverage.

## A literal valuation orbit behind the Monster labels

The prior [valuation comparison](../monster_levels_square_sum_7_21_20261022.md#3-the-exponent-cycle-numerology) gives $N=|M|/|B|$ with exponent vector $(5,7,3,4,1,2)$ at $2,3,5,7,11,13$, and four additional primes $29,41,59,71$ with exponent one. It correctly rejects a well-defined mod-9 interpretation of the prime-exponent map. There is a different, natural exact extension. For a base $b\ge2$, define
$$
f_N(b)=\max\{j:b^j\mid N\}.
$$
This is a divisor valuation, not an additive valuation for a composite base. Here $f_N(4)=\lfloor5/2\rfloor=2=f_N(13)$, so
$$
4\longmapsto2\longmapsto5\longmapsto3\longmapsto7\longmapsto4
$$
is a literal five-cycle, and 13 is a separate feeder to 2. No identification of 4 with 13 is needed; there is not a six-cycle. After one step all values lie in 0 through 7; declaring 0 and 1 terminal, the displayed cycle is the only nontrivial one because 6 maps to 5. If one instead uses the actual involution class index $N/2=[M:2.B]$, the 2-exponent drops to 4, giving the literal two-cycle $2\leftrightarrow4$, with $5\to3\to7\to4$. The distinction between order quotient and class index is visible in the dynamics.

There is an exact fiber enumerator. With $N=\prod p^{e_p}$, the number of bases of output $j\ge1$ is
$$
\prod_p(\lfloor e_p/j\rfloor+1)-\prod_p(\lfloor e_p/(j+1)\rfloor+1),
$$
by subtracting counts of $b^j\mid N$ and $b^{j+1}\mid N$. The vectors at $j=1,\ldots,7$ are $(92016,120,16,4,2,0,1)$ for $N$ and $(76656,120,16,6,0,0,1)$ for $N/2$. The order-quotient primes, the supersingular rationality filters at levels $1,2,6$ in the [Ogg/Baby Monster/Fi22 analysis](../monster_levels_square_sum_7_21_20261022.md#0-the-answers), and the active sieve primes above all arise from different maps. The genuine cross-topic principle is to compare the named filter and its fibers, not just the resulting list.

## A tournament prime-support filtration with an exact delayed mirror threshold

Any odd-order subgroup $G\le S_n$ fixes some tournament: choose one orientation on each orbit of unordered pairs and transport it by $G$; a pair stabilizer cannot exchange the endpoints, since $G$ has odd order. Conversely no tournament automorphism group has even order, since an involution would swap a pair. Thus for every odd prime $p$, the maximum $v_p(|\operatorname{Aut}(T)|)$ at order $n$ is $v_p(n!)$ (use a Sylow subgroup of $S_n$), and prime $p$ first occurs at $n=p$. At orders 23,24,25,26 the changes are respectively the first prime 23, one more factor of 3, two more factors of 5, and one more factor of 13. An *element* of order 25 first occurs at 25 (a permutation of order 25 needs a 25-cycle, and a $C_{25}$-invariant tournament exists), whereas a group order divisible by 25 already occurs at 10 via two independent five-cycles. These are distinct invariants. None is a statement about the supersingular/Ogg filtration.

There is a precise first-omitted-prime boundary for Burnside enumeration. In the ordinary sum, only all-odd cycle types can fix a tournament. If you keep all odd primes below an omitted prime $p$ as allowed cycle-length divisors, the truncation is exact before $n=p$, where the first omitted term is $2^{(p-1)/2}/p$ from a single $p$-cycle. It is a weighted term, not the number of discarded classes. In the converse-twisted sum, a returning unordered-pair orbit must fix the ordered pair after an even number of steps or swap it after an odd number. Inside a vertex cycle this requires every nontrivial length to be 2 modulo 4: generic pair returns use the full length, and opposite pairs swap after half that length. At most one singleton is allowed, and the stated conditions are sufficient also for pairs between cycles. Hence odd $p$ first enters the twisted sum at $2p$, and the first omitted term is $2^p/(2p)$, with $p$ circular-distance orbits. With the first 25 primes retained, these boundaries are 101 and 202. In the merged count both Burnside pieces have the final factor $1/2$.

The omitted ordinary part also has a uniform relative error estimate. Relative to the identity term $2^{\binom n2}/n!$, its value is at most $q_n^p/(1-q_n)$, where $q_n=n2^{-(n-1)/3}<1$ for $n\ge12$. Indeed every omitted type moves at least $p$ vertices; for $m$ moved vertices every pair touching them lies in an orbit of size at least three. Its contribution to the loss of edge-orbit count is therefore at least $\frac23[m(n-m)+\binom m2]\ge m(n-1)/3$. There are at most $n^m$ permutations with that support size. Summing the resulting geometric bound proves the estimate. This complements the older [unit-tail/prime-tail compression](../../../07-reflections/the-burnside-unit-tail-is-soluble-and-the-next-obstruction-is-the-prime-tail-ladder-codex.md): it bounds the cost of a prime filter for **vertex counts**, not the actual metagraph edges or an exact post-threshold truncation.



### Joint length with both clocks and the prime sieve

The joint statement strengthens to arbitrary units in both clocks **and** the first-25 sieve at the same length. In the padded 25-block, cuts 1, 2, 5, 12 give $(s^*,e^*)=(84,556),(56,503),(152,506),(126,481)$, whose products run through all four units mod 5. Once the 23-clock fixes the parity of $n$, choose desired base $r=4\bmod72$ for even or $r=5\bmod72$ for odd. Pick the 25-cut so its unit target is reached at $n\equiv r\bmod4$; the four different prefactors modulo 5 permit exactly this. The 23-clock period and the 25-clock period have gcd 2, and compatibility with the base modulus 72 requires only their prescribed parity and mod-4 class. Their possible extra overlaps with the sieve side registers are the primes 11 and 5, whose rules 67 and 61 are both disabled at bases 4 and 5. Choose the remaining side registers to survive, then apply CRT. This shows the two exact tournament clocks and the prime-sieve avoidance can coexist for any specified unit targets, without asserting a primality result or putting the two prescribed residues in the same tournament.


## Two six-prime filters that can actually be compared

The [level-one supersingular set](https://sigma-journal.com/2019/007/) is
\`{2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}\`; all are among the first
25 primes. The Baby Monster prime support removes exactly
\`{29,41,59,71}\`. The active divisibility sieve above meets the
level-one set in the six primes \`{5,13,17,19,29,59}\`; after redundancy
elimination only 29 drops from this intersection. The odd mirror
bridge in the companion note has norm form $x^2+2y^2$, whose
**necessary** prime-output congruence filter meets the level-one set
in the different six primes \`{3,11,17,19,41,59}\`.
The two sixes intersect in \`{17,19,59}\`. Neither is a group invariant
or a sufficiency claim for exact Hamilton counts. In particular
23 is present at both moonshine levels, inactive in the exponent
sieve, excluded by the norm congruence, absent in the self-converse
sector, and realized by an ordinary tournament. The inputs and
quantifiers, not just the resulting small sets, explain the contrast.
