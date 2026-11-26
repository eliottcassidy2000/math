# Binary cover inversion and parity at the 129-matching frontier

2026-11-26. Follows [nested unit forks](NESTED-FORK-INVERSE.md) and [pooled Moore capacity](POOLED-QUARTET-CAPACITY.md). Both results use pinned, exact replays. Neither settles a main existence question.

## 1. Binary profiles through order 23, with arbitrary ordered prefixes

For the converse binary core, let $v$ be its $s$-bit word, $x_i$ its bits from low to high, and $n=s+1$ its vertex count. A positive even $v$ has root counts one in each direction. Its ordinary path count is $F(v)$, its starting-second cover count is $v+1$, and its ending-first cover count is one plus the $s$-bit reverse complement of $v$. The last assertion follows from the standard binary starting-second cover count after taking the converse and reversing the spine order.

Let $Q(t)$ count ordered partitions of the core into $t$ path slots, with an empty slot allowed. All directed odd cycles use the root, so no two disjoint slots can each contain a directed cycle. In the [Grinberg–Stanley expansion](https://arxiv.org/abs/2307.05569), each odd cycle contributes twice, the residual vertices can occupy arbitrary slots, and the cycle itself must occupy one common slot. The extreme spine sites of a cycle are precisely an inverted pair $i<j$ with $x_i=0,x_j=1$. In its decreasing spine run one can choose any even subset of the $d=j-i-1$ interior sites. Summing the resulting binomial gives

    Q(t) = t^n + Σ_(i<j; x_i=0,x_j=1) t^(n−2−d) ((t+1)^d+(t−1)^d).

For $t=2$ write $U_a=\sum_{i<j:x_i=0}2^i a^{j-i-1}$, with $a=1,3$, before appending bit $x_j$. The update for either bit doubles $Q(2)$, then adds $2(U_1+U_3)$ if the new bit is one. The next $U_a$ is $aU_a+2^j[x_j=0]$. The analogous inverted-pair sum gives $F(v)$, so the scan does constant work per word before the inverse filters. The replay independently checks all marked counts and cover values by subset convolution on 128 literal cores, including the singleton.

In the nested inverse write $A,B$ for the effective coefficients and $M=2^{60}(2^{61}-1)$ for the midpoint. Reject $\gcd(A,B)>1$ since the odd common divisor would have to divide 2. With $k=\lfloor\log_2 A\rfloor$, the earlier square-root window and $u\le w<2u$ give

    2^((k+119)/2−2) < u < 2^((k+121)/2+1),     u=B·2^j.

Indeed $2^{120}<M+D<2^{121}$ for both offsets throughout this bank (the offsets are less than $2^{30}$); the window radius $12A$ is less than half its square-root center. Use exact integer bit lengths to find a conservative range for $j$, which is the spoke index plus any prefix length. The exact unpadded and padded inverses follow after the leading congruence test.

**No crossing** at exponent 61 among all **2,097,151** positive even binary addresses with 2–22 spine bits, the three outer cores, both outer orientations and every ordered prefix length. In the replay, 15 coarse profile/scale candidates lead to four exact padded congruence lanes but no compatible prefix/spoke split or suffix read. The previously checked arbitrary uniquely terminal cores through order 12 are a separate bank.

## 2. At 129 changed matchings: 97 candidates become 82

Here $(h,s,l,p)$ counts high branches, closed branches, their low branches, and active sheets per closed branch. These are **necessary tuples for two distinct full completions of the actual frozen stars**, not completions.

**Four-branch obstruction.** Four closed lows cannot span five or six changed links. Six would form the $K_4$ already excluded by the quartet-shadow test. For five, call the missing edge $ab$ and the other vertices $c,d$. Then $Q_c,Q_d$ both contain $a,b$, while $d$ belongs to $Q_c$ and not $Q_d$, and vice versa for $c$. The remaining quartets $Q_a,Q_b$ contain $c,d$ and omit $a,b$. They must all be distinct by the oriented-profile uniqueness at a common closed target. Enumerate these projections with multiplicity at least two, necessary for the remaining $p\ge4$. The catalogue gives 918 assignments on 219 label sets. At every target, incoming old quartets must admit at least one common residue class modulo four in the branch-56 identity gauge. Every assignment has an empty intersection at two or more targets: the counts for two, three and four empty intersections are 96, 262 and 560. Thus four closed lows have deficit at least eight. The same holds for five: each four-subset has at most four edges, so the five-set has at most six. The slack-six rows $(40,5)$ and $(41,4,4)$ lose eleven tuples in total.

**Triangle square-root obstruction.** A changed low–low matching on active sheets commutes with the half-turn $\tau$ when expressed in the two branch coordinates: both ends describe the changed matching as the old matching followed at its target by $\tau$. On a triangle with all three links changed, compose the old reverse matchings into $L$. It cannot fix an antipodal pair, since on that pair it would either fix the physical sheets or interchange them, and three half-turns would create a fixed sheet in one chart or the other, producing a triangle.

Each oriented trade gives a map between any two target sheets through its low source. Taking the catalogue routes $a\to c\to b$, $b\to a\to c$, $c\to b\to a$ in order yields $D=L^2$ on the physical sheets and their antipodal pairs. Replaying 8,914 candidate assignments gives:

| Closed pair orbits of $D$ | Physical behavior | Assignments |
|:--|:--|--:|
| None | None | 8,617 |
| Seven singletons | All seven physical pairs swapped | 122 |
| Six singletons | All six physical pairs swapped | 3 |
| One orbit of size seven | Two physical 7-cycles | 104 |
| One orbit of size seven | One physical 14-cycle | 68 |

A selected set on singleton pair orbits needs a fixed-point-free involution on those pairs, allowing only an even number of pairs; a selected seven-orbit must use the whole orbit. The 68 physical 14-cycles have no square root on those sheets, since even cycles of a permutation square come in equal-length pairs. Thus only $p\in\{4,8,12,14\}$ can pass this triangle test. In $(41,4,3)$ and $(42,3,3)$ slack six forces such a triangle, so $p=6,10$ are impossible in both. This leaves **82 necessary tuples in ten rows**; the conditional 129 lower bound is unchanged.

The same replay also checks that all 12,684 physical target-sheet pairs across the 2,114 oriented profiles are distinct. This explains why the repeated target contact test can act directly on profiles.

Replay: `node binary_covers.js` and `node moore_parity_frontier.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt`.
