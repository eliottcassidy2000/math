# Square-sum rails, the first 25 prime clocks, and the order-56 Latin grid

2026-11-26. This follows [the coupled Latin-square model](../moore57_coupled_latin_matrices_20261126/README.md), combining the old [equal-fifths rail and first-25-prime wheel](../first25_spread_target_returns_20261116/QUINTILES.md) with the [prime-value versus prime-rank ledger](../prime_lane_ports_five_groups_20261102/prime_groups.md). No degree-57 graph is asserted.

## 1. A prime clock for the actual square reflections

The old rail frame on $[1,5m]$ joins blocks 1 to 3 and 3 to 5 by the square reflections $x\mapsto 3m+1-x$ and $x\mapsto7m+1-x$. They compose to translation by $4m$. The widths are indexed by the Pell orbit with $m_0=0,m_1=5$, and satisfy the **integer** recurrence

    m_(n+2) = 23 m_(n+1) − m_n + 5.

Indeed the rail matrix has eigenvalues $\rho,\rho^{-1}$, with $\rho=(5+\sqrt{21})/2$. Its symmetric square has eigenvalues $1,\lambda,\lambda^{-1}$ for $\lambda=\rho^2$. The constant term of $m_n$ is $-5/21$: writing $a_n=\alpha\rho^n+\beta\rho^{-n}$ in the seed $(a,b)=(1,1)$ gives $\alpha\beta=1/7$, and $m_n=(a_n^2-1)/3$. Hence $u_n=21m_n+5$ obeys the trace-23 homogeneous recurrence. The endpoints $u_{-1}=u_0=5$ give

    u_n = 5(λ^(n+1)+λ^(−n))/(1+λ).

For any prime $p\nmid 2\cdot3\cdot5\cdot7$, reduce in the splitting field of $\lambda^2-23\lambda+1$. The test $p\mid m_n$ factors into

    (λ^n−1)(λ^(n+1)−1)=0.

Thus the tag is frozen exactly at **$n\equiv0,-1\pmod{\operatorname{ord}(\lambda)}$**. The choice of reciprocal root does not affect the answer. The order divides $(p-\left(\frac{21}{p}\right))/2$, because $\rho$ lies in the norm-one group in the inert case and $\lambda$ is its square. This is the split/inert character of $\mathbb Q(\sqrt{21})$, not the separate mod-four coloring of the prime square-sum graph. The exceptional prime clocks are direct: mod 2 and 3 the period is 3, mod 7 it is 7, each with zeros at 0 and −1; mod 5 all widths vanish. The replay gives every period for the **first 25 primes**, a shifted analogue of the archived multiplicative exponent wheel, not a primality test.

### Transfer to the physical 56 column labels

Modulo 56 the rail transport is $x\mapsto x+4m_n$; the exact recurrence mod 14 has period 21. Its order and the exceptions in the [actual order-56 route Latin square](../moore57_coupled_latin_matrices_20261126/README.md) are:

| Order of $+4m_n$ | Width indices $n$ mod 21 | Latin routes commuting with this transport | Exceptions |
|--:|:--|--:|:--|
| 1 | 0,6,14,20 | 56 | none |
| 2 | 7,13 | 54 | 27,49 |
| 7 | 2,3,5,8,9,11,12,15,17,18 | 54 | 10,45 |
| 14 | 1,4,10,16,19 | 52 | 10,27,45,49 |

All generators of the same translation subgroup commute with exactly the same route maps. The order-14 defect splits into disjoint order-2 and order-7 exceptions, rather than four indistinguishable errors. This is a literal comparison with the two completed rows of the existing chart. It does not assume that the other rows follow these periods. Also the two rail reflections themselves are fixed-point-free on cyclic columns only for even $m_n$ (both sum constants are then odd). In those cases their product has order 1 or 7. Width 5 gives order 14 as a transport but **not** two legal full-grid reflected matchings; the original rails only used restricted intervals.

## 2. How many whole Latin blocks can a square-sum graph supply?

Use a bipartite support matrix on two copies of the physical labels, allowing $x\to y$ when $x\ne y$ and $x+y$ is a square. A bank of $k$ pairwise disjoint whole permutation blocks is equivalent to a bipartite $k$-factor. One direction is addition of the permutation matrices; the other repeatedly extracts a perfect matching from the regular bipartite graph. Unit capacities on the physical supports and capacity $k$ at each column label give an exact flow test. The replay records both successful permutations and failing min-cuts.

| Labels | Square edges | Maximum whole supported blocks | First failing flow / demand |
|:--|--:|--:|:--|
| Integers 1–25 | 32 | **1** | 49 / 50 |
| First 25 prime **values**, through 97 | 24 | **0** | 22 / 25 |
| Fifteen Ogg-prime values | 14 | **0** | 12 / 15 |
| Integers 1–56 | 111 | **2** | 167 / 168 |

The prime-value failure is already forced by the isolated labels 37, 43 and 67; the Ogg obstruction is the old 6/9 bipartition imbalance. In contrast the integer 25-prefix supports an actual whole permutation. This is a Latin-block consequence of the old warning that prime rank cannot be replaced by prime value. For the full grid, any row of matched blocks has disjoint supports across its 55 neighbors. Thus **at most two** blocks per row can use only square sums under any common labeling by integers 1–56.

The two square blocks can actually be **fixed-point-free involutions**. The replay gives their union as an even 2-factor with cycles of lengths 10 and 46; alternate each cycle's edges. Their product has cycles $5,5,23,23$, hence order **115**, which is not a column translation order. Their opposite-parity entry counts are 24 and 40 of 56. Thus stitching different square levels really can give two parity-breaking blocks, but not four. The earlier [four-defect theorem](../moore57_coupled_latin_matrices_20261126/README.md#4-four-modular-defects-per-row-not-just-three) now implies that **every full Moore row has at least two parity-breaking matchings that contain a nonsquare contact**.

## 3. Exact loss under rigid rail symmetry

There is exactly one cyclic reflection on labels 1–56 that uses only true square sums: $y=25-x\pmod{56}$, using ordinary sums 25 and 81. For a wrapping reflection the two possible sums differ by 56. Within the available interval the only pair of squares at that distance is $(25,81)$; the other nonzero solution $(169,225)$ is too high. An unwrapped reflection would have sum 57. This one mixed-level reflection is still parity-homogeneous, so provides no modular defect.

More quantitatively, if a permutation commutes with shift by $g\in\{4,8\}$, it permutes the $g$ orbits and on each orbit selects one relative cyclic shift. For each ordered pair of orbits score the best relative shift by its count of true square contacts; exclude the identity shift when the orbit maps to itself. An exact assignment DP over the orbit permutations gives:

| Required physical commutation | Maximum square contacts in a deranged block | Minimum nonsquare contacts |
|:--|--:|--:|
| shift 4 (order 14) | **16** | 40 |
| shift 8 (order 7) | **18** | 38 |

So imposing the Pell rail symmetry on a Moore matching is stronger than just restricting its period: it forbids a completely square-supported matching and forces many nonsquare slots. In particular the explicit involutory pair above cannot have both maps commute with shift 4 or with shift 8: its product contains the prime 23 in its order, while the centralizers act on four 14-cycles, or eight 7-cycles, and contain no element of order divisible by 23. This separates actual involutory square ports from a purely periodic imitation.

The two complete physical stars in the current chart have at most 13 and 12 square-sum entries in any one block under the standard labels. Each has 222 in total, as it must: every vertex sees every other column once. No claim is made that a common relabeling improves this to a usable global construction.

Replay: `python rails_and_wheel.py` with the pinned small inputs from the [marked chart directory](../moore57_quadratic_tensors_20261125/) on the Python path and in the working directory.
