# A small exact floor for a prescribed ternary overlap

A proper frame on sorted anchors $a<b<c<d$ has $\beta=c-b-1$ triply eligible spine sites. Let $F_b$ be the least ordinary Hamiltonian count of a two-hub tournament containing a proper frame with at least $b$ triply eligible sites, for $b\ge0$. This is a genuine whole-count lower envelope, not the framed weight.

**Core reduction.** Keep both hubs, the four anchors, and exactly $b$ sites between the two inner anchors. The induced word has two outgoing-matched left anchors, $b$ free middle digits, and two incoming-matched right anchors. Each matching has seven possible ordered digit pairs, giving exactly $49\cdot4^b$ normalized words. Every such core has the required frame and is strong: the hubs reach each other via the hub arc and a path from the second hub through its left anchor, the first hub's right anchor, and the first hub. The first spine site is reached by a hub and the last beats a hub, so all spine vertices are in that SCC. Hence the core minimum is attained. Ordinary counts are monotone under induced expansion by the insertion injection, so no larger host can do better.

| $b$ | exact $F_b$ | one minimizing spine word |
|---:|---:|---|
| 0 | 17 | 2302 |
| 1 | 41 | 23202 |
| 2 | 105 | 232202 |
| 3 | 281 | 2322202 |
| 4 | 777 | 23222202 |
| 5 | 1981 | 210023312 |
| 6 | 4909 | 2100223312 |
| 7 | 12301 | 21000233312 |
| 8 | 31693 | 210002233312 |
| 9 | 82957 | 2100002333312 |
| 10 | 221341 | 21000022333312 |
| 11 | 600589 | 210000023333312 |

These are exact finite depth values but **global at those depths**, since the host order was removed from the minimization. The 36-state replay counts the literal cores and verifies the number of candidates in each row. The depth-one and depth-two starting floors agree with [Session333](https://github.com/eliottcassidy2000/math/blob/e66cab99814e6af4b22dc10e06ccbb37a2220f33/05-knowledge/results/session333_smooth_sums_and_trijection_20261114/TWO_HUB_FRAME_DEBT.md).

## The optimizer changes its core after depth four

The earlier homogeneous optimum $23\,2^b\,02$ (exponent on the repeated digit) stays valid and gives $1+8\cdot2^b+8\cdot3^b$. It stops winning at $b=5$. For $b=5,\ldots,11$ the minimizing words have the form $21\,0^p2^q3^r\,12$. There is an all-depth theorem **within this block family**: for each $b=p+q+r\ge3$ with all lengths positive, its best choice is
$$
 p=r=\lfloor(b-1)/2\rfloor,\qquad
 q=b-2\lfloor(b-1)/2\rfloor.
$$
Write $A=2^{p+r}$, $B=3^{p+r}$, $Q=2^q$, $R=3^q$, $S=3^p+3^r$ and $U=2^p3^r+3^p2^r$. The whole-count formula is
$$
 3-2Q+2QS+(4+8Q+2R)A+(4Q+2R)U+2RB.
$$
The three homogeneous lengths each have basis $1,2^s,3^s$, by the occupied-run surjection count. The replay checks the $3^3$ basis grid, proving this identity at all positive lengths. For fixed $q$ and $p+r$, the positive terms in $S$ and $U$ are minimized by balancing $p,r$. If $q\ge3$, replace $(p,q,r)$ by $(p+1,q-2,r+1)$. The change is
$$
 \tfrac{Q}{2}(3-S)+(12-10R/9)A+(2Q-2R/3)U<0.
$$
Thus leave only the one or two central 2s required by parity. This proves an exact upper envelope by a balanced two-ended layout, and explains its internal minimizer. **It does not prove that arbitrary core words keep losing to this layout beyond $b=11$.**

## Feed the small cores into the prescribed target

For 385 and 573, the next floor 777 forbids four triple sites. To test depth exactly three, start with the 3136 cores and discard any with ordinary count above the target. Reinsert sites in any of the four gaps *other than* the central gap. Retain the anchor indices and use monotonicity for pruning. For strong witnesses use the existing proven spine cap $\lfloor\log_2N\rfloor+3$; a positive-depth frame cannot depend on neutral padding. For 385 there are only 32 viable cores and **480 retained anchored states**, ending in no strong exact hit. For 573 there are 772 viable cores and 16542 retained states, including the original exact-depth-three witness. This is a substantially smaller 385 certificate than scanning the entire two-hub word tree. It is still a specific target test, not a guarantee of general coverage.

Build **g++ -O3 -std=c++17 core_census.cpp -o core_census** and run **./core_census 11** for the table. Run **python replay_core.py** for the matched-ends formula and the anchored target certificates.
