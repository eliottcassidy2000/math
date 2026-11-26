# Independent responses with prescribed amplitudes

2026-11-26. Continues the [commuting unit forks](../overnight_chanpak_moore_20261126/COMMUTING-UNIT-FORKS.md) and [nested inverse](../overnight_chanpak_moore_20261126/NESTED-FORK-INVERSE.md). Independence here means the four literal orientations have counts $N,N+a,N+b,N+a+b$. All constructed sizes are finite and the two toggled arcs are vertex-disjoint.

## 1. Why two binary spokes are not independent

Write $v=\sum_i x_i2^i$. Expanding the binary recurrence gives

    F(v) = 1 + Σ_(i<j) (2^(j−i−1)+[j=i+1])(1−x_i)x_j.

The count after joining to a fixed rooted core is affine in $F(v)$, $v+1$, and the lowest bit. The latter two have zero mixed difference. Thus any two distinct spoke toggles on the same binary arm have a nonzero mixed difference, of absolute value $h(2^{j-i-1}+[j=i+1])$, where $h$ is the coefficient of $F$. This is independent of the other bits. Observing two opposite single flips on one arm is not the same as having independent responses.

On an arbitrary nonbinary rooted core, two isolated root edits could in principle coexist. A separate exhaustive unlabeled census through order nine found none on disjoint arcs when both are required to preserve $(H,R)$ at all four corners and add their nonzero $E$ increments. At orders seven, eight and nine it tested 2, 11 and 74 disjoint candidate pairs, respectively; there were none at smaller orders. This does not exclude three independent contacts by a richer seam.

## 2. Any one even amplitude, independently of an arithmetic progression of others

To realize any $g\ge1$ as the number of spanning paths ending at a root, use an ordered tournament on $m$ vertices with just the long arc from last to first reversed. The endpoint weights, in the same vertex order, are

    1, 1, 2, 4, …, 2^(m−3), 1.

The increasing path ends at the last vertex. Paths using the reversed arc split the interior sites into a prefix before last and a suffix after first. For endpoint $j\in\{1,\ldots,m-2\}$, choosing the suffix below $j$ gives $2^{j-1}$ paths; the first vertex itself has one. Attach the marked root so it receives exactly from vertices whose endpoint weights sum to $g$. Binary expansion gives such a subset for $m=3+\lfloor\log_2g\rfloor$. The resulting inner core has exactly $g$ marked terminal paths.

With this inner core, the nested outer count for odd $c$ still has the form $A F(c-1)+B c-D$. The only change from the earlier derivation is $D=H e_+ +Rv+Eg$, and $A=Hh$, $B=Hr+RQ$ as before. Changing the isolated endpoint count $E$ of the outer six-vertex left core now changes the total by exactly $2g$ in absolute value. An inner binary spoke at bit $i\ge1$ leaves the right terminal count $g$ unchanged. It gives any prescribed even shift $\sigma$ satisfying $B2^i\equiv\sigma\pmod A$, by putting

    β=(B·2^i−σ)/A+Z+(Z mod 2),     q=β or β−1

as in the previous suffix construction. Since $A$ is odd, each eligible index repeats periodically, so choose it sufficiently large to get positive $q$. For a fixed starting index, the eligible even $\sigma$ form one full arithmetic progression of modulus $2A$. The two responses commute and their signs can be set independently by choosing the starting corner.

In addition, the literal core bank gives an exact construction of independent magnitudes **$(2g,2)$ for every $1\le g\le512$**. The replay searches the endpoint-weight representations on the smallest spine length and up to two extra sites. Two values, 260 and 426, use a third extra site. It verifies each selected terminal profile and the unit congruence; the recurring indices yield arbitrarily large heights for each pair. For example $g=7$ uses inner masks $(46,28,56,48,1,18)$ rooted at 5, outer $(H,R)=(25,28)$, $c=1681$, $i=3$. At one starting corner the four counts are **20,545,437; 20,545,451; 20,545,435; 20,545,449**. The changes are $+14$ and $-2$ with no interaction.

## 3. With one unit fixed, nearly every second magnitude is available

For a single binary right arm, take one of seven isolated left cores. For each $(H,R)$ here the set $\{\pm R2^i\pmod H:i\ge1\}$ is exactly the unit group; the small multiplicative orbits and the physical root profiles are replayed.

| $(H,R)$ | Order of $2\bmod H$ | Size with both signs |
|:--|--:|--:|
| $(23,24)$ | 11 | 22 |
| $(25,28)$ | 20 | 20 |
| $(29,31)$ | 28 | 28 |
| $(59,87)$ | 58 | 58 |
| $(71,72)$ | 35 | 70 |
| $(103,129)$ | 51 | 102 |
| $(143,173)$ | 60 | 120 |

This yields independent **$(2,2m)$** whenever $m$ is coprime to at least one of these seven moduli. Write $D=5\cdot23\cdot29\cdot59\cdot71\cdot103=1,438,942,445$. The only integers not reached by this specific criterion satisfy $D\mid m$ and $11\mid m$ or $13\mid m$. In particular every $1\le m<11D$ is reached, and the fraction not certified is $23/(143D)$. This is a theorem about prescribed **responses**, not about the set of their absolute central path counts.

Replays: `node independent_responses.js`; compile `double_arc_census.cpp` with `g++ -O3 -std=c++17` and pass it order `9`.
