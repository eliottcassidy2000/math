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

## 4. The two unfilled Moore branches: independent repairs stop at three

The marked partial board, with branches 13 and 33 unfilled, is separate from the 129-matching theorem for distinct complete charts. Keep the 50 reference contacts fixed and write $u(a)$ and $v(b)$ for the three collision counts as each of the two holes independently selects a contact. For two choices $a,b$, the actual vector is

    score(a,b) = u(a) + v(b) + W(a,b),
    W_j(a,b) = [the two chosen contacts coincide at frozen center j].

This gives an exact interaction: the vector mixed difference of a four-corner move is entirely that of $W$. Equivalently its sum is the dot product of the two signed changes of occupied center-slot indicator vectors. A zero total interaction is weaker than zero interaction at each of the three centers.

On all $32\times33=1,056$ pairs of physical choices, 985 have $W=0$ and 71 meet in one center; none meets in two. There are 261,888 four-corner choices. Of these, 196,849 commute at **each** center, as opposed to 199,715 with zero interaction in the total. For 1,563 of the former, both changes have absolute size one in the total; 56 such squares reach the floor. One example changes sheets $3\to16$ at branch 13 and $4\to19$ at branch 33, giving

    (1,2,2), (1,2,1), (1,2,1), (1,2,0).

Thus the two changes really are independent there, but the minimum is still three: the separate best old-owner collision counts are one and two, and the two holes can achieve both simultaneously without hitting each other at all. Dropping below three therefore requires changing a reference contact, not simply making the two hole choices more independent.

At each of the three minimum pairs $(16,19),(16,20),(29,16)$, the collision-owner table records what happens if a single owner moves away from both hole contacts. Every such repair hits at least two *other* original owners, and branch 38 needs at least three at $(16,19)$. For instance, in the $(16,20)$ pair, the unique least-disruptive repair for owner branch 27 selects sheet 3 there and hits owners 7 and 22. This gives concrete next branches to consider; it does not prove a global exclusion for longer simultaneous repairs.

In the original complete centerwise search, allowing up to 12 reference edits still excludes all three **concentrated** two-square caps $(2,0,0),(0,2,0),(0,0,2)$, with no node-limit truncation. The three **split** caps $(1,1,0),(1,0,1),(0,1,1)$ have verified exclusions through 11 edits; do not infer that their twelfth layer is excluded. The [depth-12 receipt](concentrated_cap12.out) records the completed cases.

The alternative MILP (squared loads at most two represented by binary overload flags) hit its 130-second limit with no feasible incumbent; that is not an exclusion. The exploratory center-preserving cycle walk likewise does not certify infeasibility. The exact result above uses actual contact rows and the complete centerwise search only.

To recreate the input, run `python nested_core_check.py --board > board.txt` in the [pinned Moore chart folder](../moore57_quadratic_tensors_20261125/). Then replay with `python two_hole_interaction.py board.txt`. See [the interaction receipt](two_hole_interaction.out).
