# Geometry, exact inverses, and the prime-to-\(7p\) arc problem

**PROVED** for the structural statements. **FINITE-EXACT** for the literal certificates and the order bound. The prime/\(7p\) existence reduction is [the prior two-hole packet](../chanpak_two_holes_atoms_20261123/README.md); the full spectrum remains open.

## 1. What transfers from Hadwiger–Nelson, and what does not

The earlier [unit-port and spindle work](../chanpak_hadwiger_planar_folds_20261122/SPINDLE-CHARTS.md) gives exact Euclidean representatives of some physical flip subgraphs, while the [unit-port packet](../chanpak_hadwiger_planar_folds_20261122/UNIT-PORTS.md) separates soft coloring bias from hard color forcing. The same warning is stronger for the present targets. Let vertices be marked tournaments, and join two if their physical single arc reversal either changes the height by two or pairs an eligible prime \(p\) with \(7p\). **This entire mixed graph is bipartite**, including any version on unlabelled classes. Both types reverse the exact sign \(\chi_4(H)\): in the first case the difference of heights is \(2\bmod4\), and in the second \(7\equiv-1\bmod4\). Thus neither union nor gluing of these edges can carry a nonbipartite Hadwiger–Nelson obstruction such as the Moser spindle. The spindle inherited from the *full* flip graph necessarily uses other kinds of edges. This is a limitation of the selected observable, not an assertion about the chromatic number of the unit-distance plane.

The transferable idea is to specify marked ports and their allowed transitions, not infer them from geometry alone. Indeed **every finite bijection of distinct planar point sets over a characteristic-zero field extends to an automorphism of Jacobian determinant one**. Choose common determinant-one linear coordinates so the source abscissae \(x_i\) are distinct and target ordinates \(v_i\) are distinct (only finitely many directions are forbidden). To send \((x_i,y_i)\) to \((u_i,v_i)\), interpolate \(P(x_i)=v_i-y_i\), \(Q(v_i)=u_i-x_i\), then use the vertical shear by \(P\) and horizontal shear by \(Q\). Inverse shears subtract the same polynomials in reverse order. For \(N\) points the map and inverse have degree at most \(\max(1,N-1)^2\). This is not the planar Jacobian conjecture: we deliberately build tame automorphisms and do not assert that general constant-Jacobian maps are invertible. A finite fitted Jacobian-one chart is not evidence that it preserves any unprescribed distance, physical arc, or integral contact. Compare the previous [bounded-degree replicated rigidity test](../chanpak_hadwiger_planar_folds_20261122/SPINDLE-CHARTS.md#3-a-finite-physical-rigidity-test-for-bounded-degree).

## 2. What transfers from the planar inverse work

The earlier [marked Jacobian bracket](../chanpak_hadwiger_planar_folds_20261122/MARKED-JACOBIAN.md) distinguishes formal prime-adic inversion from an actual Boolean contact. The new [dyadic inverse charts](DYADIC-PORTS.md) remove another ambiguity: they invert two selected integer reads on the **whole specified chart**, then still enforce that the widths are powers of two. The two-block chart has a unique positive inverse and nowhere-singular Jacobian on its physical positive quadrant, but its count is **never divisible by seven**. This is a genuine example of a global inverse failing for a target family, not just a singularity at the wrong point.

A three-run repair gives a constant determinant \(4y^2-1\) on its two-port map when the middle width \(y\) is fixed. It meets \(7\mid H\) in exactly two exponent residue classes. The determinant is a unit mod seven in both classes; nevertheless each designated flip leaves the seven lane. On 18 vertices the same strong tournament of height \(66073=7\cdot9439\) has two distinct one-arc neighbors of prime heights \(65089,66529\). Their primality and the marks are checked by exact replay. This is a physical prime–\(7p\)–prime fork, but its neighbors are **not** the cofactor prime \(9439\), and it is not a unit-height fork. The exact-factor problem calls for a separate equation.

## 3. The exact factor-seven equations for one SCC-crossing flip

Suppose \(H(T)=p\) with prime \(p\ge5,\ p\ne7\). By ordinal multiplication across SCCs, there is one nontrivial strong core \(C\) of height \(p\) and all other SCCs are singleton pads. Write \(E_C(u)\) for the number of Hamilton paths in \(C\) ending at \(u\), and
\[
U_C(u)=\sum_{S\ni u} H(C-S)\,\#\{\text{Hamilton paths in }C[S]\text{ starting at }u\}.
\]
Here \(H(\varnothing)=1\). Equivalently \(U_C(u)\) is an ordered split into two paths with the second path rooted at \(u\). We classify cross-component flips, up to taking the converse and ignoring untouched ordinal pads. The [rooted cross-component law](../chanpak_kaprekar_decimal_20261123/FIRST-JUMP-SQUARES.md#cross-component-gain-and-singleton-budget) yields exactly two possibilities:

| Layout before the flip of \(u\to v\) | Condition for \(p\to7p\) |
|---|---|
| \(C\to\{v\}\), \(u\in C\) | \(U_C(u)-E_C(u)=6p\) |
| \(C\to\{w\}\to\{v\}\), \(u\in C\) | \(U_C(u)=3p\) |

For adjacency the old paths using \(u\to v\) are precisely the \(E_C(u)\) paths removed, and the new paths using \(v\to u\) are counted by \(U_C(u)\). With \(s\ge1\) intervening singleton vertices, none of the old paths is removed, and the new ones number \(2^s U_C(u)\): each middle vertex is freely assigned to either of the two traversals. Since \(v_2(6p)=1\), this forces \(s=1\) and the second equation.

For completeness, both endpoints cannot be singleton. On the affected interval, if there are \(m\ge2\) middle vertices their ordered two-slot count is
\[
Q_2(X)=\sum_{S\subseteq X}H(X[S])H(X[X\setminus S])\equiv0\bmod4.
\]
Indeed complementary summands pair into twice odd products and the number of pairs is even. The ratio therefore cannot be seven, which would require an increase \(6H(X)\equiv2\bmod4\); any exterior ordinal factor cancels from the ratio. The cases of zero or one middle vertex have ratios one or three. These arguments also cover a flip missing the prime core entirely. A flip **internal** to the prime core is a separate architecture and is not excluded by the rooted equations.

## 4. Literal first seams and the sharp order bound

Through strong core order seven neither rooted equation works for any eligible prime. At core order eight, the adjacent chart gives precisely \(p=163\); the one-middle chart gives precisely \(73,179,191,233,307,389\). The following two row sets show the mechanisms; append the singleton(s), orient the new ordinal arcs forward, then reverse \((u,\mathrm{last})\).

| \(p\) | Core outgoing rows | \(u\) | Pads | Check |
|---:|---|---:|---:|---|
| 163 | \((222,220,104,16,36,75,152,60)\) | 0 | 1 | \(983-5=6\cdot163\) |
| 73 | \((158,12,24,16,34,79,159,62)\) | 5 | 2 | \(219=3\cdot73\) |

The first gives \(163\to1141\) on nine vertices, the second \(73\to511\) on ten. Both upper tournaments are strong.

**Nine is the global minimum order** for an eligible prime-to-its-own-\(7p\) reversal, not merely the minimum for a crossing flip. The replay exhausts every Hamilton-framed *upper* tournament of order 4 through 8, retains those with height \(7p\), and tests every arc including the frame arcs. This is complete up to relabelling because every tournament has a Hamilton path. The scan tests 4,442,879 marked arcs from eligible upper frames (counting multiple frames of one isomorphism class); none leads to its own \(p\). Orders at most three cannot reach even 35. The nine-vertex witness attains the bound. The census does not assert uniqueness among all possible internal flips on nine vertices.

There is also a different negative, **valid at every tournament order**: no \(p\to7p\) edge exists for \(p=5,11,13,17,19,23,29,31,37,41,43\). For \(p<45\) [Busch's exact minimum at strong order eight](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf) bounds its prime core by seven vertices; cross-component flips are excluded by the rooted census and the formula above. Internal flips can therefore be checked only within those cores, regardless of any pads. The [small-core replay](small-prime-ports.js) checks every internal edit as well. This rules out proposing a universal *one-flip* \(p\mapsto7p\) for the full atom reduction, even though those two endpoints could have separate witnesses.

## Replay and next target

- Compile [sevenfold-census.cpp](sevenfold-census.cpp) with `c++ -O3 -std=c++17`; it checks the unrestricted order bound and the rooted census.
- Run `node small-prime-ports.js` for the all-order obstruction at the first eleven eligible primes and the literal rows.
- Run `python3 dyadic-ports-check.py` for the two inverse charts, all 27 modular exponent phases, and the same-host prime fork.

What survives from both neighboring problems is a stricter specification: build a small exact port relation, and specify which integer targets lie in its physical image. Proper coloring and an invertible chart cannot manufacture the factor-seven endpoints. Conversely the finite exact witnesses supply real marked transitions without claiming either a renewable rule or the Chan–Pak conjecture.
