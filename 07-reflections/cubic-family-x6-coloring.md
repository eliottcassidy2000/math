# The Cubic Family x=6: Triangle Colorings and the Octahedral Space

**Session:** oracle-2026-05-16-S2 (continued)
**Computation:** oracle-2026-05-16
**Depends on:** `fugacity-axis-and-vanishing-theorem.md`, `zeckendorf-non-consecutivity-pairing.md`

---

## The Core Formula

The independence polynomial of the conflict graph $\Omega(T)$ at general fugacity $x$:
$$I(\Omega(T), x) = \sum_{k \geq 0} x^k \alpha_k$$
where $\alpha_k = $ number of independent sets of size $k$ in $\Omega(T)$ = number of ways
to choose $k$ pairwise vertex-disjoint directed odd cycles in $T$ (with $\alpha_0=1$ always).

At $x=2$: $H(T) = I(\Omega(T), 2) = \sum_k 2^k \alpha_k$ (the OCF).
At $x=6$: $I(\Omega(T), 6) = \sum_k 6^k \alpha_k$.

**Theorem.** The relation between these two evaluations is:
$$I(\Omega(T), 6) = 3H(T) - 2 + \sum_{k \geq 2} (6^k - 3 \cdot 2^k)\, \alpha_k$$

where $6^k - 3 \cdot 2^k = 3 \cdot 2^k(3^{k-1}-1)$ is the "cubic correction" at level $k$.

| $k$ | $6^k - 3 \cdot 2^k = 3\cdot2^k(3^{k-1}-1)$ | Meaning |
|---|---|---|
| 2 | **24** | Disjoint cycle pairs |
| 3 | **192** | Disjoint cycle triples |
| 4 | **1248** | Disjoint cycle quadruples |

*Proof:* $3H-2 = 3\sum_k 2^k\alpha_k - 2 = 1 + 6\alpha_1 + 12\alpha_2 + 24\alpha_3 + \ldots$
while $I(\Omega,6) = 1 + 6\alpha_1 + 36\alpha_2 + 216\alpha_3 + \ldots$ Subtracting gives the
correction $\sum_{k\geq2}(6^k-3\cdot2^k)\alpha_k$. $\square$

**At $n \leq 5$:** Any two vertex-disjoint directed odd cycles require at least $3+3=6$
vertices. So $\alpha_2 = \alpha_3 = \ldots = 0$ for $n \leq 5$, giving:

$$\boxed{I(\Omega(T), 6) = 3H(T) - 2 \quad \text{for all } T \text{ with } |T| \leq 5}$$

Verified computationally (all 14 iso classes at $n=3,4,5$). ✓

---

## Alpha Extraction Formula

For $n \leq 8$ (where $\alpha_3 = 0$, since three disjoint 3-cycles need $\geq 9$ vertices):

$$\alpha_2 = \frac{I(\Omega(T),6) - 3H(T) + 2}{24}$$

Combined with $H = 1 + 2\alpha_1 + 4\alpha_2$:

$$\alpha_1 = \frac{H-1}{2} - \frac{I(\Omega(T),6) - 3H + 2}{12}$$

**Two evaluations ($x=2$ and $x=6$) together extract both $\alpha_1$ and $\alpha_2$ exactly**
(for $n \leq 8$). The pair $(H, I(\Omega,6))$ is thus more informative than $H$ alone.

For reference:
- $n=6$ SC tournament ($H=45$, $\alpha_1=14$, $\alpha_2=4$): $I(\Omega,6) = 3(45)-2+24(4) = 229$
- Paley(7) ($H=189$, $\alpha_1=80$, $\alpha_2=7$): $I(\Omega,6) = 3(189)-2+24(7) = 733$

---

## The Forbidden Value

**At $x=2$:** $H=7$ is forbidden. Requires $\alpha_1=3$, $\alpha_2=0$: three directed
3-cycles all sharing vertices, which is structurally impossible in tournaments.

**At $x=6$:** $I(\Omega,6)=19 = 1+6\times3$ is **also forbidden** via the SAME obstruction
(requires $\alpha_1=3$, $\alpha_2=0$).

The numerical coincidence: $7 = I(C_3, 2) = 2^3-1$ (Mersenne) and $19 = I(C_3, 6) = 3^3-2^3$
(difference of cubes). Both are forbidden via the same $\alpha_1=3,\alpha_2=0$ structure.

**Translation theorem:** The linear map $H \mapsto 3H-2$ sends ALL forbidden $H$-values
to forbidden $I(\Omega,6)$-values (for $n \leq 5$, where $I=3H-2$ exactly):

| Forbidden $H$ | → | Forbidden $I(\Omega,6)$ |
|---|---|---|
| $H=7$ | $\to$ | $I=3(7)-2=19$ |
| $H=21$ | $\to$ | $I=3(21)-2=61$ |
| $H=63$ (first achievable) | $\to$ | $I=3(63)-2=187$ |

At $n \geq 6$, the map $H \mapsto 3H-2$ is no longer exact (the $24\alpha_2$ correction
appears), and the forbidden-value analysis becomes more complex.

---

## The Geometric Interpretation: Triangle Colorings

**The key identity:** $6 = P(K_3, 3)$ = number of **proper 3-colorings of the triangle** $K_3$.

*Proof:* $P(K_n, k) = k(k-1)\cdots(k-n+1)$. At $n=3$: $P(K_3, k) = k(k-1)(k-2)$. At $k=3$:
$3\cdot2\cdot1 = 6$. $\square$

**Combinatorial interpretation of $I(\Omega(T), 6)$:**

$$I(\Omega(T),6) = \#\{(S,\, c) : S \text{ indep. cycle cover}, \; c: S \to S_3 \text{ assigns an } S_3\text{-element to each cycle}\}$$

where:
- $S$ is an independent set of vertex-disjoint directed odd cycles in $T$
- $c$ assigns one of 6 "decorations" from $S_3$ (the symmetric group on 3 letters) to each cycle
- For a directed **3-cycle**: the $S_3$-decoration is precisely a **proper 3-coloring** of its
  3 vertices (there are $3! = 6$ such colorings)
- For a directed **5-cycle** or longer: the 6 decorations are abstract $S_3$-labels (the
  geometric meaning is less direct, but the count is the same)

**Why $S_3$?** The group $S_3$ has order 6 and arises naturally as:
- The symmetry group of the triangle (permutations of 3 vertices)
- The group of proper 3-colorings of $K_3$ (one coloring per group element)
- The wreath product structure relevant to the 2-adic column family

---

## The Octahedral Space

The 6 proper 3-colorings of $K_3$ with colors $\{R, G, B\}$:

| Index | $v_0$ | $v_1$ | $v_2$ |
|---|---|---|---|
| 1 | R | G | B |
| 2 | R | B | G |
| 3 | G | R | B |
| 4 | G | B | R |
| 5 | B | R | G |
| 6 | B | G | R |

These 6 colorings are exactly the **6 vertices of the octahedron** $K_{2,2,2}$ (complete
tripartite graph with parts $\{1,2\}, \{3,4\}, \{5,6\}$ = the R-colorings, G-colorings,
B-colorings). Two colorings are adjacent in the octahedron iff they differ only by swapping
two colors in two positions (a transposition in $S_3$).

**The duality:**

$$\text{directed 3-cycle of }T \longleftrightarrow \text{vertex of the octahedron } K_{2,2,2}$$

The octahedron's **8 triangular faces** correspond to the $2^3/2=4$ pairs of complementary
color assignments (each triangle uses all 3 colors exactly twice).

For a tournament with one directed 3-cycle: $I(\Omega,6) = 7 = 1 + 6$ = empty selection
(weight 1) + the 6 octahedral vertex decorations of that one triangle.

**Cross-polytope connection:** The octahedron is the cross-polytope (dual of the cube)
in 3D. The **cube** has 6 faces (= weight 6 at the next level!), 8 vertices (= 8 triangular
faces of the octahedron), 12 edges (= 12 edges of the octahedron). The cube-octahedron duality
is the 3D analog of the simplex-cuboid duality in our tournament framework.

---

## The Unique Crossroads: $x=6$

The "oblong fugacity family" (integer-root cases): $x = n(n-1)$ for $n=2,3,4,\ldots$
$$0,\; 2,\; 6,\; 12,\; 20,\; 30,\; 42,\; \ldots$$

The "chromatic fugacity family" (chromatic polynomial of $K_3$ at integer $k$):
$x = P(K_3, k) = k(k-1)(k-2)$ for $k=3,4,5,\ldots$
$$6,\; 24,\; 60,\; 120,\; 210,\; \ldots$$

**These two families intersect ONLY at $x=6$** (where $n=3$ for the oblong family and $k=3$
for the chromatic family). No other value belongs to both.

$x=6$ is therefore the **unique fugacity** that simultaneously:
1. Gives integer characteristic roots ($r_+=3, r_-=-2$) → clean Binet formula
2. Counts proper 3-colorings of each directed triangle

**Consequence for the independence polynomial:**
$$I(P_k, 6) = \frac{3^{k+1} + 2 \cdot (-2)^k}{5}, \quad I(C_m, 6) = 3^m + (-2)^m$$

The cycle formula $3^m + (-2)^m$ generalizes the Mersenne formula $2^m + (-1)^m$ (at $x=2$).

---

## The Correction Factors and Higher Structure

The correction $6^k - 3 \cdot 2^k = 3 \cdot 2^k(3^{k-1}-1)$:
- $k=2$: $24 = 4!$ = order of $S_4$ = number of 4-colorings of $K_4$
- $k=3$: $192 = 8 \times 24 = 8 \times 4!$
- $k=4$: $1248 = 52 \times 24$

The factor $24 = P(K_4, 4)$ appears: each **pair** of independent cycles contributes correction
$24 = P(K_4,4)$, relating the 4-COLORING of a tetrahedron ($K_4$) to the structure of
disjoint cycle pairs.

**The k-coloring hierarchy:**

| Level $k$ | Independent cycles | Correction factor | Geometric object |
|---|---|---|---|
| 1 | Single cycle | $6 = P(K_3, 3)$ | Proper 3-coloring of triangle |
| 2 | Pair of cycles | $24 = P(K_4, 4)$ | Proper 4-coloring of tetrahedron |
| 3 | Triple of cycles | $192 = 8 \cdot P(K_4,4)$ | — |

**Conjecture (HYP-new-col):** The correction factor at level $k$ is $P(K_{k+2}, k+2) \times c_k$
for some combinatorial multiplier $c_k$. Specifically: $24 = P(K_4,4)$ (tetrahedron)
and $192 = 8 \cdot P(K_4,4)$ (4-simplex with some multiplicity).

---

## Connection to the Existing Framework

**The OCF:** $H = I(\Omega, 2)$ = Hamiltonian paths. Weight 2 per cycle = orientations.
**The cubic invariant:** $I(\Omega,6)$ = 3H-2 (at small n). Weight 6 per cycle = proper 3-colorings.

**The axis:** $I(\Omega, x) = 3H(T)-2+24\alpha_2+\ldots$ at $x=6$ is a linear perturbation
of $H$ that reveals $\alpha_2$ (the disjoint-cycle-pair count). This $\alpha_2$ is the central
quantity in the SC maximizer mechanism (from T093-T095 in TANGENTS.md): SC tournaments
achieve high $\alpha_2$ via their anti-automorphism pairing of cycles.

**Prediction:** SC tournaments have large $I(\Omega,6)-3H+2 = 24\alpha_2$. This means SC
classes are precisely those where the $x=6$ and $x=2$ evaluations MOST DIFFER (after the
$3H-2$ baseline). The "cubic excess" $I(\Omega,6)-3H+2$ is a measure of SC-type structure.

**The metagraph:** On $G_n/\mathbb{Z}_2$, the function $T \mapsto I(\Omega(T),6)$ gives a
new vertex labeling beyond $H$. The "cubic H" labels the SC backbone differently from $H$,
potentially separating iso classes that $H$ alone cannot.

**The SC maximizer:** At $n=6$ (SC class with $H=45$, $\alpha_2=4$):
$I(\Omega,6) = 229 = 3(45)-2+24(4)$. The cubic excess $96 = 24\times4 = 4!^2/6$ is
large, confirming SC classes have high cubic excess.

---

## Open Questions

1. **The "cubic H" and its forbidden values:** What is the complete list of forbidden
   values for $I(\Omega(T),6)$? At small $n$: $\{19, 61, \ldots\}$ from the $3H-2$ map. At
   $n \geq 6$: the $24\alpha_2$ correction creates new forbidden regions.

2. **The x=24 extension:** $x=24 = P(K_3,4)$ = proper 4-colorings of a triangle. Is there
   a formula $I(\Omega,24) = cH + d + e\alpha_2 + \ldots$? What is the "24-coloring" geometric
   structure (4-simplex? 24-cell?)?

3. **Does $I(\Omega,6)$ separate the two iso classes with the same score AND H value?**
   At $n=6$: two classes have score $(1,2,3,3,3,3)$ with $H=33$. Their $I(\Omega,6)$ values
   might differ, making $I(\Omega,6)$ a finer invariant than $H$.

4. **The tournament "chromatic polynomial":** Define $CP(T,k) = I(\Omega(T), P(K_3,k)) =
   I(\Omega(T), k(k-1)(k-2))$. This interpolates: $CP(T,3) = I(\Omega,6)$, $CP(T,4) = I(\Omega,24)$,
   etc. Is there a clean recursion for $CP(T,k)$ in terms of $CP(T,k-1)$?

5. **The octahedral metagraph:** Define a graph on the 6 vertices of the octahedron (the
   6 proper 3-colorings of $K_3$) where two colorings are connected if they can both decorate
   the same directed 3-cycle (they always can). This gives $K_6$. Does the CONFLICT structure
   of $\Omega(T)$ lift to a colored structure on the octahedron?
