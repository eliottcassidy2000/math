# Exact field escape and geometric spindle recursion

**Status: PROVED elementary geometric and counting statements; FINITE-EXACT authenticated controls. No six-colour equality gadget is constructed.** Leaving a known four-colour field, removing any prescribed finite list of residue-prime screens, and forcing higher chromatic number are distinct predicates.

The positive construction below gives, for every $m\ge1$, an explicit induced unit-distance graph with $6m+1$ vertices and $10m+1$ edges. It is four-critical, and deleting any edge gives an actual unit-separated three-colour equality fragment. The family includes graphs outside the Moser field. Nevertheless it remains four-colourable even when its coordinate field escapes every prime screen in any prescribed finite collection.

## 1. Inheritance and the retained interface

The inherited [THM-4558 field screen](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md), including its correction of the rotation-rank/Heegner roadmap, rules out treating more rotations as an automatic chromatic increase. The previous [Moser-field candidate packet](unit_equality_candidates_20261011e.md) makes that obstruction computational: every omitted actual unit edge is separated by a checked four-colouring. The [two-port count algebra](existence_count_amplifiers_20261011e.md) distinguishes zero different-port support from a small relative count.

The incoming [golden/taxicab gadgets](golden_taxicab_boundary_20261011.md), G1--G3, provide genuine exact *precoloured boundary* relations. Their private leaves can restrict a core to specified lists when the leaf assignment is supplied. They do not force that external assignment in an arbitrary unprecoloured six-colouring. The all-distance equality-to-obstruction construction is treated separately in the [companion recursion packet](recursion_escape_atlas_20261011f.md); this note owns exact two-copy counts, complete cross-edge authentication, and the field-escape controls.

The classical spindle operation is inherited, not claimed as new here. Published five-chromatic constructions require actual incompatible colour constraints, rather than field degree alone: see [de Grey's paper](https://arxiv.org/abs/1804.02385) for the 1,581-vertex construction and [Madore's field-colouring paper](https://arxiv.org/abs/1509.07023) for the arithmetic viewpoint. These citations provide context; the new controls do not reproduce those larger graphs or infer a six-equality fragment from them.

The concept board is: marked two-port support; exact rotation and source copies; additional cross edges; residue-prime field screens; and a positive colouring retained under recursion. The maps below preserve only the predicates explicitly checked.

## 2. What a faithful two-copy closure actually counts

Let $F$ be a finite graph with distinct ports $o,v$. At palette $k\ge2$, let $e$ and $d$ be the numbers of proper colourings with fixed port assignments $(0,0)$ and $(0,1)$, respectively. Form two copies sharing only $o$, and join their two copies of $v$ by a new edge. If these are exactly the graph's retained constraints, its colouring count is
\[
 C_k=k(k-1)d\bigl(2e+(k-2)d\bigr).
 \tag{1}
\]
Fix the common root colour. The two tip colours must differ. One tip can have the root colour in $2(k-1)$ ways, contributing $ed$ each; both can be non-root, different colours in $(k-1)(k-2)$ ways, contributing $d^2$ each. Multiplying by $k$ proves (1).

If the seed is $k$-colourable and $k\ge3$, the closed graph fails to be $k$-colourable **iff $d=0$**. Thus this particular recursion does not turn a positive different-port count into zero. It consumes an equality interface already present at that palette. At $k=2$ the exception is genuine: $(e,d)=(0,1)$ closes to an odd cycle with no two-colouring.

For geometry, translate $o$ to zero and suppose $|v|^2=N>1/4$. The unit rotation
\[
 \rho_N=\frac{2N-1+i\sqrt{4N-1}}{2N}
 \tag{2}
\]
satisfies $|v-\rho_Nv|=1$. It makes the proposed closing edge an actual unit edge. It does **not** guarantee that copies are otherwise disjoint or that no other cross edges exist. The implementation retains both copy maps, the intended skeleton, and the entire induced unit graph.

Equation (1) is a count for the skeleton when copies are otherwise disjoint. It is a count for the full induced graph only when there are no extra edges. Ignoring an extra edge while claiming an induced count loses an actual constraint. Conversely an explicitly retained non-induced unit subgraph is legitimate; it simply has the weaker declared edge set.

A small outside-Moser control is the path $0,1,2$ and its $\rho_4=(7+\sqrt{-15})/8$-rotated copy. Their full unit graph is $C_5$: five vertices, five edges, and $5^5-5=3120$ six-colourings. Escaping the old field has not created a hard six-colour object. The analogous three-edge path gives $C_7$; its positive two-colour different-port count exhibits the exception above.

## 3. An all-length exact geometric family

Write
\[
 a=(1+i\sqrt3)/2,\qquad v=1+a,\qquad |v|^2=3.
\]
The seed $H_m$ is a chain of $m$ triangular diamonds with vertices
\[
 p_j=jv\quad(0\le j\le m),\qquad
 x_j=jv+1,\quad y_j=jv+a\quad(0\le j<m).
 \tag{3}
\]
Each diamond has edges
\[
 p_jx_j,\ p_jy_j,\ x_jy_j,\ x_jp_{j+1},\ y_jp_{j+1}.
\]
These are **all** unit edges of the seed. For example,
\[
 |p_j-x_l|^2=3(j-l)^2-3(j-l)+1,\qquad
 |x_j-y_l|^2=3(j-l)^2+1.
\]
The other cases follow by reflection or translation; distinct $p$ vertices and distinct same-type interior vertices have squared distances at least three. Hence $H_m$ has $3m+1$ vertices and $5m$ edges.

Take $N=3m^2$, $D=12m^2-1$ and $\rho=\rho_N$. Let
\[
 G_m=U\bigl(V(H_m)\cup\rho V(H_m)\bigr),
\]
where $U$ means the **full induced** unit-distance graph.

**Complete incidence theorem.** The two copies meet only at zero, and their only new unit edge joins $p_m$ to $\rho p_m$. Consequently
\[
 |V(G_m)|=6m+1,\qquad |E(G_m)|=10m+1.
 \tag{4}
\]

Here is the cross-edge proof, including the irrationality gate. We have
\[
 (6m-1)^2<3D=36m^2-3<(6m)^2,
\]
so $\sqrt{3D}$ is irrational. For two seed vectors $p,q$, expansion of $|p-\rho q|^2$ has a rational part and a rational multiple of $\sqrt{3D}$. That latter coefficient vanishes iff $p,q$ are collinear with zero. Thus a cross-copy equality or unit distance requires collinearity.

Apart from zero, the axial points $p_j$ are one collinear family. Each off-axis ray contains just one seed point: the $x_j$ slopes are proportional to $j/(3j+2)$, strictly increasing below the axial slope, whereas the $y_j$ slopes are proportional to $(j+1)/(3j+1)$, strictly decreasing above it. Therefore the only other collinear case has $p=q=x_j$ or $p=q=y_j$.

For axial points the squared distance is
\[
 |p_j-\rho p_l|^2=3(j-l)^2+\frac{jl}{m^2}.
 \tag{5}
\]
For unequal positive indices this exceeds one; for equal indices it is one exactly at $j=l=m$. For an off-axis point $p$,
\[
 |p-\rho p|^2=\frac{|p|^2}{3m^2},\qquad
 |x_j|^2=|y_j|^2=3j^2+3j+1<3m^2 .
\]
It is strictly between zero and one. Cases with one vector zero are exactly the seed unit edges already incident to the shared root. The same equations exclude any additional coincidence. This proves (4) for every $m$, not by finite extrapolation.

The coordinates lie in
\[
 L_m=\mathbf Q(\sqrt{-3},\sqrt{-(12m^2-1)}).
\]
The four basis vectors used by the script are independent: $D\equiv2\bmod3$ is not a square, and $3D$ is not a square by the displayed inequalities. At $m=2$, the graph has thirteen vertices, twenty-one unit edges and field $\mathbf Q(\sqrt{-3},\sqrt{-47})$, outside the Moser field.

## 4. Exact criticality, count, and equality export

In every proper three-colouring of a diamond, its nonadjacent ports have the same colour: the two adjacent common neighbours use two distinct colours, leaving the third colour at each port. Both arms of $G_m$ therefore force their tips to the root colour, contradicting the closing edge.

A four-colouring is explicit. Give every axial point colour zero and the two interior vertices of each diamond colours one and two. Use the same assignment on the rotated copy, except give its last tip colour three. The incidence theorem verifies every edge. Hence
\[
 \chi(G_m)=4\qquad(m\ge1).
 \tag{6}
\]

In fact every edge deletion is three-colourable. Think of the graph as a chain of $2m$ diamonds between the two tips, closed by one edge.

- Deleting the closing edge leaves the all-equal port colouring.
- Deleting any edge of one diamond permits that diamond's two ports to have different colours. For the deleted interior edge its two common neighbours may both use the remaining colour. For a deleted spoke, one interior vertex uses a port colour and the other uses the remaining colour.
- Use equal port colours along the diamonds on either side of this single change. The closing edge now has differently coloured endpoints.

This is a constructive edge-criticality proof. Since every vertex is incident to an edge, it also implies vertex-criticality. Deleting any unit edge yields a nonvacuous **three-colour forced-equality fragment** at its endpoints: a differing-endpoint colouring would restore the edge and contradict (6). The API exports the actual deletion colouring, original coordinates, and a checked diamond-chain non-three-colourability proof. This is the positive low-palette family, not a six-colour construction.

There is also a closed count. For one diamond, conditional port counts are
\[
 e=(k-1)(k-2),\qquad d=(k-2)(k-3).
\]
Its named-colour kernel has eigenvalues
\[
 R=(k-1)(k-2)^2,\qquad S=2(k-2).
\]
For a chain of $m$ diamonds,
\[
 e_m=\frac{R^m+(k-1)S^m}{k},\qquad
 d_m=\frac{R^m-S^m}{k}.
\]
Substitution in (1) gives
\[
 P_{G_m}(k)=(k-1)(R^{2m}-S^{2m}),\qquad k\ge2.
 \tag{7}
\]
In particular the six-colour count is $5(6400^m-64^m)>0$. Every omitted actual unit edge has a separating six-colouring supplied by the original full-graph four-colouring. Neither increasing $m$ nor growing the field radicand can remove that support in this family.

## 5. Finite-prime escape can coexist with chromatic number four

The companion [field/defect frontier](unit_defect_frontier_20261011f.md) supplies a local-field screen and a general finite-prime escape rule. Combining it with (3)--(6) gives a particularly sharp hostile.

For **any finite set of rational primes** $\mathcal S$, let $P$ be the product of its odd members (empty product one), and choose
\[
 m=2P,\qquad D=48P^2-1.
\]
Then
\[
 -D\equiv1\pmod8,\qquad -D\equiv1\pmod p
       \quad\text{for every odd }p\in\mathcal S.
\]
Thus every selected prime splits in the imaginary quadratic subfield $\mathbf Q(\sqrt{-D})$. Complex conjugation exchanges its two primes. It cannot fix a prime above a selected $p$ in $L_m$: a fixed upstairs prime would restrict to a fixed downstairs prime.

Accordingly the particular **conjugation-stable-prime residue-colouring mechanism** of THM-4558 Lemma R is unavailable over all selected rational primes. Yet $G_m$ is still explicitly four-colourable. This does not eliminate other colourings, other local methods, or screens at new primes. Its quantifier is
\[
 \forall\ \text{finite }\mathcal S\quad
 \exists\ \text{a specified four-chromatic geometric }G_m
 \text{ whose field escapes those screens}.
\]
It is not one field escaping every prime.

A separate input from [the equation recursion packet](divisor_equation_recursions_20261011f.md) has rotation radicand
\[
 D_*=3\cdot1237\cdot10193\cdot11329\cdot165713
     =71013535489456671.
\]
It lies outside the Moser field and $D_*\equiv7\bmod8$, so it escapes the prior prime-two obstruction. But in the joined field
$\mathbf Q(\sqrt{-3},\sqrt{-11},\sqrt{-D_*})$, prime eleven still supplies an inherited five-colour screen. Indeed $-3$ and $-D_*$ are nonsquares modulo eleven, whereas $\sqrt{-11}$ is ramified. The local decomposition group contains the common Frobenius flip of the two unramified generators and the inertia flip of $\sqrt{-11}$; their product is complex conjugation. Its residue action is Frobenius on $\mathbf F_{121}$. Lemma R and the inherited value $\kappa(11)=5$ apply. This last chromatic bound is **inherited**, with the exact Legendre markers checked here; it is not inferred from an unsuccessful finite search.

## 6. The lost-cross-edge hostile and the actual search frontier

Take the complete radius-two triangular shell (nineteen points), retain a tip of squared distance three, and rotate the shell by $\rho_3$. The resulting point set has thirty-seven distinct vertices and ninety actual unit edges. The intended two-copy skeleton with one selected closing edge has only eighty-five edges: five more unit edges join other corresponding norm-three tips.

The script gives a four-colouring of that skeleton which violates an additional actual edge. Thus the seed's two conditional counts and one marked closing edge do not determine the full induced relation. The missing information is the complete cross-edge incidence, not just a larger numerical count.

The useful recursive frontier is therefore precise: an operation must retain exact copy geometry and obtain new joint interior constraints, or it must consume a genuinely stronger equality seed. Field escape and private-copy count multiplication alone do not provide those constraints. A new six-equality claim still owes both an actual equal-colour witness and a complete exclusion of different-port colourings on the same geometric source.

## 7. Program contract and reproduction

The exact [script](../../04-computation/experiments/unit_field_escape_recursion_20261011f.py) uses rational coefficients in the independent basis $1,s,t,st$, where $s=i\sqrt3$ and $t=i\sqrt D$ (so $st=-\sqrt{3D}$ in the chosen embedding). **Spindle** retains the source **Seed**, both exact copy maps, the closing ports, intended skeleton and full unit-edge set. **audit_spindle** authenticates all of them; **faithful** requires only the pivot to be shared and no extra edges before an induced count is used. Geometry is validated before cached lookup.

**chain_three_refutation** checks actual diamond implications and their common root; **chain_four_colouring** and **deleted_three_colour** export checked literal colourings. **separated_omission** retains the original graph and omitted actual unit pair. Malformed coefficient types, Boolean aliases, incomplete incidence and lost source maps are rejected. No bound on running time for arbitrary large user inputs is asserted; the selected controls are bounded.

The declared finite universe is: all fixed-assignment count pairs $0\le e,d\le4$ at palettes two through eight; exact five- and seven-cycle controls; diamond chains $1\le m\le12$ with every pairwise distance checked and all 792 unit-edge deletions/separations; literal small seed counts; the thirty-seven-point extra-edge hostile; six specified finite-prime sets; and malformed receipts. Large graphs in the prime-escape table are represented by the proved formula rather than materialized or falsely described as enumerated.

~~~
python -B -X utf8 04-computation/experiments/unit_field_escape_recursion_20261011f.py
python -B -O -X utf8 04-computation/experiments/unit_field_escape_recursion_20261011f.py
~~~

Normal and optimized execution give **2,239 always-active checks** and agree with the [saved output](unit_field_escape_recursion_20261011f.out). Its LF-normalized SHA256 is **bcab240c1a663032abbfd6488e6826cafbb98ee88aecb464738f02c3c95d989e**. The construction and its stopping reason are exact; a six-colour equality gadget and the plane's six-versus-seven question remain open here.
