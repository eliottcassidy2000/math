# Tournament cell spaces: free line memory and profile-invisible surfaces

**Status: all-order constructions and proofs; finite exact controls.** The blue/black complement-line multigraph and the single-arc-flip class graph are different objects. This note constructs explicit cell spaces from each, rather than assigning an unstated topology to either graph. The main results are:

1. Identifying isomorphic corners of the marked chord cube creates a free fundamental group of rank `2^m − k`, where `m=(n−1)(n−2)/2` and `k` is the number of tournament classes. An incidence-preserving thickening contains the actual complement-line graph and its component fundamental groups as free factors. Thus filling chord squares does not, by itself, kill any reduced complement-line cycle.
2. One inherited equal-profile tournament edge compiles every real moment-angle cubical complex. In particular, it supplies explicit equal-profile spheres, tori, and orientable surfaces of unbounded genus, with every edge a literal designated arc reversal.
3. The same tournament vertices and the same graph edges can support radically different homotopy types when the selected higher cells differ. Complementing all marked chord bits, identifying tournament classes, and taking a homotopy quotient are also different operations.

The topology of real moment-angle complexes and the elementary graph/cofibration facts are standard and credited below. The contribution here is the explicit composition with the repository's physical tournament interfaces and the relative free-factor theorem for its complement matching. No literature-priority claim is made.

## 1. Three levels of state must be kept distinct

For `n≥3`, keep the directed path `0→1→⋯→n−1`. Its nonconsecutive pairs give `m=(n−1)(n−2)/2` chord bits. Let

- `A={0,1}^m` be the marked flags;
- `c:A→V` send a flag to its tournament isomorphism class;
- `F(t)=1−t` complement the chord bits;
- `S` reflect chord positions along the path.

Every class occurs because every tournament has a directed Hamilton path. Each pair `{t,Ft}` is one individual line in the original S151 multigraph `Γ`. Thus `|A|=2E`, where `E=2^(m−1)`, and `|V|=k`. A loop retains its two distinct flags. Blue means `St=t`; otherwise a line is black. We retain all line identities and all class vertices, without merging converse classes.

Two immediate cautions are essential:

- `F` is not tournament converse, and generally does not define a function on `V`.
- An interior point of the cube `I^m` is a formal cubical interpolation. It is not a tournament with a new class label. The cells encode allowed commuting chord edits.

The later product construction instead uses the single-arc-flip class graph `M_N`. It is not `Γ`.

## 2. The corner-identification cube has a completely explicit homotopy type

Define

`X_n = I^m / (t∼u exactly when t,u are corners and c(t)=c(u)).`

Only corners are identified. No positive-dimensional cube face is identified with another, and no class quotient is applied to the cube's interior.

**Theorem 1 (corner quotient).** `X_n` is homotopy equivalent to a wedge of `2^m−k` circles. Consequently

- `π₁(X_n)` is free of rank `2^m−k`;
- `H₁(X_n;Z)=Z^(2^m−k)`;
- `H_j(X_n;Z)=0` for `j≥2`;
- its universal cover is contractible.

**Proof.** The corner inclusion `A→I^m` is a CW cofibration. Therefore the ordinary pushout defining `X_n` has the homotopy type of the double mapping cylinder

`Y = I^m ∪ (A×[0,1]) ∪ V`,

where `(t,0)` attaches at the cube corner `t` and `(t,1)` attaches at `c(t)`. Collapse the contractible cube subcomplex of `Y` to a hub `h`. This quotient is a homotopy equivalence and produces a graph `G`: its vertices are `h` and the `k` classes, and it has one distinct spoke `a_t:h→c(t)` for every marked flag. It is connected, has `2^m` edges and `k+1` vertices, hence has cycle rank `2^m−k`. A connected graph collapses a spanning tree to give the stated bouquet. ∎

This is an identification-space theorem, not a statement that two arbitrary tournament classes have a preferred continuous deformation between them. The retained corner multiplicities carry the extra loops.

## 3. Complement-line cycles survive as free factors

Every complement pair gives a radial path `c(t)←h→c(Ft)` in `G`. The corresponding literal cube diagonals all meet at the center. They are a map of graph paths, **not an embedded copy of the complement-line graph**.

To keep the actual lines distinct, define a two-dimensional CW complex `H_n` as follows:

1. Start with the flag-spoke graph `G` of Theorem 1.
2. Adjoin every individual edge `e={t,Ft}` of `Γ`, with the original class endpoints and its original color.
3. Attach a triangular two-cell with boundary consisting of that new line and the two flag spokes `a_t,a_Ft`.

There is one triangle per complement pair. For a graph loop the triangle has two identified class corners, but its three open edges remain distinct. This is an ordinary CW attaching map. The full multigraph `Γ` is now an actual subcomplex of `H_n`.

**Theorem 2 (relative free-factor theorem).** Let `c_Γ` be the number of connected components of `Γ`. Then

`π₁(H_n) ≅ (* over components Γ_i of π₁(Γ_i)) * F_(E−c_Γ)`,

with suitable choices of paths from the common base point to those components. Each component inclusion is injective on fundamental groups. The total free rank is `2E−k=2^m−k`. Moreover `H_n≃X_n`.

**Proof.** Each triangle's newly added complement line is a free open edge in the two-dimensional complex: it belongs to exactly that triangle. Collapsing each triangle across its complement edge recovers `G`, proving `H_n≃G≃X_n`.

For the relative statement, subdivide each new line at its own midpoint and divide its triangle into two triangles by adding a hub-to-midpoint spoke. Each original flag spoke `a_t` is now a free edge of exactly one half-triangle. Collapse all those half-triangles across their flag spokes. The result is the subdivided `Γ` together with the hub and exactly one hub-to-midpoint spoke for every line. No line edge was removed. This collapse is also compatible with path reflection `S`, since its choices are identical at every flag. Every component of `Γ` contains a line, since every class has at least one marked flag and all flags are paired. Choose one spoke from each component. These `c_Γ` spokes join the components to the hub without creating a cycle. Each of the remaining `E−c_Γ` spokes adds one free generator. This proves the free-product description, including the inclusion statement. Finally

`Σ_i rank π₁(Γ_i) + E−c_Γ = (E−k+c_Γ)+(E−c_Γ)=2E−k.` ∎

An equivalent construction starts with the literal corner quotient `X_n`, adjoins separate complement edges, and fills each by a disk whose other boundary is its radial cube path. It is a homotopy-equivalent extension. One must still keep the separate edges: intersecting cube diagonals are not a substitute for the incidence graph.

**Corollary 3 (no newly killed line words).** A closed line walk in one component of `Γ` becomes null-homotopic in `H_n` if and only if it already reduces to the empty word by cancelling consecutive reverse traversals in the graph. The chord-cube two-cells introduce no new relation among the physical complement-line cycles.

This is especially useful beside S166 and S169. A torus cell imposes a commuting relation on transport. The corner-identification cube, with the specified matching thickening, does not impose a commuting relation between arbitrary complement-line cycles.

### 3.1 Every equal-degree line local system extends

A `d`-sheet graph local system assigns a `d`-element fiber to each vertex and a bijection to each oriented line, with the inverse on reverse traversal. It may have nonabelian monodromy. Assume the same degree `d` on each component of `Γ`.

**Corollary 4 (extension without a commutativity condition).** Every such local system on `Γ` extends across `H_n`. Its restriction is the prescribed local system, not just a conjugate collection of cycle signs.

**Constructive proof.** Give the hub a `d`-element fiber. For every oriented complement pair `t→Ft`, with line transport `L_e`, choose any bijection `A_t` from the hub fiber to the fiber at `c(t)`, and set

`A_Ft = L_e ∘ A_t`.

Every flag occurs in exactly one pair, so these assignments never conflict. The triangular boundary transport is the identity by construction. The cover extends over the triangles. No two unrelated line transports have been required to commute. ∎

The added hub and flag transports are auxiliary data. This is not an assertion that new sheets are tournament classes. The physical S166 binary cover is defined only on its NS-only black quotient `Q₁`; nothing here silently extends that local system through an SC branch vertex. The extension theorem can instead be applied to any specified graph with the matching incidence structure.

### 3.2 Exact finite controls

The checker reconstructs every tournament class and every fixed-path flag at `n=3,4,5,6`, retaining loops and individual lines. The actual complement graphs happen to be connected in this range:

| n | flags | classes k | lines E | rank π₁(Γ) | rank π₁(H_n) | added free rank |
|---|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 1 | 0 | 0 | 0 |
| 4 | 8 | 4 | 4 | 1 | 4 | 3 |
| 5 | 64 | 12 | 32 | 21 | 52 | 31 |
| 6 | 1,024 | 56 | 512 | 457 | 968 | 511 |

No all-order connectedness claim is inferred. The theorem includes arbitrary `c_Γ`. For `n≤5`, the checker also constructs the entire identified-corner cubical chain complex and verifies its mod-two homology independently of the graph rank formula.

### 3.3 The physical converse reflection has an actual fixed-locus topology

Unlike complement `F`, path reflection `S` does descend to `X_n`, because `c(St)=c(t)^op`. Put

`d=floor((n−1)²/4)`, and let `k_SC` be the number of self-converse classes.

The fixed cube `(I^m)^S` has one coordinate for every reflection orbit of chord positions, so it is a `d`-cube. Its `2^d` corners are precisely the reflection-fixed flags. S151 proves that every SC class has an odd number of fixed flags, so every SC class occurs among these corners.

**Theorem 4a (fixed locus and reflection character).** The fixed-point space `X_n^S` is the corner quotient of this `d`-cube by the SC class map. In particular,

`X_n^S ≃ wedge of (2^d−k_SC) circles`.

On rational first homology of the full corner quotient,

`Tr(S | H₁(X_n;Q)) = 2^d−k_SC`,

and its invariant and anti-invariant dimensions are

`b₁^± = (2^m−k ± (2^d−k_SC))/2`.

**Proof.** A point not coming from a cube corner is fixed after quotient precisely when it is already fixed before quotient, since no such points are identified. A corner image is fixed precisely when its class is SC. Each SC image has a fixed flag representative by the inherited oddness statement. This identifies the fixed-point set with the stated corner quotient; the continuous bijection is a homeomorphism because the spaces are compact and Hausdorff. Theorem 1 applies with `2^d` corners and `k_SC` classes.

For the homology character, use the equivariant flag-hub graph `G`. Its one-chain representation is the permutation representation on all flags, while its zero-chain representation is the class permutation representation plus the fixed hub. Since `G` is connected, `H₀` contributes one fixed dimension. Thus the alternating character identity gives `Tr(S|H₁)=2^d−(k_SC+1)+1`. Combining trace with total dimension `2^m−k` gives the two eigenspace dimensions. ∎

The fixed-point homology rank and the trace happen to agree in this construction; they are not identified as the same vector space. The equivariant homotopy model, rather than a generic fixed-point claim, proves the equality.

The relative collapse also identifies the extra rational homology representation precisely. If `f_Γ` is the number of components of `Γ` preserved by `S`, its extra summand has dimension `E−c_Γ` and trace

`2^(d−1)−f_Γ`.

Indeed its character is the permutation character on individual lines minus the permutation character on components; only blue lines are fixed, with no edge inversion. Thus

`Tr(S|H₁(Γ;Q)) = 2^(d−1)−k_SC+f_Γ`,

and the traces add to the one in Theorem 4a. In the verified order-six example, the full corner quotient has `(b₁^+,b₁^−)=(510,458)`, while the actual line graph has `(239,218)`.

## 4. One physical edge compiles cubical topology at constant path profile

Use S166's strong six-vertex tournament masks `18` and `22`, with lexicographic pair bits and bit `1` meaning `i→j` for `i<j`. Call them `A` and `B`. They are nonisomorphic, differ in one designated arc, and have the common path-cover polynomial

`P(z)=15z+93z²+138z³+73z⁴+15z⁵+z⁶`.

Equivalently, the shared-slot polynomial is

`F(r)=r⁶+8r⁴+6r²`,

and both have `H=15` directed Hamilton paths. These seed facts are inherited; the checker independently verifies them by full relabelling and exact path-cover enumeration.

For a bit vector `x∈{0,1}^k`, set

`T_x = A_(x₁) ⊕ ⋯ ⊕ A_(x_k)`, where `A_0=A` and `A_1=B`.

S166's ordered-strong-component theorem gives an induced cube `Q_k` in the single-arc-flip class graph `M_(6k)`. The SCC positions intrinsically retain the bit coordinates, so forgetting vertex labels does not identify these `2^k` tournaments. Every bit edge is the designated literal arc reversal inside its block. Every vertex has

`F_(T_x)(r)=(r⁶+8r⁴+6r²)^k`, and `H(T_x)=15^k`.

A collection of simultaneous coordinate flips is represented by its cube face. Its vertices and edges are physical tournament data; its higher-dimensional interior is the chosen commuting-move cell, not a continuum of additional tournaments.

## 5. The tournament real-moment-angle compiler

Let `K` be any finite abstract simplicial complex on the coordinate set `[k]`, containing every singleton. Define

`Z_K = ⋃_(σ∈K) I^σ × {0,1}^([k]\σ) ⊆ I^k`.

A cell is specified by a free-coordinate set `σ∈K` and fixed bits in all other coordinates. This is the standard **real moment-angle complex**, also called the polyhedral product `(I,S⁰)^K`.

**Theorem 5 (constant-profile topology compiler).** Label every vertex `x` of `Z_K` by the tournament class `[T_x]` above. Then:

1. All `2^k` vertices are distinct tournament classes on `6k` vertices and have the same complete path-cover profile.
2. Its one-skeleton is the induced physical cube `Q_k`, independently of `K`.
3. Every cell uses commuting designated flips in distinct SCC blocks.
4. The resulting cell space is exactly `Z_K`, so its entire homotopy type is transferred, not merely its Betti numbers.
5. If strong tournament vertices are required, S166's wrapper transfers the family to order `30k+10`, preserving its induced graph and common full profile. The wrapped Hamilton-path count is not asserted to remain `15^k`.

**Proof.** The ordered-component construction establishes the first three statements and makes the assignment injective on the cube's vertex labels. We retain precisely the cells in the displayed definition, with their original cubical face maps, which proves the fourth. The last statement is the inherited induced strong wrapper together with its module-substitution profile lemma. ∎

The selected `K` is essential extra information: it specifies which simultaneous move sets count as filled cells. All module flips commute physically, so selecting fewer cells is a declared restriction, not a theorem that the omitted moves fail to commute. The unmarked graph alone cannot recover `K`.

## 6. Spheres, tori, and an explicit unbounded-genus family

### 6.1 Spheres

Take `K=∂Δ^(k−1)`, all proper coordinate subsets. Then `Z_K=∂I^k≅S^(k−1)`. Thus for every `k≥2`, a sphere of dimension `k−1` has a cubical realization using `2^k` equal-profile tournaments of order `6k`, with physical one-arc flips as edges.

### 6.2 Tori at every dimension

Partition `2d` coordinates into `d` pairs. Allow a free set precisely when it contains at most one coordinate of each pair. Thus `K` is the join of `d` two-point complexes. Directly from the product definition,

`Z_K = (∂I²)^d ≅ (S¹)^d`.

This torus has

- `4^d` tournament vertices of order `12d`;
- common Hamilton-path count `15^(2d)`;
- `4^d binom(d,j)` cubical cells in dimension `j`;
- integral Betti number `binom(d,j)` in degree `j`.

The optional strong wrapper has tournament order `60d+10`. No arithmetic winding or event schedule is required to specify this torus: it is a finite cubical state space. A chosen motion on it would still need its own rule.

### 6.3 Surfaces of genera 0, 1, 5, 17, 49, ...

Take `K=C_k`, the cycle on coordinates `0,...,k−1`, with its vertices and edges and no higher faces; require `k≥3`. The allowed two-cells vary consecutive coordinate pairs cyclically.

**Theorem 6 (surface family).** `Z_(C_k)` is a connected closed orientable surface of genus

`g_k = 1 + (k−4)2^(k−3)`.

It uses `2^k` tournament vertices of order `6k`, all with full profile `(r⁶+8r⁴+6r²)^k` and `H=15^k`.

**Proof.** Every vertex link is the cycle `C_k`, so the square complex is a surface. Each coordinate edge belongs to exactly two squares. Its one-skeleton is `Q_k`, hence connected. The cell counts are

`V=2^k`, `E=k2^(k−1)`, `F=k2^(k−2)`.

For an explicit orientation, orient a face on cyclically consecutive coordinates `(i,i+1)` in that cyclic order and multiply its orientation by `(-1)^(sum of all fixed bits)`. At an edge varying coordinate `i`, the two incident faces are `(i−1,i)` and `(i,i+1)`. Their induced coefficients are opposite: both depend on the sum of all bits fixed along the edge, with one extra minus sign. Thus the sum of these oriented faces is an integral fundamental cycle with coefficient `±1` on each face. The surface is orientable.

Finally `χ=V−E+F=2^(k−2)(4−k)`. The equality `χ=2−2g` gives the formula. ∎

The first examples are:

| k | surface | tournament vertices | tournament order | common H |
|---|---|---:|---:|---:|
| 3 | sphere | 8 | 18 | 3,375 |
| 4 | torus | 16 | 24 | 50,625 |
| 5 | genus 5 | 32 | 30 | 759,375 |
| 6 | genus 17 | 64 | 36 | 11,390,625 |
| 7 | genus 49 | 128 | 42 | 170,859,375 |

These are unbounded genera, not every genus. The `k≥4` surfaces are aspherical. Their many nontrivial homotopy classes are invisible to every path-cover coefficient at a vertex, while still visible to the selected cubical topology.

## 7. Right-angled groups: the broad adjacent connection

For a finite flag complex `K`, let

`W_K = <s₁,...,s_k | s_i²=1, s_i s_j=s_j s_i when {i,j}∈K>`.

A standard theorem identifies `π₁(Z_K)` with the commutator subgroup `[W_K,W_K]`, the kernel of `W_K→(Z/2)^k`, and says that `Z_K` is aspherical for flag `K`. See Grbić–Ilyasova–Panov–Simmons, cited below. These are external topological inputs, not new tournament theorems.

Theorem 5 therefore realizes these standard cubical classifying spaces using a uniform equal-profile tournament hardware family. For `K=C_k`, `k≥4`, this is a surface group; for a join of `d` vertex pairs it is `Z^d`. It bridges physical commuting arc flips with geometric group theory without replacing noncommuting holonomy by a false torus model.

For an arbitrary permutation transport on the one-skeleton, extension over `Z_K` is equivalent to identity transport around every selected square. Higher cubes impose no additional obstruction to extending a covering once all their square relations hold. Omitted squares impose no relation. This is a precise choice of relations, unlike merely asserting that a graph cycle “has a homotopy.”

## 8. One graph, several incompatible homotopy types

Fix `k=4`. The same sixteen tournament classes on 24 vertices have the same 32 induced flip edges and the same full profile, with `H=50,625`. Depending on the declared higher cells, these data support:

- no two-cells: the graph `Q₄`, with free fundamental group of rank 17;
- the 16 squares specified by `C₄`: a two-torus;
- all proper cube faces: the 3-sphere `∂I⁴`;
- the full four-cube: a contractible space.

This is not an ambiguity in the definition of any one space. It shows why the higher-cell specification is indispensable. A complete path-cover profile cannot distinguish the vertex data, and the abstract one-skeleton cannot distinguish these fillings.

## 9. Complement quotients: three more different spaces

Keep the **marked** chord cube and let `F(x)=1−x`. Radial projection from the center identifies its boundary equivariantly with the antipodal sphere `S^(m−1)`. Therefore

`∂I^m / <F> ≅ RP^(m−1)`.

The full cube is the cone on its boundary with the center as apex. Hence

`I^m / <F> ≅ Cone(RP^(m−1))`,

which is contractible. Its homotopy quotient is different again:

`E(C₂) ×_(C₂) I^m ≃ B(C₂)=RP^∞`,

because radial contraction to the center is equivariant. This last construction retains stabilizer information. The center is fixed by `F`, even though `F` acts freely on the binary corners.

For `m≥3`, the boundary quotient has fundamental group `C₂`; for `m=2` it is a circle and its antipodal covering has total-space fundamental group of index two in `Z`; for `m=1` the boundary quotient is a point. Actual fixed-path tournament chord dimensions in the present `n≥3` family start at `m=1,3,6,10,...`.

None of these spaces is the corner-identification quotient `X_n`. In particular, `F` usually does not descend through `c`, so one may not place an unproved antipodal action on `X_n` and identify its quotient with a projective space.

Likewise, the full weighted-tournament cube on all `binom(n,2)` arc coordinates admits relabelling by `Sym(n)` through signed coordinate permutations. The all-one-half point is fixed, and radial contraction is equivariant. Its coarse orbit space is contractible, whereas its homotopy quotient has type `B Sym(n)`. That statement concerns the full weighted cube, including fractional coordinates; it is not the simple class flip graph and does not erase tournament stabilizer information by fiat.

## 10. Reproduction and scope

Run

```
python3 check.py > observed.json
python3 -O check.py > observed_optimized.json
cmp observed.json results.json
cmp observed_optimized.json results.json
```

The checker uses only the Python standard library and explicit failures, not `assert`. It independently verifies:

- the six-vertex seed pair's strong connectivity, nonisomorphism, one-flip relation and full path-cover polynomial;
- all original complement flags, class counts and graph components through order six;
- the thickening chain ranks, including all individual loops and line identities, and the exact reflection/fixed-class counts used in the homology character formulas;
- the complete corner-quotient cubical homology through order five;
- integer `∂²=0`, integral surface orientations, every edge's two-face incidence, cell counts and mod-two Betti numbers for `3≤k≤9`;
- torus cell and Betti counts for `1≤d≤4`, and sphere controls for `2≤k≤7`.

The integer surface orientation and manifold proof, rather than mod-two homology alone, establish orientability. The all-order proofs do not depend on extrapolating the finite checks. The compiler does not imply that every finite complex has its homotopy type among these polyhedral products, that every abstract voltage is physically realized by an original tournament metagraph, or that a static cubical space supplies an autonomous dynamical law.

## Sources and provenance

- [S151 original and merged line structure](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session151_spectrum_metagraph_20261103/metagraph/STRUCTURE.md): fixed-path flags, `F`, `S`, distinct line identities, and the warning that complement does not descend to tournament classes.
- [S166 reflection voltages](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md): the specified NS-only binary cover, its quotient markings and the failure of a local system at an SC branch vertex.
- [S166 Cartesian flip constructions](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/cartesian_flip/NOTE.md): the equal-profile masks `18,22`, induced product theorem, strong wrapper and module-substitution profile preservation. These physical ingredients are inherited.
- [S169 torus transport](https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session169_runner_jacobian_20261103/torus/NOTE.md): the exact distinction between graph monodromy and a cover extending across commuting cells.
- Allen Hatcher, [*Algebraic Topology*, Chapter 0](https://pi.math.cornell.edu/~hatcher/AT/ATch0.pdf) and [Chapter 1](https://pi.math.cornell.edu/~hatcher/AT/ATch1.pdf): standard cofibration, quotient, graph fundamental-group and covering constructions.
- Jelena Grbić, Marina Ilyasova, Taras Panov and George Simmons, [*One-relator groups and algebras related to polyhedral products*](https://arxiv.org/abs/2002.11476), published version 2022: standard real moment-angle, right-angled Coxeter and asphericity background. Section 7 uses these established results; the explicit surface proof above is included to make that family self-contained.
