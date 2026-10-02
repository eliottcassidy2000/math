# Is the product of two Petersen graphs 4-polytopal? An attempt: structure, a snark lemma, intrinsic linking, exact SAT exclusions, and Petersen × triangle is not polytopal

**opus S15, 2026-10-01.**

- Script: `04-computation/experiments/petersen_product_4polytope_20261001.py` (+ `.out`).
- Previous note: [`odd_zeta_parallels_petersen_20261001.md`](odd_zeta_parallels_petersen_20261001.md) §1 (status of
  Ziegler's question; the "≤ 8 mixed 2-faces per vertex" remark).

**Status.**
- **Ziegler's question in dimension 4: still OPEN.** No proof that `P□P` is not 4-polytopal is claimed.
- PROVED (elementary): the structure lemmas of §2 and the snark lemma of §3.
- PROVED (from the Conway–Gordon–Sachs theorem): the linking lemma of §4.
- COMPUTED (exact SAT; necessary conditions only, with sound lazy cuts):
  - **`P×C₃` (Petersen × triangle) is not polytopal** (§6.2). The run uses all 7681 induced cycles as candidate
    2-faces. PPS give this graph as an example that local data cannot decide: it is the graph of a cellular
    `S¹×RP²`.
  - Validation: the same pipeline finds PPS's cellular `S¹×RP²` with graph `K₃,₃×K₃` when the homology condition
    is dropped, and rejects it when the condition is kept, as PPS announce.
  - Two sub-cases of `P□P` are impossible (§6.3).
- COMPUTED (negative): the purely local 2-face conditions are satisfiable, and so are the full local conditions
  imposed on 9 of the 10 fibres. The obstruction, if any, is global (§6.4).
- Independent audit: PENDING (running).

## 0. The owner's prompt

"Try proving the product of two Petersen graphs is not 4-polytopal."

## 1. Setting

- `P` is the Petersen graph. `G = P□P` has 100 vertices `(u,v)` and 300 edges, and is 6-regular.
- An edge is **horizontal (H)** if it changes `u` and **vertical (V)** if it changes `v`. The fibres are `P×{v}` (H) and
  `{u}×P` (V).
- Suppose `G` is the graph of a 4-polytope `Π`. Standard facts used throughout:
  - Every face's graph is an induced subgraph. So 2-faces are induced cycles, and facets are induced 3-polytopal
    subgraphs.
  - Two faces meet in a common face. Two 2-faces share nothing, a vertex, or an edge.
  - The vertex figure at `x` is a 3-polytope on the 6 neighbours of `x`. Its graph `Γ_x` has an edge `ab` iff the
    angle `a–x–b` lies in a 2-face; each angle lies in at most one 2-face. Its faces are the facets at `x`.
  - The edge figure at `xy` is a polygon. The 2-faces through `xy` are cyclically ordered, the same way (up to
    reversal) as seen from `x` (the rotation at `y` in `Γ_x`) and from `y`.
  - `∂Π` is a 3-sphere, so `H₁(∂Π; Z/2) = 0`, and the 2-faces' boundaries span the cycle space of `G`
    (dimension 201).
- Pfeifle–Pilaud–Santos (PPS, Israel J. Math. 2012) excluded dimension 6 and left 4 and 5 open. Their remark that
  local data cannot decide refers to the 4-manifold `RP²×RP²`, i.e. to dimension 5. **In dimension 4 the vertex
  links must be 2-spheres on 6 points; that is a strong local constraint, and this note tests it.**

## 2. Elementary structure of a hypothetical realization

- **L1 (no facet lies in a fibre).** A proper induced subgraph of `P` has a vertex of degree ≤ 2, and `P` itself is
  non-planar.
- **L2 (both directions turn somewhere at every vertex).** `Γ_x` has at least one pure H edge (two H-neighbours),
  at least one pure V edge, and at most 8 of the 9 mixed edges.
  - Without a pure V edge, each V-neighbour needs all three H-neighbours, so `K₃,₃ ⊆ Γ_x`.
  - The census of the 1227 labelled polyhedral graphs on `{H₁,H₂,H₃,V₁,V₂,V₃}` gives (mixed, pure H, pure V) from
    (3,3,3) to (8,3,1). Pure H ≥ 1, pure V ≥ 1 and mixed ≤ 8 always hold.
- **L3 (pure 2-faces are long).** A 2-face through a pure angle has ≥ 5 vertices, because two vertices at distance 2
  in `P` have a unique common neighbour.
- **L4 (squares and their corners).** Let `s = e×f` be a square (4-cycle) with opposite corners `x, x'`.
  - If the angle of `s` at `x` lies in a 2-face `Q ≠ s`, then the angle of `s` at `x'` lies in no 2-face.
  - Reason: a 2-face `Q'` there would share two non-adjacent vertices with `Q`; and `Q = Q'` would force `Q = s`.
  - So a non-face square has at most two realized corner angles, at adjacent corners, both on bent faces. A face is
    **bent** if it is neither a square nor inside a fibre.
- **L5 (2-faces inside one fibre).** `P`'s induced cycles are its 12 pentagons and 10 hexagons. Families of
  pairwise properly meeting ones are contained in one of four maximal families (computed):
  - the 6 face pentagons of either hemi-dodecahedral embedding `E₁`, `E₂` of `P` in `RP²`;
  - a complementary (vertex-disjoint) pair of pentagons;
  - a hexagon `C_a = P − N[a]` with the 3 pentagons through `a` from one embedding.
  - In particular, at most one hexagon per fibre. Each complementary pair has one pentagon in `E₁` and one in `E₂`.

## 3. A snark lemma: the prism scenario is the 3-edge-colouring problem

**Lemma 1.** Suppose every vertex figure is a triangular prism whose triangles are `{H₁,H₂,H₃}` and `{V₁,V₂,V₃}`,
and every mixed angle lies in a square 2-face. Then the Petersen graph is 3-edge-colourable.

*Proof.*
- At `x = (u,v)` the three rungs of the prism pair the edges at `u` bijectively with the edges at `v`.
- For an edge `f` of the second factor, let `M_f = {e : e×f is a 2-face}`. At every `x = (u,v)` with `v ∈ f`,
  exactly one edge `e ∋ u` lies in `M_f`, so `M_f` is a perfect matching.
- For the three edges `f₁, f₂, f₃` at `v`, the matchings are disjoint (each `e ∋ u` has exactly one rung), so
  `M_{f₁} ∪ M_{f₂} ∪ M_{f₃} = E(P)`. That is a 1-factorization. ∎

So the most symmetric local picture — every fibre a hemi-dodecahedron of pentagon 2-faces, glued by "colour-matched"
squares — is exactly what Petersen's being a snark forbids. For a 3-edge-colourable cubic graph `H` that tiles a
surface, the same recipe gives a candidate cell structure for `H□H`.

## 4. Intrinsic linking

**Lemma 2 (Conway–Gordon–Sachs in a fibre).** In a 4-polytope with graph `P□P`, every fibre `P×{v}` has a
complementary pentagon pair `(C, C')` with `lk(C×{v}, C'×{v})` odd. Hence neither `C×{v}` nor `C'×{v}` is a 2-face,
and **no fibre has all six pentagons of `E₁` (or of `E₂`) as 2-faces.**

*Proof.*
- `∂Π ≅ S³` contains the fibre as a subgraph of its 1-skeleton. By Sachs and Conway–Gordon, in every embedding of
  `P` in `S³` the linking numbers over all pairs of disjoint cycles sum to an odd number. The only disjoint cycle
  pairs of `P` are the 6 complementary pentagon pairs.
- A 2-face is an embedded disc meeting the 1-skeleton only in its boundary. So its boundary has linking number 0
  with every disjoint cycle of the 1-skeleton.
- If all six `E₁`-pentagons were 2-faces, every pair would contain one, and the sum would be 0. ∎

More generally, a pair is unlinked as soon as one member bounds a 2-chain of 2-faces avoiding the other. An example
is an annulus of squares plus a cap one level up; `cgs`-clauses of this kind are valid extra constraints.

**Observation (computed).** Let `β` be a symmetric bilinear form on `Z₁(P; Z/2)`. Then the sum of `β(C, C')` over the
6 complementary pairs is 0.
- `β_{v,w}(Z, Z') = lk(Z×{v}, Z'×{w})` is a bilinear form for distinct fibres.
- Suppose all 15 squares `e×vw` are 2-faces (a "fully faced ladder"). Then each `C×{w}` can be slid to `C×{v}`
  without meeting `C'×{w}`, so Lemma 2 forces the sum `Σ β_{v,w}(C, C')` to be odd.
- Hence **the linking form across a fully faced ladder is never symmetric.** The CGS parity lives entirely in the
  antisymmetric part. (Heuristically, that part records how the three squares are ordered around each vertical
  edge; this is not proved here.)

## 5. The SAT model (all constraints are necessary conditions)

**Variables.**
- `F_Q`: candidate 2-face `Q` (an induced cycle of `G`) is chosen.
- `A(x;a,b)`: the angle `a–x–b` is covered by a chosen face.
- `M(x,y;z,w)`: the 3-path `z–x–y–w` lies in a chosen face.
- `C(x,y;z,z')`: `z, z'` are consecutive around `y` in `Γ_x`.

**Static constraints.**
- **(I) Proper intersection.** At most one chosen face contains any given *non-adjacent* vertex pair. This exactly
  captures improper meeting of two induced cycles.
- **(V) Vertex figures.** `Γ_x` is 3-connected and planar. The CNF (separations of ≤ 2 vertices; no `K₅`, subdivided
  `K₅` or `K₃,₃` on 6 vertices) was checked against brute force: 1, 25, 1227 labelled graphs on 4, 5, 6 vertices.
- **(R) Edge figures.** `C` is defined from the `A`'s by Tutte's theorem: in a polyhedral graph the faces are the
  induced non-separating cycles. The resulting DNF was checked against planar embeddings on all 1227 + 25 + 1 graphs.
  Two faces through `xy` must be consecutive around `y` in `Γ_x` iff they are consecutive around `x` in `Γ_y`.

**Lazy constraints (exact cuts).**
- **(H) Homology.** If the chosen faces do not span `Z₁(G; Z/2)`, find a cochain `w` that vanishes on them but not on
  `Z₁`. Add the clause "some face with odd `w`-weight is chosen". It is valid because a sphere's 2-faces span `Z₁`.
- **(F) Facets.**
  - A *germ* is a face of `Γ_x`. Germs at `x` and `y` glue when they contain the same two consecutive 2-faces through
    `xy`. A facet is a glued class.
  - Each class must be a 3-polytope: one germ per vertex (F1), induced (F4), Euler characteristic 2 (F3),
    3-connected (F5). Two facets must meet in a common face (F6).
  - A germ is encoded exactly by "induced and non-separating in `Γ_x`".
- **Soundness of the facet cuts.**
  - **F1, F4.** The chosen faces along a path of glued germs force both ends into one facet, so a forbidden
    repetition (F1) or a missing edge `ab` with `a ~ b` in `G` (F4) excludes the path's literals.
  - **F6.** Paths in both facets from `a` to `b` (`a, b` non-adjacent) force two distinct facets through `a` and
    `b`. They are legal only if they share a 2-face containing `a` and `b`. That face's angle at `a` is the edge
    common to the two germs, so the clause carries those faces as escape literals.
  - **F3, F5.** These exclude the whole facet (all its 2-faces and germs). That is equivalent to the facet being present.

**Sanity checks (all pass, §6.1).**
- Known 4-polytopes are accepted with genuine facet structures, e.g. pentagon×pentagon: 10 pentagonal prisms;
  rectified 5-cell: 5 tetrahedra and 5 octahedra.
- `K₃,₃×K₂` (PPS Prop. 2.11) and `P×K₂` (PPS Thm 2.3) are rejected.

## 6. Results

### 6.1 Calibration

| graph | expected | result |
|---|---|---|
| `C₅×C₅` | 4-polytope | SAT, 25 squares + 10 pentagons, 10 prism facets |
| `J(5,2)` (rectified 5-cell) | 4-polytope | SAT, 30 triangles, 5 tetrahedra + 5 octahedra |
| icosahedron × `K₂` | 4-polytope | SAT (passes every check) |
| `K₃,₃×K₂` | not polytopal | UNSAT |
| `P×K₂` | not polytopal | UNSAT |
| `K₃,₃×K₃`, without (H) | PPS: unique manifold is a cellular `S¹×RP²` | SAT: 6 triangles + 36 squares; 6 triangular prisms + 9 cubes (PPS's complex) |
| `K₃,₃×K₃`, with (H) | PPS: not 4-polytopal | UNSAT (all 456 induced cycles) |

### 6.2 Petersen × triangle is not polytopal

**Computed result.** The graph `P×C₃` (30 vertices, 5-regular) is not the graph of any polytope.
- **Dimension 4.** All induced cycles of `P×C₃` were enumerated: 7681, with at most 12 vertices; the longest wind
  around the triangle. They are all the candidate 2-faces, so the model is complete. It is UNSAT after 263 lazy
  iterations (21 minutes): 27 homology cuts and 1916 facet cuts (235 F1, 10 F4, 1671 F6).
- **Dimension 5.** A 5-polytope with a 5-regular graph is simple. PPS Thm 2.3 (a product is simply polytopal iff its
  factors are) excludes it, since `P` is not polytopal.
- **Dimensions ≤ 3.** Excluded because the graph is non-planar; no dimension above 5 is possible because the graph is
  5-regular.

PPS cite `P×C₃` as an example where local data cannot decide polytopality, since it is the graph of a cellular `S¹×RP²`.
- That complex has every fibre a full hemi-dodecahedron.
- Here it is cut by the homology condition: `H₁(S¹×RP²) ≠ 0`. It is also forbidden by Lemma 2.
- What decides the question is the global structure (homology and the facets as 3-polytopes), exactly as PPS
  predicted.

The same run restricted to the 3301 cycles of length ≤ 10 is UNSAT in about 1–2 minutes. Its record is in the `.out`.

**Cross-checks.**
- A second run with the CGS clauses added follows a different search path. It is also UNSAT (281 iterations).
- The final formula, with the base clauses and every lazily added cut, is re-solved from scratch by a second solver
  (Glucose 4). For faces ≤ 10 it agrees: UNSAT, 419026 clauses (`pc3 10 --crosscheck`).

### 6.3 Two excluded sub-cases of P□P

Both are UNSAT with the static constraints alone, with no lazy cut needed. In both, "no bent faces" means every
2-face is a square or a fibre cycle (665 candidates).
1. **Every vertex lies in exactly 8 square 2-faces** (the 25 non-face squares partition the vertex set): UNSAT in
   about 2 minutes. Without (R) it is SAT, so the edge-figure orientation is what kills it.
2. **Every vertex figure is a triangulation all of whose 8 triangles contain both an H- and a V-neighbour** (all
   facet germs are prism germs): UNSAT in about 2 minutes.

### 6.4 What does not work, and where the difficulty sits

- **Local 2-face data is satisfiable.** With no bent faces, (I)+(V) alone has solutions, e.g. 180 squares, 65
  fibre pentagons and 11 fibre hexagons with polyhedral vertex figures everywhere (the `.out`'s solution). So `P□P` passes every purely local test in dimension 4
  at the level of 2-faces.
- **The full no-bent problem is hard but globally rigid.**
  - With (R), imposing (V)+(R) only at the vertices of 8 of the 10 H-fibres is SAT in 10–20 s; with 9 fibres, SAT
    in 260 s.
  - The full instance did not finish within 1–3 hours under CaDiCaL 1.5.3/1.9.5, Lingeling, Glucose and
    MapleChrono, nor with lex-leader symmetry breaking (20, 60 or 200 group elements), nor split into the 30
    vertex-figure orbits at one vertex.
  - Capping the number of non-face squares at 30, 40 or 50 had not finished after an hour either.
  - Such a jump from "nine fibres easy" to "ten fibres hard" is the signature of a global parity obstruction of the
    snark type in Lemma 1, which CDCL solvers handle badly.
- **Bent faces.** `P□P` has 225 squares, 240 + 200 fibre cycles, and then 3600, 6300, 14400 and 172440 induced
  cycles of lengths 7–10. The bounded explicit model with faces ≤ 8 did not finish its first SAT call.

## 7. Next steps

- **A human proof of the no-bent packing case** (§6.3.1). There the non-face squares define perfect matchings `M_v`
  (one per level) and `N_u` (one per position).
  - Each vertex figure is `K₃,₃` minus the edge `H_{M_v(u)}V_{N_u(v)}`, plus pure edges at both ends of that
    missing edge (25 types, listed by the script).
  - The rotation at a non-matching H-edge is forced, e.g. "the fibre face through `(g, M)` is opposite the square
    `(g, N_u(v))`". That propagates `N_u(v) = N_w(v)` along edges.
  - This is the natural route to a 1-factorization contradiction as in Lemma 1.
- **The full no-bent case.**
  - A parity invariant generalising Lemma 1 from prism vertex figures to all 1227 types would settle it.
  - The four fibre-family types of L5 and Lemma 2 cut the search further.
- **Bent faces need a structural lemma** (e.g. bounding the number of turns) before any exhaustive method can
  apply.
- **Dimension 5 is untouched.** There the `RP²×RP²` structure shows that local arguments must be global.

## 8. Reproduction

`python 04-computation/experiments/petersen_product_4polytope_20261001.py all` takes about 5 minutes. It runs:
- the encoding verification;
- the calibration suite;
- the `P□P` structure facts;
- the `K₃,₃×K₃` validation;
- the two excluded sub-cases;
- the 8-fibre relaxation;
- `P×C₃` with faces ≤ 10.

`... all-long` runs `P×C₃` with all 7681 induced cycles (20–40 minutes). Both records are in the `.out`.

Requirements: `python-sat` (CaDiCaL 1.9.5) and `networkx`.
