# Is the product of two Petersen graphs 4-polytopal? An attempt: structure, intrinsic linking, exact SAT exclusions, and Petersen × triangle is not polytopal

**opus S15, 2026-10-01.**

- Script: `04-computation/experiments/petersen_product_4polytope_20261001.py` (+ `.out`).
- Previous note: [`odd_zeta_parallels_petersen_20261001.md`](odd_zeta_parallels_petersen_20261001.md) §1 (status of
  Ziegler's question; the "≤ 8 mixed 2-faces per vertex" remark).

**Status.**
- **Ziegler's question in dimension 4: still OPEN.** No proof that `P□P` is not 4-polytopal is claimed.
- PROVED:
  - L1, L3, L4: elementary.
  - L2 and L5: by hand plus exhaustive enumeration, independently reproduced.
  - Lemma 1: elementary. Its hypotheses are already inconsistent with L1, which retracts the first version's "snark
    forbids it" reading (§3, MISTAKE-558).
  - Lemma 2: from the Conway–Gordon–Sachs theorem.
- **COMPUTER-ASSISTED: `P×C₃` (Petersen × triangle) is not polytopal** (§6.2, THM-4535).
  - It is an exact SAT search over all 7681 induced cycles as candidate 2-faces.
  - The cut layer is custom code whose validity is argued in §5.
  - The UNSAT verdict was **reproduced by an independent re-implementation** (audit, §9). Three SAT solvers agree on
    the final formulas.
  - There is no DRAT/LRAT certificate.
  - PPS cite this graph as one that local data cannot decide: it is the graph of a cellular `S¹×RP²`.
- Validation of the pipeline: without the homology condition, it finds a cellular `S¹×RP²` with graph `K₃,₃×K₃`
  (6 triangular prisms + 9 cubes); it rejects that complex when the condition is kept. This is consistent with PPS's
  announced, unpublished claim that this `S¹×RP²` is the unique strongly regular combinatorial manifold with graph
  `K₃,₃×K₃`.
- COMPUTED: two sub-cases of `P□P` are impossible (§6.3).
- COMPUTED (negative, no-bent model only):
  - (I)+(V) is satisfiable.
  - (I)+(V)+(R) imposed at the vertices of 8 or 9 of the 10 H-fibres is satisfiable, also with orientability (O)
    added.
  - So any obstruction in the no-bent model must involve vertex conditions in every H-fibre (§6.4).
- Independent audit: DONE (2026-10-01, blind, own code; §9). The headline is reproduced. Lemma 1's commentary is
  corrected, wording overclaims are removed, and a non-termination bug in the F6 cut is fixed. MISTAKE-558.

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
  - The edge figure at `xy` is a polygon. The 2-faces through `xy` are cyclically ordered the same way, up to
    reversal, as seen from `x` (the rotation at `y` in `Γ_x`) and from `y`. "Consecutive" means "in a common facet".
  - `∂Π` is a 3-sphere, so `H₁(∂Π; Z/2) = 0`, and the 2-faces' boundaries span the cycle space of `G`
    (dimension 201). `∂Π` is orientable.
- Pfeifle–Pilaud–Santos (PPS; arXiv:1009.1499v1, whose theorem numbering is used here; Israel J. Math. 192 (2012))
  excluded dimension 6 and left 4 and 5 open.
  - Their remark that local data cannot decide refers to the 4-manifold `RP²×RP²`, i.e. to dimension 5.
  - **In dimension 4 the vertex links must be 2-spheres on 6 points. That is a strong local constraint, and this note
    tests it.**

## 2. Elementary structure of a hypothetical realization

- **L1 (no facet lies in a fibre).** A proper induced subgraph of `P` has a vertex of degree ≤ 2, and `P` itself is
  non-planar.
- **L2 (both directions go straight somewhere at every vertex).** `Γ_x` has at least one pure H edge (a 2-face goes
  straight through `x` along the H-fibre), at least one pure V edge, and at most 8 of the 9 mixed edges.
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
  pairwise properly meeting ones lie in one of **28 maximal families, of three kinds** (computed):
  - the 6 face pentagons of either hemi-dodecahedral embedding `E₁`, `E₂` of `P` in `RP²` (2 families);
  - a complementary (vertex-disjoint) pair of pentagons (6 families);
  - a hexagon `C_a = P − N[a]` with the 3 pentagons through `a` from one embedding (20 families).
  - In particular, at most one hexagon per fibre. Each complementary pair has one pentagon in `E₁` and one in `E₂`.

## 3. The prism scenario and the snark property

**Lemma 1 (a 2-face-level obstruction).** Suppose:
- every vertex figure is a triangular prism whose triangles are `{H₁,H₂,H₃}` and `{V₁,V₂,V₃}`;
- every rung (every mixed angle that lies in a 2-face) lies in a square 2-face.

Then the Petersen graph is 3-edge-colourable, which is false.

*Proof.*
- At `x = (u,v)` the three rungs pair the edges at `u` bijectively with the edges at `v`. By L4, the 2-face through
  the rung `H_e–x–V_f` is the square `e×f`.
- For an edge `f` of the second factor, let `M_f = {e : e×f is a 2-face}`. At every `x = (u,v)` with `v ∈ f`,
  exactly one edge `e ∋ u` lies in `M_f`, so `M_f` is a perfect matching.
- For the three edges `f₁, f₂, f₃` at `v`, the matchings are disjoint (each `e ∋ u` has exactly one rung), so
  `M_{f₁} ∪ M_{f₂} ∪ M_{f₃} = E(P)`. That is a 1-factorization. ∎

**Remark (correction of the first version).** For 4-polytopes, the hypotheses of Lemma 1 already contradict L1.
- Let `F` be the facet whose germ at `x` is the triangle `{H₁,H₂,H₃}`. Its two 2-faces through `xH₁` are not squares.
- Hence their angles at `H₁` are pure, and `F`'s germ at `H₁` is again the H-triangle. By connectivity, `F`'s graph
  is the whole fibre `P`, contradicting L1.
- So the snark property is **not** what excludes this picture for 4-polytopes. Lemma 1 shows that the picture already
  fails at the level of 2-faces, before facets are considered.
- The same argument shows: for a cubic graph `H` tiling a surface other than the sphere, this "colour-matched" picture
  is never the 2-skeleton of a cellular 3-manifold.
- The first version of this note said the picture "is exactly what Petersen's being a snark forbids". It also said
  the recipe gives a candidate cell structure for `H□H` for any 3-edge-colourable `H` tiling a surface. Both statements
  are withdrawn (MISTAKE-558).

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

More generally, a pair has even linking number as soon as one member bounds a `Z/2` 2-chain of 2-faces avoiding the
other. An example is an annulus of squares plus a cap one level up.

**Observation (computed).** Let `β` be a symmetric bilinear form on `Z₁(P; Z/2)`. Then the sum `Φ(β)` of `β(C, C')`
over the 6 complementary pairs is 0.
- `β_{v,w}(Z, Z') = lk(Z×{v}, Z'×{w})` (mod 2) is a bilinear form for distinct fibres.
- Suppose all 15 squares `e×vw` are 2-faces (a "fully faced ladder", `vw` an edge). Then each `C×{w}` can be slid to
  `C×{v}` without meeting `C'×{w}`, so Lemma 2 makes `Φ(β_{v,w})` odd.
- Hence **the linking form across a fully faced ladder is never symmetric.** Since `Φ` vanishes on symmetric forms,
  it factors through `β ↦ β + βᵀ`, an alternating form. Over `Z/2` there is no symmetric/antisymmetric splitting.
  (Heuristically, `β + βᵀ` records how the three squares are ordered around each vertical edge; this is not proved
  here.)

## 5. The SAT model (all constraints are necessary conditions)

**Variables.**
- `F_Q`: candidate 2-face `Q` (an induced cycle of `G`) is chosen.
- `A(x;a,b)`: the angle `a–x–b` is covered by a chosen face.
- `M(x,y;z,w)`: the 3-path `z–x–y–w` lies in a chosen face.
- `C(x,y;z,z')`: `z, z'` are consecutive around `y` in `Γ_x`.

**Static constraints.** All encodings assume maximum degree ≤ 6 (asserted in the code).
- **(I) Proper intersection.** At most one chosen face contains any given *non-adjacent* vertex pair. This exactly
  captures improper meeting of two induced cycles.
- **(V) Vertex figures.** `Γ_x` is 3-connected and planar. The CNF (separations of ≤ 2 vertices; no `K₅`, subdivided
  `K₅` or `K₃,₃` on 6 vertices) was checked against brute force: 1, 25, 1227 labelled graphs on 4, 5, 6 vertices.
- **(R) Edge figures.** `C` is defined from the `A`'s by Tutte's theorem: in a polyhedral graph the faces are the
  induced non-separating cycles. The resulting DNF was checked against planar embeddings on all 1227 + 25 + 1 graphs.
  Two faces through `xy` must be consecutive around `y` in `Γ_x` iff they are consecutive around `x` in `Γ_y`.
- **(O) Orientability (optional; added after the audit, used only in §6.4 and the validations).**
  - Every `Γ_x` carries an oriented rotation system (successor variables). It must be traced coherently along every
    face of `Γ_x`.
  - Across every edge the cyclic order of the 2-faces seen from one end is the reverse of the order seen from the
    other.
  - Validation: genuine 4-polytopes are accepted, and an independent orientation check of their face lattices finds
    no conflict. PPS's `S¹×RP²` for `K₃,₃×K₃` is non-orientable, and with (O) the `K₃,₃×K₃` instance is UNSAT even
    without (H).

**Lazy constraints (exact cuts).**
- **(H) Homology.** If the chosen faces do not span `Z₁(G; Z/2)`, find a cochain `w` that vanishes on them but not on
  `Z₁`. Add the clause "some face with odd `w`-weight is chosen". It is valid because a sphere's 2-faces span `Z₁`.
- **(F) Facets.**
  - A *germ* is a face of `Γ_x`. Germs at `x` and `y` glue when they contain the same two consecutive 2-faces through
    `xy`. A facet is a glued class.
  - Each class must be a 3-polytope: one germ per vertex (F1), induced (F4), Euler characteristic 2 (F3),
    3-connected (F5). Two facets must meet in a common face (F6).
  - A germ is encoded exactly by "induced and non-separating in `Γ_x`".

**Soundness of the facet cuts.**
- **Glue step.** Suppose the two 2-faces `Q₁, Q₂` and the germs at `x` and `y` that they glue are all present. Then
  the facets of the two germs both contain `Q₁ ∪ Q₂`. Two 2-faces through an edge lie in at most one common facet,
  so these two germs belong to one facet.
- **F1, F4.** The literals along a path of glued germs force both ends into one facet. So a repeated vertex (F1), or a
  missing edge `ab` with `a ~ b` in `G` (F4), excludes the path's literals.
- **F6.** Paths in both classes from `a` to `b` (`a, b` non-adjacent) force two distinct facets through `a` and
  `b`. These are legal only if they share a 2-face containing `a` and `b`. That face's angle at `a` is the edge common
  to the two germs, so the clause carries those faces as escape literals.
  - The witness pair `a, b` is chosen outside every 2-face shared by the two classes. Otherwise the cut can fail to
    exclude the current assignment.
  - The first version did not do this. The audit reproduced a stall: the same assignment repeated, with no effect on
    soundness. If no such pair exists, the whole pair of classes is cut.
- **F3, F5, and whole-class cuts (closure argument).** Suppose all 2-faces and germs of a class hold.
  - By the glue step along a spanning tree, all its germs lie in one facet `F'`.
  - The class is closed under gluing and `F'`'s graph is connected. So `F'` has exactly the class's vertices and
    2-faces.
  - Hence `F'` would have the forbidden property: `χ ≠ 2`, not 3-connected, not induced, or meeting another such
    class improperly.
- An iteration cap is reported as INCONCLUSIVE, never as SAT.

**Sanity checks (all pass, §6.1).**
- Known 4-polytopes are accepted with genuine facet structures. For example, pentagon×pentagon gives 10 pentagonal
  prisms, and the rectified 5-cell gives 5 tetrahedra and 5 octahedra.
- `K₃,₃×K₂` (PPS Prop. 2.11) is rejected. So is `P×K₂`, by PPS Prop. 2.13, or by Thm 2.3 plus non-planarity.

## 6. Results

### 6.1 Calibration

| graph | expected | result |
|---|---|---|
| `C₅×C₅` | 4-polytope | SAT, 25 squares + 10 pentagons, 10 prism facets |
| `J(5,2)` (rectified 5-cell) | 4-polytope | SAT, 30 triangles, 5 tetrahedra + 5 octahedra |
| icosahedron × `K₂` | 4-polytope | SAT (passes every check) |
| `K₃,₃×K₂` | not polytopal | UNSAT |
| `P×K₂` | not polytopal | UNSAT |
| `K₃,₃×K₃`, without (H) | PPS: a cellular `S¹×RP²` | SAT: 6 triangles + 36 squares; 6 triangular prisms + 9 cubes (PPS's complex) |
| `K₃,₃×K₃`, with (H) | PPS: implied by their announced (unpublished) uniqueness claim | UNSAT (all 456 induced cycles) |
| `K₃,₃×K₃`, without (H), with (O) | the complex above is non-orientable | UNSAT |

### 6.2 Petersen × triangle is not polytopal

**Computer-assisted result.** The graph `P×C₃` (30 vertices, 5-regular) is not the graph of any polytope.
- **Dimension 4.** All induced cycles of `P×C₃` were enumerated: 7681, with at most 12 vertices.
  - Of the 2580 longest, 2220 wind around the triangle (2160 once, 60 twice) and 360 do not.
  - The induced cycles are all the candidate 2-faces, so the model is complete.
  - Recorded run (`.out`, corrected code): UNSAT after 184 lazy iterations (10 minutes). Cuts: 26 homology cuts and
    1241 facet cuts (157 F1, 6 F4, 1078 F6).
  - The final formula (1,394,307 clauses) is re-solved from scratch by Glucose 4: UNSAT (15 minutes).
- **Dimension 5.** A 5-polytope with a 5-regular graph is simple. PPS Thm 2.3 (a product is simply polytopal iff its
  factors are) excludes it, since `P` is not polytopal.
- **Dimensions ≤ 3.** Excluded because the graph is non-planar; no dimension above 5 is possible because the graph is
  5-regular.

PPS cite `P×C₃` as an example where local data cannot decide polytopality, since it is the graph of a cellular `S¹×RP²`.
- That complex has every fibre a full hemi-dodecahedron, with 18 pentagonal prisms as facets.
- It passes all local and facet conditions and is cut only by the homology condition: `H₁(S¹×RP²) ≠ 0`. It is also
  forbidden by the argument of Lemma 2 applied to a fibre `P×{c}`, and it is non-orientable.
- **Both global layers are indispensable**, as the audit checked:
  - static + (H) without facet checks is satisfiable (1 triangle, 45 squares, 13 pentagons, 1 hexagon, spanning
    `Z₁`);
  - static + facets without (H) is satisfiable (the `S¹×RP²` complex and several variants).
- So homology and the facets together decide what local data cannot.

The same model restricted to the 3301 cycles of length ≤ 10 is UNSAT after 175 iterations. Its final formula
(418,764 clauses) is re-solved by Glucose 4: UNSAT. Both records are in the `.out`.

**What has been checked, and what the checks cover.**
- **The solver.** Two SAT solvers agree on the final formula. This checks the solver, not the encoding or the cut
  generator.
- **The encoding and the cut generator.** These were checked by the audit's **independent re-implementation** (§9).
  - It shares no code with the script: a different Petersen labelling, its own enumeration of the 1227 polyhedral
    graphs, selector-based (V) and (R), its own homology and facet code.
  - All 7681 candidates: UNSAT after 225 iterations. Glucose 4 and MapleChrono agree on its final formula
    (1,088,001 clauses).
  - It is also UNSAT when every F6 cut is a whole-class cut, i.e. without the escape-literal argument: 457 iterations.
  - All 10,280 clauses the two pipelines generated on `P×C₃` with (H) off are satisfied by PPS's `S¹×RP²` complex.
  - On 29 genuine 4-polytopes, the true face lattices satisfy every recorded clause.

### 6.3 Two excluded sub-cases of P□P

Both are UNSAT with the static constraints alone, with no lazy cut needed. In both, "no bent faces" means every
2-face is a square or a fibre cycle (665 candidates).
1. **Every vertex lies in exactly 8 square 2-faces** (the 25 non-face squares partition the vertex set).
   - UNSAT in about 1 minute.
   - Without (R) it is SAT, e.g. 200 squares, 48 pentagons and 15 hexagons (`pp-packing-norot` in the `.out`). So
     the edge-figure condition is what kills it.
2. **Every vertex figure is a triangulation all of whose 8 triangles contain both an H- and a V-neighbour** (all
   facet germs are prism germs): UNSAT in about 1 minute.

The audit reproduced both independently, and also the SAT result without (R).

### 6.4 What does not work, and where the difficulty sits

All statements here are about the no-bent model.
- **(I)+(V) alone is satisfiable.** For example, 180 squares, 65 fibre pentagons and 11 fibre hexagons give polyhedral
  vertex figures everywhere (the `.out`'s solution). (R) is also a local condition, and with (R) the question is open.
- **The full no-bent problem with (R) is open: no solver finished.**
  - Imposing (V)+(R) only at the vertices of 8 of the 10 H-fibres is SAT in about 10 s. With 9 fibres it is SAT in
    125 s (both in the `.out`). With 9 fibres plus 3 vertices of the tenth it is still SAT, but took about an hour
    (author's run, reproducible as `pp-levels 9 3`).
  - Adding orientability (O) keeps these relaxations satisfiable:
    - 8 fibres: SAT in 136 s (in the `.out`; an independent orientation check of the witness finds no conflict);
    - 5, 6, 7 and 9 fibres: SAT in 36 s, 4 s, 50 s and 235 s (author's runs of `pp-levels-orient k`, not in the
      `.out`).
  - The full instance did not finish within 1–3 hours (over an hour with (O)) under any of these:
    - CaDiCaL 1.5.3/1.9.5, Lingeling, Glucose or MapleChrono;
    - lex-leader symmetry breaking (20, 60 or 200 group elements);
    - a split into the 30 vertex-figure orbits at one vertex;
    - a cap of 30, 40 or 50 on the number of non-face squares;
    - added (O).
  - Speculation: the jump from "nine fibres easy" to "ten fibres hard" could come from a global parity obstruction,
    which CDCL solvers handle badly. It is equally consistent with a satisfiable but hard instance.
- **Bent faces.** `P□P` has 225 squares, 240 + 200 fibre cycles, and then 3600, 6300, 14400 and 172440 induced
  cycles of lengths 7–10. The bounded explicit model with faces ≤ 8 did not finish its first SAT call.

## 7. Next steps

- **A human proof of the no-bent packing case** (§6.3.1). There the non-face squares define perfect matchings `M_v`
  (one per level) and `N_u` (one per position).
  - Each vertex figure is `K₃,₃` minus the edge `H_{M_v(u)}V_{N_u(v)}`, plus at least one pure edge at each end of
    that missing edge, and possibly others.
  - There are 25 labelled types per missing edge, 16 of which have further pure edges. The script does not list them.
  - The rotation at a non-matching H-edge is conjecturally forced, e.g. "the fibre face through `(g, M)` is opposite
    the square `(g, N_u(v))`". That would propagate `N_u(v) = N_w(v)` along edges.
- **The full no-bent case.**
  - The three fibre-family kinds of L5 and Lemma 2 cut the search further.
  - Orientability (O) is a further necessary condition already implemented.
- **Bent faces need a structural lemma** (e.g. bounding the number of turns) before any exhaustive method can
  apply.
- **Next products.** `P×C₄` has 132,578 induced cycles and `P×K₄` has 104,498, so the complete model is within reach.
  `P×K₄` was still running at the time of writing.
- **Dimension 5 is untouched.** There the `RP²×RP²` structure shows that local arguments must be global.

## 8. Reproduction

`python 04-computation/experiments/petersen_product_4polytope_20261001.py all` takes about 10 minutes. It runs:
- the encoding verification;
- the calibration suite;
- the `P□P` structure facts;
- the `K₃,₃×K₃` and orientability validations;
- the two excluded sub-cases and the packing case without (R);
- the 8-fibre relaxation;
- `P×C₃` with faces ≤ 10, with the second-solver re-solve.

`... all-long --crosscheck` takes about an hour. It adds the 9-fibre relaxation, the 8-fibre relaxation with (O), and
`P×C₃` with all 7681 induced cycles, with the Glucose 4 re-solve. Both records are in the `.out`.

Requirements: `python-sat` (CaDiCaL 1.9.5, Glucose 4) and `networkx`.

## 9. Audit record (2026-10-01)

A blind auditor subagent wrote its own code. The record is in the session scratchpad: `audit22/AUDIT.md`, all
scripts and outputs.
- **SOUND:**
  - the §1 standard facts, including the reading of PPS's `RP²×RP²` remark as concerning dimension 5;
  - L1, L3, L4;
  - L2 (1227 graphs and the census reproduced exactly);
  - Lemma 2 (CGS correctly stated; parity invariance checked; 2000/2000 random embeddings odd);
  - (I), (V), (R), (H) and all facet cuts. The CNF and DNF match the auditor's brute force; the closure argument
    above was requested and added.
- **SOUND WITH CORRECTION:**
  - L5: 28 maximal families of three kinds.
  - Lemma 1: the hypothesis says "every rung", and the scenario is already excluded by L1. The snark commentary is
    withdrawn.
  - The Observation: no "antisymmetric part" over `Z/2`.
- **REPRODUCED independently:**
  - the `P×C₃` headline: enumeration identical under an explicit isomorphism; UNSAT after 225 iterations; Glucose 4
    and MapleChrono agree; also UNSAT without the F6 escape argument;
  - the `K₃,₃×K₃` validation;
  - both `P□P` sub-cases, and the packing case without (R);
  - the levels-8/9 relaxations, whose witnesses were checked by the auditor's own checker.
- **BUG (not soundness).** The F6 witness pair could lie inside a shared 2-face, so the loop stalled. It is fixed.
  An iteration cap was printed as SAT; that is fixed too.
- **Wording corrections applied:**
  - §6.4 overclaims ("passes every purely local test", "globally rigid", "signature of a parity obstruction");
  - unrecorded runs are either now recorded or removed: the CGS-clause runs are dropped;
  - THM-4535 no longer calls reruns "independent";
  - the "longest cycles wind" claim is corrected;
  - the timings are corrected;
  - the PPS numbering source is stated.
