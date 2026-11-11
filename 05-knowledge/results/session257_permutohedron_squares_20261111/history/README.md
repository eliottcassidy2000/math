# Recovered source atlas: the actual cycle, colors, primes, and permutohedron

Session257, 2026-11-11. Inheritance and synthesis proposals only; not new theorem claims. The shared checkout was read without editing it. `SOURCE_PINS.json` and `sources/` retain exact source versions. Source paths below are repository-relative. Startup inheritance: AGENTS.md, START-HERE, ACTIVE-GUARDRAILS and RESEARCH-PROTOCOL; corrected claims outrank historical prompts. Session256/packet/SYNTHESIS.md was read in full.

## 1. Exact meanings recovered

### The cycle is a valuation display, not a square-graph vertex set

`05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`, §§0,1,3, records the original October23 owner prompt and the exact integer

N = 2^5 3^7 5^3 7^4 11 13^2 = 39,046,339,332,000.

The arrows are 13→2→5→3→7→4 and 11→1, namely p→v_p(N). Identifying the terminal4 with initial13 produces **five** quotient states, not six. This is a valid directed-graph quotient but not integer equality or a polynomial dynamical system. Session32 `VALUATION_CYCLE_AND_QUOTIENT.md` proves no nontrivial quotient ring identifying13 and4 supports the whole map. Even the selected cycle is not an integer-polynomial map mod9. Under the actual 2A stabilizer normalization N/2 the five-cycle becomes a two-cycle. The exact group valuation comparison is meaningful; dynamical significance of the identification remains an analogy.

### The two prime sets are genuine group/level shells

The fifteen Ogg primes split into Fi22 support {2,3,5,7,11,13}, Baby-minus-Fi22 {17,19,23,31,47}, Monster-minus-Baby {29,41,59,71}. They are also respectively positive valuation excess, equal nonzero valuations, and Monster-only support in the two group orders. The modular realization is the rationality filtration of X_0(1), X_0^+(2), X_0^+(6), adjoining primes dividing the level. This is CITED mathematics with independently checked finite computations in the source, not an inference from numerical pattern. Baby Monster is not itself a subgroup of Monster: the stabilizer is 2.B. Four absent primes each have one Frobenius pair of nonrational coarse level2 points; their underlying j vertices are rational and Frobenius swaps parallel 2-isogeny edges.

### The requested three colors are chronological birth classes

The precise original user edge list is reproduced and decoded in `05-knowledge/results/session40_marked_ports_square_sum_20261023/square-sum/README.md`, §2. For edge a–b, a<b, birth means arrival b:

1. Traceable birth orders15,16,17,23,25: edges (1,15),(10,15),(9,16),(8,17),(2,23),(13,23),(11,25),(24,25).
2. Failed double-neighbor births13,14,19,20,21,22,24: (3,13),(12,13),(2,14),(11,14),(6,19),(17,19),(5,20),(16,20),(4,21),(15,21),(3,22),(14,22),(1,24),(12,24).
3. All remaining edges, including7–18 and the earlier scaffolding.

This is neither a proper edge coloring nor forced/optional/unused. At25, color1 alone contains unused1–15, optional13–23, forced11–25. The original literal3–12 is a nonedge; source explicitly does not guess its replacement. Extension beyond25 is a declared model, not a newly attributed instruction: traceable birth / failed degree2 birth / remainder.

Do not confuse this palette with the square labels q=a+b (partial-reflection matchings), the three Q12 components, or the distinct Fibonacci three-hidden-ones palette. An additional October23 prompt explicitly proposes island allegiance colors and Berggren branches; that is a proposed investigation, not an established identification.

## 2. Source/target/map/invariant/loss/sidecar atlas

### A. Exact edge-polynomial birth transfer [PROVED inherited]

Source: Session40 square-sum README §6. Source object: weighted graph G plus neighbor set N of a new vertex v. Target: Hamilton-path generating polynomial of G+v. Map:

H_(G+v) = Σ_(a∈N) x_va E_a(G) + Σ_({a,b}⊂N) x_va x_vb F_ab(G).

Invariant: every actual path and its multiaffine edge monomial, with unit multiplicity. E_a counts old paths exposing a; F_ab counts spanning two-path forests with a,b endpoints in DIFFERENT components. Lost by scalar H: all resurrection terms, notably H(Q24)=0 but F_11,24(Q24)=10. Sidecar: component-owned named endpoints, singleton convention, and higher-component forest states for future iteration. This is the immediate predecessor of any new colored seam polynomial; do not reinvent it as new.

### B. Named block transfer [PROVED + FINITE-EXACT inherited]

Source: Session32 `ports/PORT_TRANSFER.md`. Source: Q25's induced blocks A={2,7,9,16,18}, B={4,5,11,12}, C={1,3,6,8,10,17}. Target: boundary continuation tensors. Map: enumerate internal path forests, residual endpoint demands and component pairings, then glue actual boundary edges. Invariant: exact Hamiltonian feasibility and counts. Data: 5,101,149 local states;1,543 endpoint-budget triples;19 degree-feasible assemblies;9 rejected cyclic assemblies;10 paths survive. Lost by whole-block contraction: Hamiltonicity (quotient has no path although Q25 does), despite preserved multigraph cycle rank8. Sidecar: named ports, pairings, global connectedness. Block A is branched, not a path. This is the closest prior for two-/four-seam gluing.

### C. Color specialization and exact reconstruction [PROVED/FINITE-EXACT inherited]

Source: Session40 square-sum README and COLOR_FIRST_FAILURE.md. Source: full edge polynomial or actual path; target: three-color monomial, histogram, or ordered word. Map: substitute x_e=z_color(e), or read colors in path order. Invariant: color-use statistics. Lost: actual incidence and endpoints' component ownership. Exact endpoints+histogram first fail26; endpoints+ordered word first fail36 (12,960 paths,12,957 labels). Sidecar: full numerical square-sum word plus initial vertex, or endpoints plus cotree edge bits and degree/connectivity checks. Cheap positive consumer: specialize inherited edge polynomial only AFTER gluing, rather than glue color marginals.

### D. Matrix-tree and covariance interfaces [PROVED inherited]

Source: Session40 `audit-first/MATRIX_TREE_PATH_PROJECTION.md` and `PATH_SCATTER_AND_MARKING.md`. Source: edge-weighted spanning-tree polynomial. Target: Hamilton-path polynomial by degree projection with endpoint marking. Invariant: connectedness and shared-edge derivatives; positive path ensemble gives covariance/conditional Hessian. Lost by separate marginals: joint realizability (at25 edges1–24 and4–5 appear individually but never together). Sidecar: full multiaffine polynomial or realized shared-incidence state. A new Gram comparison must distinguish this positive ensemble covariance from an indefinite gluing bilinear form.

### E. Coherent permutohedral orientation [PROVED inherited]

Source: `session210/coherent-permutohedron/PROOF.md`, with older THM-207 `01-canon/theorems/THM-207-permutohedron-transitive.md`. Source: labelled tournament. Target: permutohedron with orientations coherent across each parallel edge class. Map: orient adjacent swaps using the tournament's comparison. Invariant: Hamilton paths correspond exactly to sinks; faces obey product sink laws. Lost by score projection or arbitrary orientation: literal comparison coherence. Sidecar: comparison attached to each unordered label pair. Generic convex-linear objectives have one sink; arbitrary acyclic orientations admit many more counts and are not the tournament carrier. Most promising square-sum consumer: keep numerical labels and square adjacency as a filter on permutation vertices, while keeping coherent tournament orientations separately.

### F. Score permutohedron and strong gate [PROVED inherited]

Source: `session211/score-permutohedron/PROOF.md`. Source: tournament plus literal two-vertex gate. Target: score-zonotope point and strong-component face. Map: take scores / gate translation. Invariant: ordered strong-component face recovered from tight inequalities; gate moves source/sink face into interior with margin1. Lost: internal arcs and gate identity. Equal-score strong tournaments can have gated-cycle counts11 versus7,37 versus21, without violating the Hamilton-PATH forbidden counts. Sidecar: coherent orientation and recognizable gate. Do not infer Chan–Pak from interior score geometry.

### G. Group valuation/action interface [CITED + elementary consequences]

Source: Monster wave2 note and Session32 valuation note. Source: Monster action on2A involutions, stabilizer2.B. Target: prime fixed-point strata / valuation vector. Map: compare Sylow orders and stabilizer divisibility. Invariant: full-Sylow fixed, partial order-p fixed, or no order-p fixed behavior. Lost by prime support: valuation multiplicities; lost by exponent-label quotient: prime/composite distinction and normalization. Sidecar: group action, stabilizer, actual valuations. Promising consumer: induced action on labelled permutations/permutohedral faces, but no map to square-sum adjacency is inherited.

### H. Pythagorean lift and contextual response [Session256 inherited]

Source: Session256 spinor/PROOF.md, compiler/PROOF.md, gram/PROOF.md. Source: primitive Euclidean coordinates plus Berggren address. Target: actual apex-over-transitive-word tournament, or sorted C3[TT_x,TT_y,K1]. Map: A/B/C substitutions O→OI,I→I; O→OI,I→OIO; O→O,I→OOI. Invariant: counts of letters recover spinor/triple; full positive matrix state computes H. Lost by sorting: cut order and attainable residues (morphic CA has73, sorted107). Sidecar: word/base order/apex or full matrix response. For cyclic gluing with fixed observerU, truncated ordered path-cover profiles give a positive Gram kernel rank|U|+1. AtU=K1 norm is x²+2y², not Gaussian x²+y². This offers actual operation-compatible comparison to B/D, not merely three branches versus three colors.

## 3. Current status and nonduplicating next tests

- Chan–Pak completeness remains OPEN. Prior sources already establish the absence7,21; translating it to permutohedra is not a new proof.
- The old Session40 README §7 says prescribed-edge eventual behavior remains open, but its **later SYNTHESIS §5 and LATE_INHERITANCE.md supersede that coarse question**: chord insertion plus the cited all-n≥32 cycle theorem makes every edge optional for n≥max(32,a,b). Retain finer endpoint/paired-incidence questions rather than rediscover this coarse result.
- Best new target: lift the birth-color edge polynomial to coherent permutohedral sinks with numerical square adjacency retained, compare operation-specific two-seam/four-seam response forms, and test whether group fixed-point action data supplies a lawful symmetry on that same labelled carrier.
- Keep the 4↔13 quotient as a finite-state clock only if useful; no canonical arithmetic or square-graph transport has been established.

## 4. Provenance and reproducibility scope

Source documents are pinned in SOURCE_PINS.json with commit URLs; this atlas does not rerun their extensive historical censuses. Explicit current user history identifies October23 original graph/valuation prompt, October23 birth-color list, and November7 micro/macro permutohedron/gate request. History confirms intended objects; historical mathematical status comes from the repaired source documents, not from the prompts. No private attachment-derived Session32 lane is imported or published here.
