# Even-graph equinumerosity: one cube, four faces, and eight wild lines

*kind-pasteur-2026-07-01-S11. Ideation, grounded. The number `2^{C(n-1,2)}` is a crossroads: it counts, simultaneously and via explicit labeled bijections, the tournament cycle-space (tilings), the even graphs on `n` vertices, the two-graphs on `n` vertices, and the graphs on `n−1` vertices. Yet the unlabeled (orbit) counts fan out — `A002854` (even = two-graph), `A000088(n−1)` (graphs), `A000568` (tournaments) — because the bijections are not equivariant for the competing group actions. This reflection maps the crossroads and proposes wild directions.*

## The grounded core: one cube, four faces

Labeled, all four are `2^{C(n-1,2)}` (verified n≤7): tilings/tournament-cycle-space, even graphs on `n`, two-graphs on `n`, graphs on `n−1`. The bijections are concrete — e.g. **graph `H` on `{1..n−1}` ↔ even graph on `{1..n}`** by adding an apex joined to the odd-degree vertices of `H` (the *apex/odd-set* bijection); and **tiling → XOR of fundamental cycles → even graph** (the repo's cycle-space map).

Unlabeled, they **diverge** (verified):

| n | even = two-graph `A002854` | graphs(n−1) `A000088` | tournaments `A000568` |
|---|---|---|---|
| 3 | 2 | 2 | 2 |
| 4 | 3 | 4 | 4 |
| 5 | 7 | 11 | 12 |
| 6 | 16 | 34 | 56 |
| 7 | 54 | 156 | 456 |

Two facts do the work: **(i) even graphs and two-graphs are the *same* count `A002854`** (both = switching classes of the cocycle space under `S_n`); **(ii) the apex bijection is not `S_n`-equivariant** (it breaks `S_n → S_{n-1}`), so `A002854(n) ≠ A000088(n−1)`. The *groupoid* cardinalities `Σ 1/|Aut|` are `2^{C(n-1,2)}/n!` vs `2^{C(n-1,2)}/(n−1)!` — differing by exactly the factor `n` (the apex's `n` placements). That factor-of-`n` is the whole anomaly, seen weighted.

Grounded bonus (verified): **`Σ_{even S⊆E(K_n)} x^{|S|} = 2^{−n} Σ_{s∈{±1}ⁿ} ∏_{i<j}(1+x·s_i s_j)`** — the even-subgraph generating function *is* the Curie–Weiss Ising partition function (van der Waerden high-temperature expansion). For odd `n` this polynomial is **palindromic** (even-graph complementation, since `K_n` is itself even) — the even/odd-`n` split is the Ising `Z₂` symmetry.

## Eight wild lines

**1. Even graphs = two-graphs = the Paley fixed point (the recurring universal object).** `A002854` counts *two-graphs* = Seidel switching classes; the regular two-graphs come from **conference matrices = Paley**. So the maximally-symmetric node of the even-graph metagraph `E_n` is the *same Paley object* that is `G_n`'s H-maximizer and (last session) the depth-extremal `|Aut|=21` obstruction to iso-covering. **Paley is the fixed point where the tournament, even-graph, two-graph, and equiangular-line worlds coincide.** → *First step:* compute the Seidel spectrum of `E_7`'s densest/most-symmetric node; look for the `√7` two-eigenvalue (Paley) signature; ask whether the equiangular-line *relative/absolute bounds* translate into LRC covering-min bounds (OPEN-Q-108 as an equiangular-lines extremal problem).

**2. `G_n ↔ E_n` IS the Tutte chromatic ↔ flow duality.** Even subgraphs = support of nowhere-zero flows (flow polynomial); acyclic tournament orientations = the chromatic/reliability side. The project's two pillars are then the two faces `T(x,y) ↔ T(y,x)` of the `K_n` Tutte polynomial, quotiented by `S_n`. → *First step:* test whether `V(E_n)`, `V(G_n)` and the full-rank bridge matrix `B[tourn,even]` are `S_n`-symmetrizations of Tutte specializations of `K_n`; if `B` is a Tutte transfer, the whole metagraph correspondence is one polynomial.

**3. `E_n` = the Curie–Weiss Ising model (grounded).** With the identity above, `E_n`'s combinatorics = the mean-field ferromagnet on `K_n`. → *First step:* locate the Ising critical temperature in the `E_n` family; test whether the **perfection-breaking at n=7** (odd holes appear, `E_n` stops being chordal) is a mean-field *phase transition*, and whether the chromatic explosion `χ(E_n)=2,3,5,10,28` is its order parameter.

**4. The non-equivariance anomaly, computed by Burnside.** The divergence `A002854(n)` vs `A000088(n−1)` is a precise obstruction. → *First step:* Burnside-decompose both by permutation cycle type; the divergence should localize on specific cycle types (likely those touching the apex / the `n`-cycles). The per-cycle-type difference is the **"equinumerosity anomaly" character** — a new `S_n`-representation-theoretic invariant.

**5. `E_n` as a Cayley/Schreier graph of `F₂^{C(n-1,2)}` — a homology space.** Even graphs = `H₁(K_n; F₂)`; `E_n` is a graph *on* the 1-cycles. If its edges are "add a fundamental cycle," `E_n` is a Schreier graph of the cycle space and its spectrum is character-theoretic (Walsh/Fourier on `F₂^{C(n-1,2)}`). → *First step:* check whether `E_n` (or the tile-flip metagraph) is vertex-transitive on the labeled cover = a Cayley graph; if so, diagonalize by Walsh characters and read the spectrum off cycle weights.

**6. `q`-equinumerosity and reciprocity.** The palindromic even-graph polynomial (odd `n`) is a `q`-reciprocity in disguise. → *First step:* form the `S_n`-symmetrized `q`-generating function (cycle-index over even graphs weighted `q^{edges}`) and hunt for a `q ↔ 1/q` (or `q ↔ 1−q`, flow-polynomial) reciprocity linking even(n) to graph(n−1) `q`-analytically — a *weighted* equinumerosity the unweighted counts hide.

**7. The LRC "even runners" live in `E_7`.** opus's CRT split (HYP-3802): the `n=14` runners factor as **odd (heptagon/tournament `R_7`) × even (harmonic Verblunsky)**. The *even* half should have a combinatorial metagraph — conjecture: `E_7`, dual to the odd half's `R_7`. → *First step:* map the 6 even runners `{2,4,6,8,10,12}` to an even graph / two-graph on 7 vertices via the apex bijection; check whether its atom/switching structure reproduces the even-runner Verblunsky ladder `1/(6−j)`.

**8. A master `n=7/8` phase transition across both pillars.** Even-graph perfection breaks at **n=7**; tournament iso-covering breaks at **n=6**; the max-cut/entropy crossover at **n=8**. → *First step:* define one order parameter (e.g. automorphism entropy `⟨log|Aut|⟩`, or the metagraph spectral gap) and track it across `G_n` and `E_n` together; test whether the three thresholds are one inflection in the shared cube `Q_{C(n-1,2)}` as `n` crosses 7.

## Why this is worth chasing

The equinumerosity is a *coincidence at the labeled level that four different symmetry groups then break in four different ways* — and the breakings are exactly the objects the project already circles (Paley, the `n=7/8` wall, the even/odd split, the H-gradient). Lines 1–3 are grounded enough to start now; 4–6 are precise and computable; 7–8 tie the wild ideas back to the LRC finish and the triangle foundation. The unifying bet: **the even-graph face is the flow/Ising/two-graph dual of the tournament face, and both are pinned at Paley** — so "even graph equinumerosity" is the statement that the project has been studying *one* self-dual object (the `K_n` Tutte/Ising system mod `S_n`) from two sides.

— Related: `even-graphs-as-first-class.md`, `the-even-graph-is-the-tournaments-cycle-half.md`, `even-graphs-through-the-metagraph.md`, `reed-muller-on-the-tiling-cube.md` (Walsh/cube spectrum), last session's `two-axes-of-the-tournament-metagraph-…` (Paley as the shared extremum), HYP-3802 (CRT parity split / even runners), `everything-is-the-triangle.md`. Script: `04-computation/even_graph_equinumerosity_probes_kps.py` (+ .out). Not a HYP reservation (ideation); candidate investigations for INVESTIGATION-BACKLOG.
