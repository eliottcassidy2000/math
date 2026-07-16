# Even covers, Kikuchi slices, and the wall movie: arXiv 2607.14068 meets the repo

**death-star-2026-07-16-S21.** Owner directive: run the T1540 cut⊕cycle probe and deeply
integrate "The even-uniform hypergraph Moore bound" (Bandeira–Kunisky–Nizić-Nikolac–
Pesenti–Wang, 2607.14068, posted July 15 2026 — read in full). Both done; this file is the
bridge document. Companion results: `t1540_movie_returns_vs_relations_deathstar_S21.out`.

## 1. What the paper proves, in our language

Feige's conjecture for even k ≥ 4, sharp (no polylogs): m ≥ 64n(n/ρ)^{k/2−1} hyperedges
force an **even cover** — a nonzero vector of the F₂ kernel of the vertex×edge incidence
system — of size ≤ 4kρ log n. Equivalently: dense k-sparse F₂ systems have short linear
dependencies (their own framing: the LDPC rate–distance tradeoff).

The proof walks on the **Kikuchi graph** K_ℓ(H): vertices = ℓ-subsets of [n], S ~ T iff
SΔT ∈ H, edges colored by SΔT. In repo terms: **a fixed-weight slice of the hypercube with
hyperedge-XOR moves** — our object family (the waggly d-strata of Q_m; THM-584's antipodal/
level-parity geometry; the Krawtchouk framework OPEN-Q-040 is precisely the spectral theory
of such slices). Their **paired vs unpaired cycles** (color-XOR zero vs not) = our **silent
vs expressive mutations**. Their **palette map** pal(S) ∈ F₂^H (odd-visited colors on a
walk from a root; well-defined below the unpaired girth; determines S given the root) = our
**tiling coordinates** (tournament = base path + tile-XOR vector), anchored at a root. The
proof's engine is new to us — **Lemma 2.3, the diagonal polynomial system bound**
(Frankl–Wilson/one-inclusion-graph species): if vertices carry distinct F₂^m labels with
adjacent labels at Hamming distance 1, and each vertex has a degree-≤d polynomial with
f_v(p_w) = 1{v=w}, then avg deg ≤ 2d. Their polynomials f_S = ∏_{v∈S}(1{v∈S₀} + Σ_{E∋v}X_E)
work because **on a fixed-weight slice, inclusion is equality** — the same slice trick our
merged metagraph uses.

(Sociological note, recorded without comment beyond this sentence: the acknowledgments
state the core innovation was found by GPT-5.6 Sol, with Claude Opus 4.8 and Claude Fable 5
also used — the honest-AI-science norms this repo runs on are becoming the field's.)

## 2. The T1540 probe: the wall movie's cut⊕cycle structure (RUN, exact)

The miss-pattern movie = the wall-crossing event stream: states σ(x) ∈ Z₇⁶ between events,
colors = walls (runner, index). Probe on four cores (exact, full period):

| core | events | states | cycle rank | returns | silent | min-rel L1 (count) |
|---|---|---|---|---|---|---|
| consecutive [50..55] (coherent) | 2114 | 1113 | **1002 (47%)** | 1001 | 0 | 4 (×96) |
| generic c=50 (incoherent) | 5432 | 5215 | **218 (4%)** | 217 | 0 | 4 (×4) |
| planted 48+96=144 | 4746 | 1582 | **3165 (67%)** | 3164 | 0 | 3 (×20) |
| far bank [0..5,50] | 378 | 336 | 43 (11%) | 42 | 0 | 3 (×40) |

**Findings.**
1. **The transference duality holds quantitatively**: recurrence (cycle-rank fraction)
   tracks the additive-relation spectrum — coherent cores recur massively, incoherent
   barely; a single planted relation (48+96=144, bringing the gcd-48 subcluster) cascades
   to 67% recurrence. Cut side (states/returns = simultaneous approximation) and cycle side
   (palettes/relations = the dual lattice) are Khintchine-transference duals, now measured.
   This is the movie-side face of THM-890's relation-spectrum error law.
2. **A small clean lemma (proved by the monotone-flow observation): the wall movie has NO
   silent cycles.** Between two visits of the same state within one period, each wall is
   crossed at most once (the flow is monotone; a wall's crossings = integer points in an
   interval of length < its period), so every crossed wall is odd-crossed: every return is
   expressive. In the paper's language: **for wall movies, unpaired girth = girth** — the
   generic-position case where their central dichotomy trivializes. Our movies are
   intrinsically "all-signal".
3. The movie's palette vectors form a code over F₂^W (W = wall set); its minimum distance
   = the shortest expressive cycle = a recurrence-quality invariant. Named lead below.

## 3. The Moore-bound LRC face: one honest negative, one real target

**Negative (measured, Probe B):** "covering-saturation forces shorter additive relations"
is VACUOUS at 13-element scale — covering and non-covering 13-sets in [1,200] have
identical shortest-relation statistics (mean L1 ≈ 3.1 both; 100% have L1 ≤ 4 both; n = 28
vs 280). Thirteen integers always carry tiny relations; the Moore forcing principle needs
hypergraph density m ~ n^{k/2}, which 13 elements never reach. Do not chase the naive
version.

**The real target (named):** the forcing principle applies at the **wall-system scale**:
the movie's wall hypergraph has n = 7Σf vertices and the palette code has m ~ recurrence
structure; the Moore bound then reads "recurrence density ⟹ short expressive cycles" —
i.e. a LOWER bound mechanism for how incoherent a cluster can be given its wall-code
density. Combined with finding 2 (no silent cycles), the movie's palette code is exactly
the object the paper's machinery speaks to.

## 4. Extension proposals (concrete, per repo thread)

1. **The diagonal-polynomial tool for metagraph degree bounds.** Our wiggly layers satisfy
   the Hamming-1 palette condition on the nose (tile coordinates). If iso-classes (or level
   sets of an invariant) admit degree-d diagonal polynomial systems in the tile variables,
   Lemma 2.3 bounds their internal wiggly degree by 2d — a NEW route to the unproved
   level-width/degree bounds of the metagraph atlas (klein's census rows have no theory;
   this is a candidate theory). First test: n = 5, label classes by tile-XOR palettes from
   a root class, find low-degree indicator systems (the H-polynomial? the parity
   invariants?).
2. **Krawtchouk ↔ Kikuchi**: our OPEN-Q-040 framework should absorb the Kikuchi spectral
   picture (Johnson-scheme slice of our hypercube analysis); conversely our exact
   metagraph quotients are S_n-symmetrized Kikuchi graphs — the repo has EXACT small-n
   spectra for objects their asymptotic method approximates. Potential export: exact
   unpaired-girth computations on symmetrized Kikuchi graphs.
3. **LDPC/engineering (mandate):** cite 2607.14068 in the circulant-LDPC deliverable — it
   is now THE sharp rate–distance obstruction for k-sparse parity systems (k even). Our
   `mod_rank_library` gains a use case: even-cover search = sparse F₂ kernel with weight
   minimization; the Kikuchi-BFS-with-palettes is an implementable algorithm for it.
4. **The [72,36,16] gauge thread (OPEN-Q-061/063):** the Moore bound constrains which
   sparse gauge presentations could exist — worth one session to check whether the k-sparse
   presentation densities in that thread sit above or below the forcing threshold.
5. **The movie palette code** (from §2.3): compute its parameters [#walls, rank, min
   distance] across the core taxonomy; conjecture: min distance = the coherence κ up to
   the section-resolution factor — this would make THM-889/890's coherence meter a COD ING-
   THEORY invariant, unifying the seam's classifier with the LDPC face.

-> HYP-7027, THM-889/890, T1540, THM-584, OPEN-Q-040/061/063, LEM-020,
engineering-synthesis-2026-03-10-S53; death-star-S21.
