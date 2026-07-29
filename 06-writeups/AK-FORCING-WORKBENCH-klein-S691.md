# The arithmetic-Kakeya forcing-graph workbench (klein-S691, 2026-07-28)

**Scope:** in-repo workbench for the external benchmark problem "find an
X-constructible graph with a forcing pair of score ≤ 1.675" (arithmetic
Kakeya certificates). Engine: `04-computation/ak_forcing_engine.py`.
Exploit: `04-computation/ak_exploit_cert.py`. Search:
`04-computation/ak_strict_search.py`. Status labels below are literal.

## 0. The benchmark in repo language

AK(α) is a family of additive-combinatorics statements strong enough (α → 1)
to give Kakeya-dimension bounds via Bourgain's slice method; the human record
is Katz–Tao 2002's α = 1.67513… (the (1,2)-root of x³ − 4x + 2). The
benchmark encodes "proofs" as finite certificates: a grid graph
d₁×…×d_k with axis-labels f_i into X ⊂ Z² ((0,0) allowed; nonzero labels
avoid the line a+b=0), plus initial seeds R (singleton x·δ_v) and initial
wildcards T₀, with closure rules; success = every vertex acquires a
(1,−1)-functional; score = (m+r)/(n−t); the target ≤ 1.675 would beat the
24-year-old record. This is certificate culture — the same shape as our
B5/Bonferroni and deck/fan gates: a finite decidable object whose validity
implies an infinite statement (here CITED from the benchmark's framing, not
re-proved).

## 1. FINDING (verified): the literal spec is unsound — score → 1

Rule (1) as literally written constrains only the FIRST i coordinates of the
two endpoints (suffixes free). Differencing two loose edge-vectors sharing a
far endpoint yields within-layer transporters x(δ_{(p,s)} − δ_{(p,s')}),
unavailable in the glued-copies (equal-suffix) reading. Consequence,
machine-verified both ways (`ak_exploit_cert.py`):

- dims [3,3], X = {0,(1,0),(0,1),(1,1)}, m=10, r=4, t=0: **score 14/9 ≈ 1.556,
  forced under the literal reading; stuck at 3/9 under equal-suffix.**
- family [D,n₀]: score = 1 + 1/D + 1/n₀ − 1/(Dn₀) → 1. A 6×6 member scores
  47/36 ≈ 1.306.

Since certificates of score → 1 would "prove" the full arithmetic Kakeya
conjecture by a 36-vertex object, the literal reading cannot be the intended
one; any sound verifier must implement equal suffixes. Both modes are in the
engine; every claim below is about the STRICT (equal-suffix) game.
*Actionable: report upstream to the benchmark authors.*

## 2. The strict game is a species-flow / illegal-junction game

Reformulation (exact, implemented): lattice vectors are Z²-species flows on
the labelled graph; a vertex v ∉ T fires iff some flow supported on T ∪ {v}
deposits a nonzero (a,−a) at v.

- **Telescopes:** flow passes through a non-T vertex only with exact species
  cancellation ⟹ monochromatic runs transport one direction end-to-end at
  interior-blind cost.
- **Junctions:** turning at a non-T vertex needs ≥3 incident species with a
  Z-relation (e.g. (1,1) = (1,0) + (0,1)); a seed at a vertex doubles as a
  junction-enabler there.
- **Wildcards (T-vertices):** universal junctions.
- **Firing = illegal junction:** the deposited species (a,−a) is exactly the
  direction EXCLUDED from X (labels have a+b ≠ 0). A fired vertex is where
  an x-in/y-out turn happens with no balancing species — the game places
  illegal junctions, one per vertex.

## 3. Two proved floors and an empirical plateau

- **Slot bound (PROVED, elementary):** a firing flow's v-column combines only
  v-incident generators; no single X-direction spans (1,−1); so every fired
  vertex needs ≥2 incident generators (edge-endpoints/seeds) of independent
  directions: 2m + r ≥ 2(n − t).
- **Rank floor (PROVED):** let Φ(T) = dim_Q(L_Q + ⊕_{v∈T} Q²). A legal firing
  adds ≤ 1 to Φ (the (1,−1)-direction at v was already present mod T);
  Φ(T₀) ≤ (m+r) + 2t; counting fired vertices gives **m + r ≥ n − t, i.e.
  score ≥ 1** — the format cannot certify below AK(1), reassuringly.
- **The 2-plateau (REFUTED for k ≥ 2 — session headline):** every natural
  HAND design (paths, cycles, ladders, rook-rail grids, wildcard skeletons,
  both-end pincers) lands on exactly 2, and exhaustive k=1 search (n ≤ 5,
  pool {(1,0),(0,1),(1,1)}, all ≤2-wildcard placements, greedy-pruned seeds)
  finds nothing below 2 — the path conjecture stands. But randomized search
  over deeper products BROKE the plateau (all FINITE-EXACT, re-verified from
  scratch):
    * **score 13/7 ≈ 1.857**, dims [2,2,2], m=10 r=3 n=8 t=1, X ⊆
      {(1,0),(0,1),(1,1),(1,2),(2,1)}: fs = [{(1,)→(1,2)},
      {(1,1)→(1,0),(2,1)→(1,0)}, {(2,1,1)→(2,1),(2,2,1)→(0,1)}],
      T₀={(1,1,2)}, R={(1,2,2)→(0,1),(2,2,1)→(1,1),(1,1,1)→(1,0)}.
      Fires in 2 rounds; round 2 fires five vertices at once — genuinely
      collective lattice forcing, not sequential percolation.
    * score 15/8, dims [2,2,2], t=0; score 19/10, dims [4,3] junction strip;
      score 31/16, dims [2,2,2,2].
  So the strict game's floor is < 2 and label-shared deep products are where
  the descent lives; the proved floor remains 1 (rank bound). Next rungs:
  11/6, 9/5, 7/4 (the KT99 exponent — reaching it would validate the game's
  faithfulness to the human proof ladder), then 5/3 < 1.675.

- **Mechanism of 13/7 (extracted, exact):** the round-1 firing of (2,1,1) —
  no neighbor in T — is the flow (×2, depositing 2·(1,−1)):
  `−6·y@(1,2,2)seed + 6·w@(2,2,1)seed + 1·u[(1,1,2)−(2,1,2)] +
  3·u[(1,2,2)−(2,2,2)] + 3·x[(1,1,2)−(1,2,2)] + 6·x[(2,1,1)−(2,2,1)] −
  3·x[(2,1,2)−(2,2,2)] − 2·p[(2,1,1)−(2,1,2)] − 6·y[(2,2,1)−(2,2,2)]`.
  Reading: the fired vertex's own two species combine (6·(1,0) − 2·(2,1) =
  2·(1,−1)); everything else is a circulation through five non-T vertices
  with exact cancellation, enabled by FOREIGN seeds acting as junction fuel
  (a seed is a lattice generator — infinitely reusable at any coefficient)
  and closed through the wildcard. Design principles extracted:
  (i) seeds ≠ local slots only — each seed is a permanent junction-enabler
  for ALL later firings; (ii) pass-through requires ≥3 species or a
  T-vertex, so hub vertices (degree ≥ 3, mixed species) are the scarce
  resource; (iii) axis-shared rails serve simultaneously as slots and as
  return-circulation, which is why deep products beat every k ≤ 2 design.

## 4. Where < 2 must live, if anywhere

The three structural loopholes the plateau argument does not close:

1. **Label-sharing economies (k ≥ 2 deep products):** an axis-i label is one
   choice fanning into d_{i+1}···d_k parallel edges — the constructible
   fingerprint; iterated-doubling shapes [2]^k are exactly the KT proof
   skeletons. Search class D.
2. **Junction fabrics:** middle rows with two shared rung species make
   4-direction junctions; staircase cancellations could route second
   directions backward. Search class C.
3. **Identification-gluing, independently motivated:** a merged class
   inherits its members' incident species and can manufacture a cheap
   junction.  The smaller denominator alone is harmful; the possible gain
   is the larger numerator saving.  Measured on twin `[2,2,2]` witnesses,
   merge gives `12/7` versus `7/4` for an edge-plus-seed substitute.  This
   mechanism must stand on its own: the owner's `3+3+1=7` fragment has now
   been decoded as THM-2473's Keller fibre, where two sheets escape to
   infinity rather than two finite vertices being identified.  See
   `05-knowledge/results/FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md`.

If the plateau is a theorem for ALL constructible graphs, then the strict
reading is unwinnable below 2 and the benchmark's own soundness/encoding
claim (KT proofs ⟹ certificates at 1.75, 1.675…) fails for a THIRD reading
we have not reconstructed — most likely the intended semantics includes the
intuitive version's IDENTIFICATIONS (gluing copies along vertex sets), which
strictly enriches the strict game.  Mode3-v2 now implements a permissive
version of that operation and reaches `9/5`; its witness still needs the
same-H legality/compilation theorem described below.  Identification remains
the leading structural hypothesis for where `1.75` and `1.675` certificates
could live, but THM-2473's `9 -> 7` escape is no evidence for it.

## 4b. The format hierarchy (corrected reading, late-session)

Close re-reading of the intuitive definition shows the two benchmark
formats disagree in BOTH directions, and the session's games sit as:

    strict (= verifiable minus the loose bug)  ⊆  mode3-v2  ⊇  intuitive

- The INTUITIVE recursion takes k+1 copies of ONE H per step: per-(step,
  suffix) gluing choices INCLUDING identifications and per-suffix labels —
  but identical inner structure across copies (same H).
- The VERIFIABLE format allows f_i to depend on the whole prefix (inner
  labels varying per outer position — not directly intuitive-legal;
  plausibly emulable via deeper same-H recursion) but has NO merge device
  and suffix-uniform step labels ((0,0) ∈ X exists only so seed functions
  are total — not an identification marker).
- MODE3-V2 (this session's explorer) = verifiable + per-suffix labels +
  legal-slot merges: a common superset. Its records are upper bounds for
  the intuitive game modulo a same-H legality pass.

Current records: **13/7 strict** ([2,2,2], robust attractor, two
independent witnesses); **9/5 mode3-v2** ([3,3] with a 5-fold merged hub;
m=5, r=4, n=5, t=0 — merges manufacture free species-rich junctions,
collective bootstrap from single seeds, no wildcards). The 9/5 witness's
rows differ internally, so intuitive legality is OPEN; the compilation
questions (merges → verifiable; prefix-varying labels → same-H depth) are
the benchmark's unstated equivalence lemmas and the workbench's top
theory targets.

## 5. Session artifacts

- Engine + verifier (both rule-1 readings), exact rational closure, 6-line
  submission emit/parse, baselines (trivial=2, rail-grids=2).
- Machine-verified unsoundness exploit family for the literal reading
  (14/9 at n=9; → 1).
- Exhaustive small-path scan (min = 2) and randomized classes B/C/D runs
  (log in session scratch; results to be appended on completion).
- HYP-9046 / HYP-9047: external-import reframes for the LRC(14) wall and
  budget layers (independent of this workbench).

## 6. Honest boundaries

- No AK statement is proved or claimed here; scores are certificate
  arithmetic, and the certificate⟹AK theorem is the benchmark's, unaudited.
- The 2-plateau conjecture is open even for paths (exhaustive only to n=5,
  greedy seeds only).
- Mode3-v2 implements permissive legal-slot merging, but its `9/5` record is
  not yet compiled into the intuitive same-H recursion.  The strict and
  intuitive formats remain incomparable until the missing compilation
  lemmas are proved.
