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
- **The 2-plateau (EMPIRICAL + exhaustive small):** every natural design —
  paths, cycles, ladders, rook-rail grids, wildcard skeletons, both-end
  pincers — lands on score EXACTLY 2, and exhaustive k=1 search (n ≤ 5,
  pool {(1,0),(0,1),(1,1)}, all ≤2-wildcard placements, greedy-pruned seeds)
  finds nothing below 2. Everything sits on the line m + r = 2(n − t).
  **CONJECTURE (paths):** strict k=1 admits no score < 2. The mechanism is a
  frontier law: a straight frontier hands each new vertex one backward
  direction; all cancellation routes for a second direction exit forward.

## 4. Where < 2 must live, if anywhere

The three structural loopholes the plateau argument does not close:

1. **Label-sharing economies (k ≥ 2 deep products):** an axis-i label is one
   choice fanning into d_{i+1}···d_k parallel edges — the constructible
   fingerprint; iterated-doubling shapes [2]^k are exactly the KT proof
   skeletons. Search class D.
2. **Junction fabrics:** middle rows with two shared rung species make
   4-direction junctions; staircase cancellations could route second
   directions backward. Search class C.
3. **Symmetric-branch degeneration (the owner's μ₃ fragment):** a depth-2
   3-branch recursion whose μ₃-fixed branch carries a coincidence counts
   3+3+1 = 7 not 9 — an identification that shrinks the denominator's
   effective n. The verifiable format has no identifications, but T-seeding
   and shared labels can emulate; see
   `05-knowledge/results/FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md`.

If the plateau is a theorem for ALL constructible graphs, then the strict
reading is unwinnable below 2 and the benchmark's own soundness/encoding
claim (KT proofs ⟹ certificates at 1.75, 1.675…) fails for a THIRD reading
we have not reconstructed — most likely the intended semantics includes the
intuitive version's IDENTIFICATIONS (gluing copies along vertex sets), which
strictly enriches the game beyond both implemented modes. Next session
should implement identification-gluing as mode three: the intuitive
definition allows "identify x's copies in H_{i−1} and H_i", which merges
vertices and lowers n without touching m, r — exactly the 9→7 move of the
μ₃ fragment. **This is the leading hypothesis for where 1.75 and 1.675
certificates actually live.**

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
- The identification-gluing mode is unimplemented; until it exists, "strict
  min = 2?" remains a question about OUR two modes, not the benchmark.
