# THM-471: the B3 general-shape enumerator — first Schipperus-forced cutoffs, the XOR quotient, and the m=3 (open-case) probe

**Status:** parts A–C PROVED (tuple grammar = faithful miniaturization; finder
completeness; XOR-quotient linearity); part D COMPUTED (statuses inline; UNSAT
verdicts sound independently of finder completeness — see C2). Implements
POKE Task 1.2 (THM-460 line 150's "general-shape recursive enumerator").
**Source:** kind-pasteur-2026-06-11-S1 (HYP-2394/2395). Scripts
`04-computation/erdos592_shape_miniatures_kp0611.py`,
`04-computation/erdos592_shape_xor_quotient_kp0611.py` (+ .out files).

## A. The tuple grammar (PROVED faithful to THM-460 B2/B3)

Exponents are tuples δ ∈ N^m (δ[i] = CNF coefficient of ω^i). Ambient
A(m; s, c) = functions [s]^m → [c]; positions sorted most-significant-first by
reversed-tuple lex, so value-tuple lex order = largest-disagreement order and
list index order = ambient order (the m=1 invariants of THM-460 D survive).

  Bin(0) = point.
  FINITE peel (least nonzero index 0): two Bin(δ−e₀) halves A < B, ONE common
    cross-split position for all cross pairs, internal splits strictly less
    significant.
  LIMIT peel (least nonzero index i ≥ 1): M order-separated parts
    Bin((δ−e_i) + j·1 on positions 0..i−1), j marching.
  BT(m,M) = order-separated parts Bin((j,…,j)), j marching.

Faithfulness: ω^δ = ω^(δ−e_{i₀}) · ω^(ω^{i₀}) by CNF peeling; B2 at rank i₀
turns the index structure into a tower whose part exponents embed in positions
0..i₀−1 — exactly the limit production. The finite production is THM-453's
ordered-sum pruning. March conventions: **j1** (j = 1..M; faithful to THM-460
D's m=1 games) and **j0** (j = 0..M−1; the truncated march, with the j=0 part a
single point). At m=1 the checker coincides with is_binary_grid (verified: 0
mismatches over all 4-subsets of [3]³).

**Size analysis (the miniature-design table; m≥2 sizes compound through limit
peels):** |BT| at (m,M): m=2: j1 = 4 (M=1), 156 (M=2); j0 = 7 (M=2). m=3:
j1 = 16 (M=1), **3,506,256 (M=2 — vacuous in every feasible ambient**, the m≥3
analogue of THM-460 D's vacuousness guard — this is WHY the probe ladder must
use M=1-j1 and the j0 truncation); j0 = 43 (M=2).

## B. Finder architecture and completeness (PROVED)

ShapeFinder: recursive generator (propose halves/chunks with order separation +
independence threading via bitmasks; dispose with is_bin_shape) — the
TowerFinder architecture of THM-460 D, shape-generic.

**C1 (completeness).** Any independent Bin(δ) embedding, sorted in ambient
order, decomposes uniquely: its fin-halves / lim-chunks occupy consecutive index
ranges of grammar-determined sizes and are themselves sorted Bin embeddings of
the sub-exponents. gen() enumerates all candidate tuples with exactly these
order constraints and accepts precisely those passing is_bin_shape. ∎
(Empirics: cross-validation vs brute-force combination enumeration on random
graphs, 0 disagreements in 100 trials over five m ≤ 2 ambients — m=1 M=2 j1 at
(2,3),(3,2); m=2 M=1 j1 at (2,2),(2,3); m=2 M=2 j0 at (2,2); the m=3 structural
case has no feasible brute partner, so its completeness rests on C1's proof.
The finder also carries a constructive cross-split pruning — enumerate the
cross-split position p above A's internal splits, restrict the B-half to points
agreeing with A[0] above p and differing at p; complete since a sorted Bin's
halves agree per-half at p and above — and a node budget that converts search
explosions into honest TIMEOUTs, never SAT certificates.)

**C2 (soundness asymmetry — MISTAKE-067 discipline).** UNSAT verdicts depend
only on CLAUSE soundness (every blocked set passes is_bin_shape; triangle
clauses definitional), NOT on finder completeness. Completeness gates only SAT
verdicts ("no independent shape remains"). All m=1 SAT cells additionally match
THM-460 D's independent results.

## C. The XOR quotient at c=2 (PROVED + the THM-465 shortcut ported)

At c=2 the ambient is {0,1}^(s^m) and the disagreement set D(x,y) ⊆ positions
is an XOR-invariant feature with **F₂-linear triangle composition**:
D(x,z) = D(x,y) Δ D(y,z). Hence an XOR-measurable rule T is triangle-free ⟺
**T is a sum-free set in the group F₂^(s^m)** — THM-469's sum-free seam
reappearing one level up, now in an elementary-abelian 2-group (cap-set/Schur
world). Every group triple {d₁, d₂, d₁⊕d₂} is realizable (x = 0, y = 1_{d₁},
z = 1_{d₁⊕d₂}-shifted), so the triangle clause set is purely group-theoretic;
hitting clauses come from CEGAR against the complete ShapeFinder.

## D. Computed cells (Glucose3 + CEGAR; UNSAT sound by C2; SAT cells re-verified)

Calibrations (all match THM-460 D / THM-453):
* m=1 M=1: SAT (1,2); UNSAT (1,3), (2,2) — the R(1,2)=3 Ramsey shadow.
* m=1 M=2 j1: SAT (2,3) 5 edges, (3,2) 5 edges, (2,4) 31 edges.
* XOR consistency: m=2 M=1 (2,2) is raw-UNSAT and xor-feature-UNSAT (no clash).

**New numbers (first of their kind):**
* **m=2 M=1 j1: UNSAT at (2,2) (N=16, 2.2s) and (2,3) (N=81, 9611 shape clauses,
  23s)** — the first finite cutoff data points forced by Schipperus's positive
  ω^(ω²) theorem (THM-460 C2 predicted their existence; these are the first ever
  computed). Sweep at (2,4), (2,5), (3,2), (2,6): see .out (in flight).
* **m=2 M=2 j0: SAT at (2,2)** via a 5-class XOR table (40 edges, re-verified) —
  after M=1 dies, the M-ladder (longer towers = weaker demands) reopens, so the
  (M; s, c) cube has genuine structure at m=2.
* **m=3 probes at (2,2)** (M=1 j1, 16-leaf shape, N=256; M=2 j0, 43-leaf shape):
  raw + XOR-quotient runs in flight at write-up — outputs stream to
  `05-knowledge/results/erdos592_shape_{miniatures,xor_quotient}_kp0611.out`.
  These are the first finite data of any kind on the ω^(ω³) open case's miniature.

## Honesty

- The j0 march is a TRUNCATION design choice (the infinite tower's cofinality
  only constrains the tail); j0 results are labeled as such and not mixed with
  j1 rows.
- m=2 M=1 cutoffs are cutoffs of THIS miniature family; C2-forcing guarantees
  cutoffs exist for every M, but says nothing about WHERE, and the m=2/m=3
  differential is evidence-generating only (THM-460's honesty note inherited).
- m=3 M=1 UNSAT (if it lands) would be consistent with both signs of the open
  case; only persistent SAT across growing (s,c) is asymmetric evidence.

**Cross-refs:** THM-460 (grammar source, C1/C2 bridges), THM-469 (sum-free
seam; the XOR quotient is its group version), THM-465 (feature-quotient method),
HYP-2394/2395.
