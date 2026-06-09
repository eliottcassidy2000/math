# Erdős Problem 592 — survey, repo reframes, and the finite tree-grid core

**Author:** mac-mini-2026-06-09-S1
**Status:** working draft (results sections updated as the lab runs)
**Companion:** T768, HYP-2338..2340, THM-449,
`04-computation/erdos592_shuffle_pattern_lab_macmini_s1.py`,
`04-computation/erdos592_pattern_killers_macmini_s1.py`,
`04-computation/erdos592_treegrid_dichotomy_macmini_s1.py`

---

## 1. The problem, exactly

**Erdős Problem 592** (erdosproblems.com/592; Erdős 1987, $1000):
> Determine which countable ordinals β have the property that, if α = ω^β, then in any
> red/blue colouring of the edges of K_α there is either a red K_α or a blue K_3.

In partition-calculus notation: characterise the countable α with **α → (α, 3)²**
(red set of FULL order type α, or a blue triangle). Asked in general form (α, m) by
Erdős–Rado 1956; the (α,3) case re-posed with the prize in 1987.

### Verified state of the art

(Cross-checked against erdosproblems.com/592 + Džamonja–Koutsoukou-Argyraki–Paulson,
arXiv:2011.13218, which we read in PDF; their Theorem 2.2 and Observation 2.2.)

| α | α → (α,3)² ? | who |
|---|---|---|
| α not a power of ω | **NO** (easy split) | folklore/Galvin (Obs. 2.2 in 2011.13218) |
| ω, ω² | **YES** (indeed → (·,m) ∀m) | Ramsey; Specker 1957 |
| ω^n, 3 ≤ n < ω | **NO** | Specker 1957 |
| ω^β, β ≥ 3 decomposable | **NO** | Galvin (& Galvin–Larson 1974, pinning) |
| ω^(ω^γ), γ = ω^δ or ω^δ1+ω^δ2 (≤ 2 CNF summands) | **YES** | Schipperus 1999 thesis / APAL 2010 (γ=1: Chang 1972, Milner ∀m, Larson 1973 short proof; γ=2 also Darby independently) |
| ω^(ω^γ), γ with ≥ 4 CNF summands | **NO** | Schipperus (refutes the Galvin–Larson conjecture) |
| ω^(ω^γ), γ with exactly 3 CNF summands | **OPEN** | smallest open case: **α = ω^(ω³)** |

Related: the m-spectrum is finer than the (·,3) story — Larson, *An ordinal partition
avoiding pentagrams* (JSL 65 (2000)) and C. Darby, *Negative partition relations for
ordinals ω^(ω^α)* (JCTB) give negative relations against larger cliques for small
summand counts. (We could not access the full texts this session; the exact m-vs-summand
table is an import TODO — do not cite specifics beyond the titles.)

Erdős's general (α,m) problem and this problem are "very much open" as of the 2021
formalization survey; we found no 2022–2026 movement on the countable case (the
Raghavan–Todorčević "proof of a conjecture of Galvin" is the *uncountable/reals* Galvin
conjecture — different problem).

### Where the difficulty sits

Both directions are about ONE object class. Say a **witness** for α is a graph G on α with
1. G triangle-free, and
2. no independent set of G has order type α.

Then α → (α,3)² ⟺ **no witness exists** (take G = blue graph; red K_α = independent set
of full type). The open case asks: does ω^(ω³) carry a witness?

---

## 2. Reframes through this repo's objects

### 2.1 Triangle-free graphs with no full independent set (the witness frame)

The problem is literally about **triangle-free graphs on ordinals** — the repo's
triangle obsession lands on the other side of the board: THM-A (SC–NS ribs bipartite ⟹
triangle-free) is the *cheap* mechanism for triangle-freeness, and Galvin's negative for
decomposable α IS exactly that mechanism:

> α = α₁ + α₂ (both < α) ⟹ G = complete bipartite graph between the initial segment
> and the tail. Triangle-free needs... no wait, K_{β,γ} is triangle-free, every
> independent set lies in one side, sides have type < α. QED.

So **additive indecomposability of α = "no bipartite witness"**, and the entire
difficulty of 592 is about how much better than bipartite a triangle-free graph can do.
The known witnesses for ω^n (Specker) and ω^(ω^γ), γ ≥ 4 summands (Schipperus) are
"hierarchical multipartite-with-twists" constructions; the positive theorems say that
below 3 summands no twist works.

Z₂-grading view (repo dialect): a bipartite witness is a graph whose edges all cross a
2-coloring of α with both parts small. Indecomposable α admit no such grading; deeper
witnesses = deeper "iterated gradings" along the CNF — the levels of the grading tree
are exactly the CNF summand structure. The 2-vs-3-vs-4 summand boundary is a statement
about how many grading levels triangle-freeness can usefully exploit.

### 2.2 The tournament dictionary (THM-449 part 2)

A 2-coloring c of [α]² ⟺ a tournament T_c on α: for x < y, put x→y iff c{x,y} = red,
else y→x. Then:
* a red K_X (X of type α) = an **order-faithful transitive subtournament** of full type
  (T_c agrees with the ordinal order on X);
* a blue K_3 = an **order-reversed transitive triple** (z→y→x for x<y<z).

So 592 reads: *for which α does every tournament on α contain an order-faithful
transitive copy of α or an order-reversed transitive triple?* Rédei's theorem lives on
finite restrictions of exactly these objects (every finite restriction has a Hamiltonian
path = finite transitive "thread"); 592 asks for transfinite *order-aligned* transitivity.
This is the Erdős–Szekeres / monotone-substructure face of the problem.

### 2.3 Pairs are tiles: Larson's interaction schemes = staircase geometry (HYP-2338)

Larson's machinery (the proofs of ω² → (ω²,m), Specker, and ω^ω → (ω^ω,m), Milner;
formalized in Isabelle by Paulson) classifies pairs by **interaction schemes**:

* For ω²: U = {(a,b) : a<b<ω} ordered lexicographically — *this is the infinite
  staircase*: the repo's tile set {(x,y) : x ≥ y+2} is exactly such a pair set. Her four
  forms of pairs-of-pairs (form 0: a<b<c<d, form 1: a<c<b<d, form 2: a<c<d<b, form 3:
  a=c) are the repo's tile-pair classes DISJOINT / CROSS / NEST / SAME-VERTEX.
  **Verified mechanically** (lab part 5, 0 mismatches over [0,8)²-pairs): the tile-pair
  relation classifier and the pattern classifier agree class-by-class.
* For ω^ω: W = finite increasing sequences, pairs classified by forms 2k-1/2k — the
  number k of interleaved blocks, with the block-length partial sums c,d woven into the
  scheme. k = an **oscillation count** (Todorčević's osc, avant la lettre), and the
  layer-by-k filtration is the analogue of the repo's waggly Hamming-distance layers
  d=1..m. The partial-sum data c,d = the **cut-space/staircase-profile** data of the
  repo's cut ⊕ cycle split (CNF split level = cut data; within-level shuffle = cycle
  data).

### 2.4 The CNF ultrametric = the strip/level structure

For x ≠ y < ω^n, the **split level** ℓ(x,y) = largest exponent where the CNF
coefficients differ is an ultrametric valuation; pairs of ordinals carry a canonical
level filtration = the repo's strips. All known witness colorings are built from split
levels + within-level comparisons. (Specker positive at n=2 = "no witness can be built
from these"; see §3.)

---

## 3. The finite core (this session's computational results)

### 3.1 The pattern calculus (lab part 1–2)

Model ω^n as W(n) = increasing n-tuples, lex. The **pattern** of a pair = the
order/equality pattern of the 2n values. Any pattern set S ⊆ patterns defines a
"pattern-measurable" graph G_S. Triangle-freeness of G_S is a **finite, exact** check
(triple table over [0,3n)).

**Results, n=2** (6 patterns: DISJOINT, CROSS, NEST, SAME-LEFT, SAME-RIGHT, ADJACENT):
* 5 of 6 patterns have monochromatic triangles by themselves; the **unique** maximal
  triangle-free pattern set is **{ADJACENT}** — i.e. **the shift graph** Sh([ω]²)
  (edges (a,b)–(b,c)): the classical triangle-free, infinitely-chromatic graph appears
  canonically as *the only* triangle-free pattern graph on ω².
* {ADJACENT} is grid-avoidable: max S-free grid branching grows with the value universe
  (t ≈ 2,2,3,3,4 at V = 6,9,12,16,20 — a √V-ish profile, since dodging b=c burns values
  quadratically). **Specker's positive theorem in miniature**: no pattern witness, and
  the only triangle-free candidate is avoidable.

**Results, n=3** (31 patterns incl. shared-value ones; triple table = 2649 realizable
pattern-triples over [0,9)):
* exactly **13 maximal triangle-free pattern sets** (sizes 3–5);
* **every one of them is grid-avoidable** — exact backtracking finds S-free 3-grids by
  V = 14–24 for all 13 (early "killers" at V ≤ 14 were finite-size artifacts).

**Consequence (computational, V ≤ 24):** *Specker's witness for ω³ ↛ (ω³,3) is NOT
measurable in the pure order-pattern algebra of pairs.* Any witness must use structure
beyond the order pattern of the six values — value/length arithmetic (Larson's scheme
data: partial sums, i.e. staircase profiles) or deeper. This sharpens HYP-2339: the
"patterns suffice" intuition is FALSE at the very first negative level. The
uniformization that powers the positive proofs (Ramsey/Nash-Williams over forms) does
not capture the negative witnesses.

### 3.2 The tree-grid dichotomy Q(n,t) (lab part 3)

Strip away values entirely. On the leaves [t]^n of the complete t-ary tree of height n,
say G is a **finite witness** if G is triangle-free and every embedded complete binary
subtree (binary subgrid, 2^n leaves) contains an edge.

* Fact (used as the bridge): X ⊆ ω^n has order type ω^n ⟺ X contains a full
  ω-branching grid. So an infinite witness restricts to finite witnesses on finite
  grids; UNSAT of Q(n,t) at any t is finite, certified evidence for the positive
  direction at ω^n.
* Calibration n=1: Q(1,t) ⟺ complete triangle-free ⟺ t ≤ 2. ✓ (the Ramsey shadow —
  verified by exact DPLL).
* n=2 (ω² positive expected to shadow): witnesses exist at t = 3 (8 edges) and t = 4
  (32 edges, greedy); t = 5: [RESULT PENDING — exact DPLL running].
* n=3 (ω³ negative): [RESULT PENDING].

The n=2 witnesses die or survive by a Zarankiewicz-type mechanism: with all rows
independent, hitting all 2×2 rectangles across a row pair forces the bipartite
complement to be K_{2,2}-free (≤ z(t;2,2) ≍ t^{3/2} non-edges), so cross graphs are
dense, and dense 3-partite triangle-free is impossible — within-row edges only shift
the load. A finite Q(2,t) cutoff would be a clean, new, small Ramsey-type number
sitting directly under Specker's theorem.

---

## 4. Honest scoping

* Nothing here touches the open case ω^(ω³) yet; the session built the *calculus* at
  the ω^n levels where the answers are known, to expose the mechanism (pattern algebras
  too weak; witness structure = the thing to extract from SAT solutions).
* The "grid ⟺ full type" characterization is proved for the models used here (ordered
  sums argument, induction on n); the finite-shadow direction (UNSAT ⟹ evidence for
  positive) is a necessary-condition statement only — the infinite positive theorems do
  NOT follow from finite UNSAT (and finite SAT does not give infinite witnesses).
  Labels kept explicit in the lab outputs.
* Literature gaps to import: exact Darby/Larson m-spectrum statements; Schipperus's
  ≥4-summand witness construction (APAL 2010, paywalled); EHMR's presentation of
  Specker's ω³ witness. The lab is designed to *rediscover* the latter's structure
  independently (greedy/DPLL witnesses + structure analysis).

## 5. Next moves (handoff if session ends)

1. Finish Q(2,t) exact cutoff; push Q(3,t) witnesses to t=3,4 and CLASSIFY their edge
   structure (meet-level histograms, within-level patterns) — read Specker's coloring
   off the SAT solutions.
2. Lift the calculus one floor: ω^(ω^k) miniature (elements = finite sets of k-tuples;
   grids = trees-of-trees). Recompute the dichotomy; the 2/3/4-summand boundary should
   appear as a SAT/UNSAT phase line in (k = summand count, structural budget).
3. The join algebra: patterns of (values ∪ partial sums) — test whether Specker's
   witness becomes pattern-measurable there (Larson's schemes suggest yes); if so, the
   "scheme algebra" is the right finite presentation for the open case.
4. Oscillation parity: check whether the SAT witnesses' blue edges concentrate on a
   single oscillation parity class (repo GF(2) instinct; HYP-2339 refinement).
