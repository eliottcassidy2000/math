# THM-453: Erdős 592 — the triangle-free reframe, the tile dictionary, and the finite tree-grid bridge

**Status:** PARTIAL — parts A–D are PROVED (proofs written out below); part E is exact
finite computation (verified, scripts + outputs in repo); part F is a proof sketch at
the infinite level for the translation-invariant case (labelled as such). The ordinal
open problem itself (ω^(ω³)) is untouched.
**Source:** mac-mini-2026-06-09-S1 (T768, HYP-2344..2346)
**Companion code:** `04-computation/erdos592_shuffle_pattern_lab_macmini_s1.py`,
`erdos592_pattern_killers_macmini_s1.py`, `erdos592_treegrid_dichotomy_macmini_s1.py`,
`erdos592_treegrid_pysat_macmini_s1.py` (+ .out files in 05-knowledge/results/)

Erdős Problem 592 ($1000, OPEN): characterise countable α with α → (α,3)².

**Witness frame.** A *witness for α* is a triangle-free graph G on α such that no
independent set of G has order type α. Then α → (α,3)² ⟺ no witness for α exists
(G = the blue graph; a red K_α = an independent set of full type).

---

## A. Tournament dictionary (PROVED)

Identify a 2-coloring c of [α]² with the tournament T_c on α: for x < y, x→y iff
c{x,y} = red, else y→x. Red-homogeneous X of type α = a set on which T_c is transitive
and agrees with the ordinal order; blue triangle = order-reversed transitive triple. So:

> α → (α,3)² ⟺ every tournament on α has an order-faithful transitive subtournament
> on a set of type α, or an order-reversed transitive triple.

Proof: unfold both dictionaries; the correspondences are bijective. ∎

## B. Decomposable ordinals have bipartite witnesses (Galvin's observation, with proof)

If α = α₁ + α₂, α₁, α₂ < α: let G be the complete bipartite graph between the initial
segment I (type α₁) and the final segment F (type α₂). Bipartite ⟹ triangle-free; an
independent set lies inside I or inside F, so has type ≤ max(α₁,α₂) < α. ∎
(Repo dialect: a Z₂-grading witness. "Additively indecomposable" = "no bipartite
witness"; all known deeper witnesses are iterated gradings along the CNF.)

## C. Grid characterization of full-type subsets (PROVED)

Model ω^n as N^n, lex order. A **full grid** = an ω-branching height-n tree of value
choices; its leaf set ⊆ N^n. A **binary subgrid** picks 2 children at every node
(2^n leaves). [Same statements hold in the increasing-tuple model W(n).]

**Lemma C.** X ⊆ N^n has order type ω^n ⟺ X contains the leaf set of a full grid.

Proof. (⟸) By induction the leaf set of a full grid has type ω^{n-1}·ω = ω^n, and
every subset of N^n has type ≤ ω^n. (⟹) Induction on n. X = Σ_a X_a (ordered sum over
first coordinates), otp(X_a) ≤ ω^{n-1}. If only finitely many sections had type
exactly ω^{n-1}, then beyond a bound every section has type < ω^{n-1}, i.e.
≤ ω^{n-2}·m_a (m_a < ω), and an ω-sum of such is ≤ ω^{n-1} — total otp(X) < ω^n,
contradiction. So infinitely many sections have full type; choose a₀ < a₁ < ⋯ among
them, recurse inside each, assemble the grid. ∎

## D. The finite tree-grid bridge (PROVED — compactness; corrects MISTAKE-066)

On the ABSTRACT complete t-ary tree of height n with leaves [t]^n define:

> **Q(n,t):** ∃ triangle-free graph on [t]^n meeting every binary subgrid
> (equivalently: with no independent binary subgrid).
> **R(n,2)** := least T with Q(n,T) false (∞ if none).

Q(n,t) is antitone in t: a (t+1)-grid witness restricted to an embedded t-subgrid is a
t-grid witness, because the hitting edge of each binary subgrid lies among that
subgrid's own leaves. So Q(n,t) ⟺ t < R(n,2).

**Theorem D1.** If Q(n,t) holds for all t < ω then ω^n ↛ (ω^n,3)², witnessed by a
*strong witness*: triangle-free with no independent binary subgrid of the full grid.

Proof. Restrictions are coherent (above), each level set {witnesses on [t]^n} is
finite and nonempty, so König's lemma gives a coherent sequence G_t; let
G = ∪_t G_t on N^n (leaves of the full grid; [t]^n ↑ N^n using children 0..t−1).
G is triangle-free (any triangle lives in some [t]^n) and meets every binary subgrid
of the full grid (it lies in some [T]^n, and G_T ⊆ G meets it — the meeting edge is
inside the subgrid). If X ⊆ N^n were independent with otp(X) = ω^n, Lemma C gives a
full grid inside X, hence a binary subgrid inside X, hence an edge of G inside X —
contradiction. ∎

**Corollary D2.** ω^n → (ω^n,3)² ⟹ R(n,2) < ∞.

So Specker's ω² → (ω²,3)² FORCES a finite cutoff **R(2,2)**: a previously
uncontemplated (to our knowledge) finite Ramsey-type number sitting directly beneath
Specker's 1957 theorem. The cutoff exists by D2 + Specker, but no explicit bound falls
out of the compactness proof — quantifying it is a fresh finite combinatorics problem.

**Asymmetry note.** ω^n ↛ (ω^n,3)² does NOT formally imply Q(n,t) SAT for all t
(a witness may admit independent binary subgrids while killing all full-type sets).
"R(n,2) = ∞" = existence of *strong* witnesses, a priori stronger than the negative
relation. (MISTAKE-066 was asserting the false converse direction in an early
docstring; caught and corrected the same session.)

## E. Computed values and structure (exact; Glucose3 + CEGAR, and exhaustive DPLL)

* **R(1,2) = 3** (calibration: "triangle-free + meet every pair" dies at 3 vertices —
  the Ramsey ω → (ω,3) shadow).
* **R(2,2) > 10** (witnesses found at every t ≤ 10; run continuing): despite Specker
  FORCING R(2,2) < ∞, finite witnesses persist with edge counts
  1, 6, 18, 40, 77, 143, 208, 313, 420 at t = 2..10. The hidden finite threshold is
  large — the Mantel/Zarankiewicz pressure (see F) has constant-factor slack ≈ 2.
* **Q(3,t) SAT for t ≤ 5** (witness edge counts 1, 18, 90, 223), consistent with ω³'s
  negative relation having strong witnesses.
* Witness structure: at n=2 edges concentrate on cross classes (meet level 0) with
  second-coordinate increase, spread over ALL root-gaps with declining frequency; at
  n=3 dominated by meet-level-0, suffix-profile ('<','<') (domination-type) plus
  meet-level-1 ('<') classes, again gap-graded. The witnesses are NOT order-pattern
  functions: gap magnitudes are used essentially (cf. part E.3 below).

### E.1 Patterns = tile-pair geometry (HYP-2344, VERIFIED)
The order/equality patterns of pairs of increasing 2-tuples are exactly the repo's
staircase tile-pair classes {DISJOINT, CROSS, NEST, SAME-LEFT, SAME-RIGHT, ADJACENT};
mechanical check over all pairs in [0,8)²: 0 mismatches. Larson's interaction forms
0,1,2,3 (her ω²/Specker-positive proof) = {DISJOINT, CROSS, NEST, SAME-LEFT}.

### E.2 n=2 pattern algebra (exact)
Of the 6 patterns, 5 admit monochromatic triangles; the unique maximal triangle-free
pattern set is **{ADJACENT} — the shift graph** Sh([ω]²), and it is grid-avoidable
(S-free t-grids exist with t ≈ √V). Hence: no pattern-measurable witness at ω², with
the canonical triangle-free high-chromatic graph appearing as the lone candidate —
Specker positive in miniature.

### E.3 n=3 pattern algebra (exact)
31 patterns, 2649 realizable pattern-triples, exactly **13 maximal triangle-free
pattern sets** (sizes 3–5) — and ALL 13 are grid-avoidable (exact backtracking;
S-free 3-grids by V ≤ 24 in every case). Consequence (computational, V ≤ 24):
**Specker's ω³ witness is not measurable in the pure order-pattern algebra of pairs**;
it must use value/length arithmetic (gap magnitudes / Larson-scheme partial sums =
staircase profiles). This corrects the naive form of HYP-2345.

## F. The graded-relation (additive-ladder) reformulation — invariant witnesses

Write ω² = N² lex; classify a graph's edges by row structure: within-row graphs R_a on
columns, and cross-row bipartite relations B_{a,a'} ⊆ N×N. Call a witness
**translation-invariant** if R_a ≡ R and B_{a,a'} = B_{a'−a} (gap-graded).

**Proposition F1 (PROVED, bookkeeping).** A translation-invariant strong witness for
ω² is exactly a family (R, {B_g}_{g≥1}) with:
 (i) R triangle-free;
 (ii) ∀g₁,g₂: B_{g₁} ∘ B_{g₂} ∩ B_{g₁+g₂} = ∅ (no cross-cross-cross triangles —
      composition-freeness graded by (N,+));
 (iii) forward and backward B_g-neighborhoods are R-independent (no cross-cross-row
      triangles);
 (iv) ∀g: every 2×2 rectangle {b₁,b₂}×{c₁,c₂} with both column-pairs R-independent
      meets B_g (hitting all binary subgrids).
Proof: direct translation of triangle-freeness and subgrid-hitting through the
row/gap classification; every triangle uses gaps with g₁+g₂ = g₃ (or 2,1 or 0 cross
edges), every binary subgrid is a (root-pair, column-pairs) rectangle. ∎

**F2 (proof sketch, infinite level).** No such family exists — i.e. even before
quoting Specker, the invariant case dies by an additive/density argument:
(iv) forces each B_g to be rectangle-dense, so its bipartite complement is K_{2,2}-free
on the R-independent pair structure (Zarankiewicz: complements are O(t^{3/2})-sparse on
[t] blocks); then typical forward neighborhoods are co-small, any two co-small sets
meet, so B_{g₁} ∘ B_{g₂} contains all pairs with typical endpoints — but then (ii)
forces B_{g₁+g₂} inside an atypical (sparse) set, contradicting ITS rectangle-density.
Edge cases (vertices with small B-degree hiding behind partners, load carried by R)
are exactly what makes the finite threshold R(2,2) large; turning this sketch into an
explicit bound on R(2,2) is the sharpest finite handoff of the session.

**F3 (why n=3 breathes).** In ω³ = N³ the "columns" of the same decomposition carry
their own level structure: condition (iv) weakens from "hit every 2×2 rectangle" to
"hit every binary subgrid of the column-ω²" — a much sparser demand — so the B_g can
afford the composition-freeness (ii). This is the mechanism by which the dichotomy
flips between n=2 and n=3, made quantitative by the Q(n,t) computations.

**F4 (the ladder reading — the repo connection).** Condition (ii) is a *Schur/cauldron
condition* (avoid g₁+g₂ = g₃ relations) on a gap-graded family of relations; (iv) is a
density floor. Erdős 592's mechanism at level ω^n is thus "dense graded relation
families with no additive composition" — the same density-vs-additive-relation tension
as the repo's Sidon/B_h ladder (THM-446/HYP-2314: Sidon ⟺ C4-free ⟺ minimal additive
energy; cauldron = Schur 3-term). For ω^(ω^γ) the gap grading lives over the CNF
semigroup of γ; **the number of CNF summands of γ = the rank of the grading
semigroup**, and Schipperus's dichotomy (≤2 summands positive, ≥4 negative, =3 OPEN)
becomes: at which grading rank does density-vs-composition-freeness flip? The open
case ω^(ω³) is the rank-3 instance. (Framing — precise at n=2,3 via F1–F3; the
ω^(ω^γ) lift is the program of HYP-2346.)

## Known results recorded (literature, not ours)

Specker 1957: ω²→(ω²,m)², ω^n↛(ω^n,3)² for 3≤n<ω. Chang 1972/Milner: ω^ω→(ω^ω,m)²
(Larson 1973 short proof). Galvin–Larson 1974: β≥3 with the property ⟹ β=ω^γ.
Schipperus 1999/2010: ω^(ω^γ)→(ω^(ω^γ),3)² for γ a sum of ≤2 indecomposables
(γ=2 also Darby); ↛ for ≥4 summands (refutes the Galvin–Larson conjecture). OPEN:
exactly 3 summands; smallest open case α = ω^(ω³). m-spectrum refinements: Larson
"An ordinal partition avoiding pentagrams" (JSL 2000), Darby "Negative partition
relations for ordinals ω^(ω^α)" (JCTB) — exact statements not yet imported (paywalled).
(Sources verified this session: erdosproblems.com/592 via reader proxy;
arXiv:2011.13218 read in PDF, Thm 2.2 / Obs 2.2.)
