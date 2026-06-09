# The additive ladder reaches the ordinals

**Source:** mac-mini-2026-06-09-S1 (Erdős 592 session; T768, THM-453, HYP-2344..2346)

## The arc

Erdős Problem 592 asks which countable α satisfy α → (α,3)². It looks like pure
infinitary combinatorics — the answer known so far is indexed by how many Cantor
normal form summands an exponent has (≤2 yes, ≥4 no, =3 open). Nothing about it
advertises a connection to additive combinatorics.

But push the witness objects (triangle-free graphs on α with no full-type independent
set) through a finite shadow — Q(n,t): triangle-free graphs on the t-ary grid of
height n with no independent binary subgrid — and ask what a *translation-invariant*
witness is. The answer (THM-453 F1) is a family of bipartite relations B_g graded by
the row-gap g, satisfying

    B_{g₁} ∘ B_{g₂}  ∩  B_{g₁+g₂}  =  ∅        (triangle-freeness)

against a Zarankiewicz-type density floor (every 2×2 rectangle of independent column
pairs must be hit). The first condition is a **Schur condition**: the family must
avoid the additive relation g₁ + g₂ = g₃ *at the level of relation composition*. The
second says each grade must be dense.

This is the repo's additive-relation ladder again (THM-446/HYP-2314, the Erdős-64
session): Sidon sets avoid 4-term relations a+b=c+d and are exactly the C₄-free
extremal objects; the cauldron game's "boils" are Schur 3-term relations a+b=c;
density-versus-relation-avoidance is the whole game. Here the same tension reappears
one category up — not sets of integers avoiding additive relations, but *graded
families of relations* whose composition must avoid the grading's addition. A triangle
in the ordinal graph IS an additive relation among gaps, decorated with column data.

## Why the dichotomy can flip with dimension

At n=2 (ω²) the columns are structureless, so the density floor is "hit every 2×2
rectangle" — brutally strong; dense relations compose to nearly-everything, and the
Schur condition strangles the family. Outcome (computed exactly this session, after
fixing the verifier): **R(2,2) = 5** — by t=5 no triangle-free graph survives the
binary-subgrid demand. Specker's 1957 positive theorem ω² → (ω²,3)² *forces* this
finite death (THM-453 D2), but the compactness proof gives no bound; the actual value
is 5, sitting next to R(1,2)=3 (the Ramsey shadow).

At n=3 (ω³) the columns are themselves a level deep, the density floor weakens to
"hit the column-side binary subgrids", and the graded family can afford
composition-freeness: Q(3,4) is SAT (346 edges, complete verification) where Q(2,5)
is UNSAT. The dimension does not add power to the witness; it *weakens the density
demand* — that asymmetry is the entire Specker dichotomy in finite miniature.

Climbing to ω^(ω^γ), the gaps live in ω^γ and the grading semigroup acquires the CNF
structure of γ; the number of summands of γ is the **rank** of the grading. So
Schipperus's theorem ("≤2 summands positive, ≥4 negative") reads: the
density-vs-Schur tension kills rank-1 and rank-2 graded families and is survivable at
rank ≥4. Erdős 592's open case is precisely **rank 3**. The problem that looked like
transfinite exotica is a threshold question in a graded additive-relation algebra —
the same species of question as "is there a dyadically-Sidon min-degree-3 graph"
(Erdős 64) one floor up the repo's ladder.

## The smaller echoes

* The unique triangle-free *pattern* graph on ω² is the **shift graph** — the
  canonical triangle-free high-chromatic object appears with no design choices, and
  Specker-positivity at ω² is exactly its grid-avoidability.
* Larson's interaction forms (her ω² and ω^ω proofs) are the staircase tile-pair
  classes of this repo, verified mechanically class-by-class. The partial-sum "scheme"
  data her proof carries — and which our pattern-killer computation showed is
  *essential* (no pure order-pattern witness exists at ω³) — is staircase-profile/
  cut-space data. The cut ⊕ cycle split shows up as: CNF split level = cut, in-level
  shuffle = cycle.
* K₃ is the first obstruction on both sides of the repo now: blue triangles here,
  Schur boils and C₃/C₄ rungs in the additive sessions, the triangle-free ribs of
  THM-A in the metagraph. "Everything is the triangle" acquires a transfinite face.

## The procedural lesson (MISTAKE-067)

The first CEGAR implementation reported R(2,2) > 14. The number was wrong because the
final "no independent subgrid" check committed greedily to first subtrees — and the
falsity became visible only when the *witness structure was printed and read*: its B₁
visibly failed to hit rectangles it had to hit. Read the witness, not the verdict.
With the complete verifier the answer collapsed to the clean R(2,2) = 5. The pattern
(a structureless persistence that evaporates under exact verification, leaving a small
clean constant) is the same shape as MISTAKE-062/063's lesson from the monad-explorer
sessions: structurelessness at small scale is usually an artifact of a wrong lens or a
broken instrument, not a fact about the mathematics.

## What would make this load-bearing

If the rank-3 graded formulation can be made *exact* for ω^(ω³) — the gap semigroup
and the correct density floor written down at that level the way F1 does for ω² —
then the open case of a $1000 Erdős problem becomes a (transfinite-indexed but
structurally finite) question about composition-free graded relation families, where
both Schipperus's positive machinery (≤2 summands) and his negative constructions
(≥4) should appear as the two regimes of one density/composition trade-off. That is
the handoff (HYP-2346).
