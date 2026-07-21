---
id: THM-1850
title: "DIRECTED WOWII: a PROVED domination–transitivity inequality γ(T) + tr(T) ≤ n+1, plus three refuted directed-inequality conjectures with explicit witnesses. Formulating tournament analogs of Written-on-the-Wall-II invariant inequalities (α↦tr = largest transitive subtournament, χ↦dichromatic, plus γ, fas, #C₃, H, diam) and testing them exhaustively over all iso classes n=4..7. MAIN RESULT (proved, elementary): γ(T) + tr(T) ≤ n+1, i.e. the domination number plus the largest transitive subtournament is at most n+1 — take the source a₁ of a maximum transitive subtournament A and adjoin V∖A; a₁ dominates the rest of A, so {a₁}∪(V∖A) is a dominating set of size n−tr+1. Tight exactly at the transitive tournament (γ=1, tr=n). COROLLARY (with tr ≥ n−fas): γ(T) ≤ fas(T)+1, tight on the 3-cycle. REFUTED with smallest explicit counterexamples: (C) dichr ≤ ⌈n/tr⌉ fails at n=7 (tr=4, dichr=3); (G) the naive directed-103 tr ≤ ⌊n−log(diam)⌋ fails at n=4 (tr=3, diam=3); (J) H ≤ 2^{n−tr} fails at n=4 (tr=3, H=3)."
status: >
  γ + tr ≤ n+1: PROVED (elementary, two lines) and VERIFIED — the witness set {source of a
  max transitive subtournament} ∪ (V∖A) is dominating and of size n−tr+1 for all 528 iso
  classes n ≤ 7; tight (γ = n−tr+1) in 50 of them (the near-transitive end). γ ≤ fas+1:
  PROVED corollary (tr ≥ n−fas is standard, kind-pasteur THM-1390). Three refutations:
  VERIFIED-EXACT smallest counterexamples over all classes n ≤ 7. The elementary inequality
  may be folklore; presented with the directed-WOWII framing and an in-repo proof. The other
  confirmed holds (Erdős–Moser tr ≥ ⌊log₂n⌋+1, domination γ ≤ ⌊log₂n⌋+1, tr ≥ n−fas) are
  standard and used as anchors, not claimed.
source: klein-2026-07-21-S397 (owner: work on the directed analogies of the WOWII inequalities)
depends_on: []
related:
  - THM-1805  # the 3-cycle = intransitivity atom (the directed odd-cycle obstruction)
  - THM-1390  # kind-pasteur: tr ≥ n − fas (the feedback-arc bound, used in the corollary)
  - "07-reflections/the-wowii-103-refutation-and-what-it-lends-the-repo-klein-S395.md (the idea source)"
script: 04-computation/directed_wowii_klein_S397.py (+ .out)
---

# THM-1850 — directed WOWII: a domination–transitivity inequality, and three refutations

## The directed invariant zoo

Tournament analogs of the WOWII (undirected) invariants:

| undirected | directed (tournament) |
|---|---|
| independence `α` (no edge) | **`tr(T)`** = largest **transitive** subtournament (no 3-cycle) |
| chromatic `χ` | **dichromatic `dichr(T)`** = min colours, each class acyclic |
| clique `ω` | `tr(T)` |
| domination `γ` | `γ(T)` (min dominating set) |
| — | `fas(T)` (feedback arc set), `#C₃`, `H(T)` (Rédei), `diam(T)` |

Ten candidate directed inequalities were tested **exhaustively over all iso classes `n = 4..7`**
(528 classes total).

## The proved result

> **Theorem. `γ(T) + tr(T) ≤ n+1`** — equivalently `γ(T) ≤ n − tr(T) + 1`.

*Proof.* Let `A` be a maximum transitive subtournament, `|A| = tr`, and let `a₁ ∈ A` be its
**source** (beats every other vertex of `A`, by transitivity). Put `S = {a₁} ∪ (V∖A)`. Then `S`
dominates `T`: every vertex of `V∖A` is in `S`, and every other vertex of `A` is beaten by `a₁`.
So `γ(T) ≤ |S| = 1 + (n − tr) = n − tr + 1`. ∎

Verified: the witness set is dominating and of size `n−tr+1` for all 528 classes; **tight**
(`γ = n−tr+1`) in 50, all at the near-transitive end. It is tight at the transitive tournament
itself (`γ=1`, `tr=n`, `1 = n−n+1`).

**Corollary (`γ ≤ fas + 1`).** With the standard `tr ≥ n − fas` (THM-1390): `γ ≤ n − tr + 1 ≤
n − (n−fas) + 1 = fas + 1`. Tight on the 3-cycle (`fas = 1`, `γ = 2`). So the domination number
is bounded by the feedback arc set plus one — a clean directed inequality chaining transitivity,
domination, and acyclicity: `n − fas ≤ tr ≤ n − γ + 1`.

The interesting regime is **near-transitive** `T` (large `tr`), where `γ ≤ n − tr + 1` is a
genuine constraint; for intransitive `T` the domination bound `γ ≤ ⌊log₂ n⌋ + 1` is far smaller.

## Three refuted directed conjectures — explicit witnesses

The WOWII discipline: state the inequality, search for the smallest counterexample.

| conjecture | smallest counterexample |
|---|---|
| **(C)** `dichr ≤ ⌈n/tr⌉` | `n=7`: `tr=4`, `dichr=3`, but `⌈7/4⌉ = 2`. A tournament with a transitive quarter still needs 3 acyclic colours. |
| **(G)** directed-103 `tr ≤ ⌊n − log(diam)⌋` | `n=4`: `tr=3`, `diam=3`, `⌊4 − log3⌋ = 2 < 3`. The naïve `α↦tr`, `ecc↦diam` port of WOWII-103 fails immediately. |
| **(J)** `H ≤ 2^{n−tr}` | `n=4`: `tr=3`, `H=3`, `2^{1}=2 < 3`. The near-transitive tournament already has more Hamiltonian paths than `2^{n−tr}`. |

Each is a genuine directed-WOWII refutation with an 11-vertex-style explicit witness — (C) shows
the dichromatic number is *not* controlled by `n/tr`; (J) shows the Rédei path count is *not*
capped by the transitivity defect exponentially.

## What this demonstrates

The directed-WOWII pipeline **works and produces content**: one proved domination–transitivity
inequality (with a corollary tying in the feedback arc set), and three clean refutations, all
exhaustive over `n ≤ 7`. The confirmed classics (Erdős–Moser `tr ≥ ⌊log₂n⌋+1`, `γ ≤ ⌊log₂n⌋+1`,
`tr ≥ n−fas`) served as anchors. The next pass should widen the invariant set (out-domination,
kings, the arborescence ranking THM-1750) and push the confirmed holds toward `n=8,9` before
conjecturing, per the repo's own small-case discipline.

*Files: `04-computation/directed_wowii_klein_S397.py` (+ `.out`).*
