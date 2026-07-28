---
id: HYP-9046
title: "Signed contractive-idempotent decomposition as the LRC(14) wall-closing template"
status: >
  WILDCARD HYPOTHESIS — external-import reframe, no in-repo claim. The cited
  external theorem is CITED-ABSTRACT (fetched 2026-07-28, proof not audited
  in-repo). Decisive cheap tests specified below; nothing here modifies the
  status of any THM.
source: klein-2026-07-28-S691
related:
  - THM-636
  - HYP-6130
  - THM-731
  - THM-538
  - THM-2349
  - HYP-9040-cancellation-completeness-synthesis
external:
  - "C. Beke, M. K. Goh, H. Hatami, S. Jaffe, D. Naylor, *A characterization of idempotent Schur multipliers*, arXiv:2607.14316 (2026)."
  - "B. Green, T. Sanders, quantitative idempotent theorem (Ann. of Math. 2008) — the abelian ancestor of the mechanism."
---

# HYP-9046 — signed contractive-idempotent decomposition for the signed wall

## The reframe in one paragraph

The large-diameter LRC(14) route died at a *cancellation wall*: an unsigned
`|OffLine|` bound provably fails (THM-636 context; HYP-6130 wider-band floor),
and the recorded lesson is that **the closing estimate must be signed**.
Independently, the relation-lattice consolidation (THM-538/THM-729/THM-731)
showed both covering and density routes are sums of `Ĝ` over a lattice
`L = {a : a·E' = 0}` — i.e. the wall functional is a pairing of `Ĝ` against
the indicator of a *coset-structured boolean pattern*, which is an idempotent
in the multiplier sense. The new external theorem (arXiv:2607.14316) says:
every boolean matrix pattern of bounded Schur-multiplier norm `γ` is a
**finite SIGNED sum of at most `2^{Cγ⁶} contractive idempotents** (blocky
patterns — disjoint unions of combinatorial rectangles). Blocky patterns pair
against `Ĝ`-data as *products of one-dimensional sums*, which is exactly the
shape the repo's positive tools (Fejér bulk, THM-731's positive-definite
x-integral, THM-2048's quantized peel) know how to bound WITH SIGNS RETAINED.
So the proposal: stop trying to bound the wall pattern absolutely; decompose
the pattern itself into signed blocky pieces and bound each piece by a
rank-one (product) estimate.

## Typed connection (per protocol)

- **Source object:** the wall functional `W = Σ_{a ∈ L\{0\}} w(a)·Ĝ(a)`
  restricted to the off-line/large-diameter regime where the unsigned bound
  fails; equivalently the boolean incidence pattern
  `M[a, ξ] = 1[a·E'(ξ) = 0]` in (relation, window-position) coordinates.
- **Target object:** `W = Σ_{j≤L} σ_j ⟨B_j, Ĝ-data⟩`, `σ_j ∈ {±1}`, each
  `B_j` a disjoint-rectangle (blocky) pattern.
- **Map:** the signed decomposition supplied by the external theorem, applied
  to `M` as a Schur multiplier pattern.
- **Preserved predicate:** the exact signed value of `W` (the decomposition is
  an identity, not an inequality — no positivity is discarded, which is
  precisely what the wall requires).
- **Destroyed/lost:** piece count `L ≤ 2^{Cγ⁶}` enters multiplicatively; all
  structure of `M` beyond its coset-ring skeleton; effectivity of `C` (the
  external bound's constant is not explicit for us).
- **Needed sidecar (the real obligation):** a proof that the Schur/multiplier
  norm `γ` of the LRC pattern is bounded **uniformly in the window/height
  parameter**. If `γ` grows with `V`, the piece count explodes and the route
  is dead. This is a self-contained, testable question.
- **Cheapest decisive test:** compute a computable proxy for `γ` (the
  factorization/γ₂ norm via small SDP, or the exact Schur norm for tiny
  cases) for the pattern `M` of the known exact extremals: the deep well, the
  tight AP at its resonant ruler (hostile control — non-covering, must not be
  misused), and two or three of the 165 frontier rows at small height.
  Plot γ against V over heights ≤ 55 (THM-1290's exact zone). Monotone
  growth ⟹ kill and record. Bounded ⟹ promote to a worked lane: bound one
  blocky piece's signed contribution on THM-636's ≤6-lift family and check it
  beats the recorded unsigned failure there.

## Controls

- **Positive control:** a pure single-coset pattern (one clock) — its Schur
  norm is 1 and the decomposition is itself; the pairing must reproduce the
  known exact clock sums (THM-2057-style closures).
- **Hostile control:** the cancellation-wall witness where `|OffLine|`
  provably fails. The signed decomposition must strictly beat the unsigned
  bound there, or the hypothesis is refuted *in the only regime it was
  invented for*.

## Why this is not another syntax bridge (MISTAKE-230–235 guard)

The map is not an analogy: `M` literally *is* a boolean pattern acting by
Schur multiplication in the pairing that defines `W`, and the decomposition
is an exact identity in that algebra. The two failure modes are honest and
named: unbounded `γ`, and per-piece bounds that don't sum below the wall.
No LRC statement is claimed; this is a template plus two finite tests.

## Note on provenance

Import found during the 2026-07-28 external-puzzle intake (klein-S691); the
owner's puzzle bundle pointed at arXiv:2607.14316 in proximity to an
arithmetic-Kakeya benchmark text. The Kakeya text itself is handled in
`04-computation/ak_forcing_engine.py` and its workbench note.
