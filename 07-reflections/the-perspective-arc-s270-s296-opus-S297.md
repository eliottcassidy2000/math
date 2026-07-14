# The perspective arc (S270–S296): one counting identity, twenty theorems, one Lean file

> **Frontier correction (2026-07-14):** the analytic identities and Lean ledger described here
> remain useful at their stated per-core scopes, but the original closing interpretation is false.
> HYP-6780 proves the THM-755 cutoff is scale-covariant and gives an unbounded `f=13` covering
> ray below it. Consequently this arc does not assemble LRC(14); see the corrected
> `00-navigation/LRC14-FRONTIER-2026-07-15.md` and THM-758.

*opus-2026-07-14-S297. The completed-arc synthesis. Companion to
00-navigation/LRC14-FRONTIER-2026-07-15.md.*

## The seed

The owner's prompt (S270): *a tournament is made of n−1 perspectives, each of n−1 arcs;
among these (n−1)² arcs, T(n−2) are double-counted.* For n = 14: 169 = 13 + 2·78 — the
13 single-counted **spokes** to the origin and the 78 = C(13,2) double-counted **pairs**:
the GF(2) cut ⊕ cycle split of K₁₄, the repo's oldest object (base path = cut, tiles =
cycle), landing on the Lonely Runner problem.

Everything that followed was this one decomposition applied to LRC(14), over and over:

- **the origin band** — the origin's own danger cell [−1/14, 1/14]: the incompressible
  frame; the fence that made no-wrap unconditional (S280); the trivial Fourier cap
  ‖ĉ‖ ≤ |G| (S289); the 7-clock's half-cell (S288).
- **the spokes** — the single-counted sector: per-runner peels; the jump envelope
  ‖ĉ_m‖ ≤ r/(πm); the strand crossings of the budget.
- **the pair sector** — the double-counted T(12): coherence (packs, clusters), vertex
  events at difference-runner grids, the mirror pairing, the orientation symmetrization,
  the slot-occupancy grid (13 speeds × 12 clocks).

## The theorem ledger (opus lane)

| theorem | content | frame sector |
|---|---|---|
| THM-737 | pack-clock sampling (multiplicative coherence; d ≥ 2 detuned) | pairs → one clock |
| THM-739 | cluster-clock Area bound (additive coherence; klein's shapes ∀W) | pairs → one clock |
| THM-740 | hierarchical clocks; the product-area (separation ⟹ factorization) | cross-cluster pairs |
| THM-742 | exact-fiber strand identity; arrangement constants (20× on C₁) | spokes vs strips |
| THM-743 | vertex cost = pair difference δ; buried events free | pair events |
| THM-744→745 | F-telescoping; the PAIRING THEOREM (mirror = time reversal; ρ ≡ 0) | mirror |
| THM-746/747 | the exact 3-term phase expansion; the phase sum triangulated | all three |
| THM-750 | THE CLOSED BUDGET: W(L−Area) = Φ+Q+κ exact; U1 discharged | all three |
| THM-752 | fine-comb witness at the TOOTH threshold; ratio-13 cascade | origin band duty cycle |
| THM-754 | the 7-CLOCK PARTITION (k=7 self-dual; tight witnesses = corners) | origin band = pair cell |
| THM-755 | the CAPPED ENVELOPE: (H) for v > r/(π\|G'\|) — origin cap × spoke envelope | origin × spokes |
| THM-756 | the (H)-bands closed (4,032 pairs; the AP/GW corners) | assembly |

Plus the classification runs (S286/287: the slot-occupancy zeroes and the k=7 scarcity
collapse), the frontier syntheses, and the library (lrc14_certificates.py).

## The Lean file as the arc's monument

LRCClosedBudget.lean (47 declarations, 0 sorries, kernel-pure) is organized BY the frame:
origin band (strong_no_wrap), spokes (per_crossing_exact), telescoping (F-atoms, mirror),
pair events (κ forms), grid/mirror lemmas, the budget spec — then the full analytic chain:
Raabe (the finite Poisson atom), the grid deficit, the capped-envelope kernel, the Fourier
envelopes, the spectral theorem, the overlap identity, the family assembly, and the capstone
geometric_disc_eq_discB. The last unformalized sentence in the (H)-edge was about how two
intervals overlap on a circle; it is now a theorem.

## The method patterns (what actually worked)

1. **Exact referee first, formalize second.** Every statement was Fraction-verified before
   proving (the S295 orientation catch: 8,392 configs found V(t) = RHS(−t) — an error the
   machine-checked diagonal instances were structurally blind to). The referee also caught
   MISTAKE-142 (my own unsound C₂ bookkeeping) and three budget bugs (S283).
2. **The finite reformulation.** At every "analytic" wall the arc found the exact finite
   object underneath: Poisson → Raabe; the E-functional → constant absorption; the drift →
   the strand identity; the wrap fluctuations → the Euclid tower (which then never entered).
3. **Honest negatives, logged.** The Ξ refinement that gained nothing (S277 net-zero); the
   naive fine-comb threshold missing by 1.15×; the k=7 "lemma" that turned out to be the
   problem itself (THM-754 — the discovery was WHY, not a lemma).
4. **Claim-first under fleet velocity.** ~25 ID collisions in the sprint; the first-pusher
   protocol + stub-checkpoint-then-work kept every one a rename, never a loss.
5. **The fleet compounding.** klein consumed the capped envelope within hours (S310's
   finite-decidability); mac-mini's dichotomy consumed the U-map; kps's exact-ℚ forms became
   the Lean discB definitionally. The arc's outputs became the fleet's inputs mid-flight.

## The meta-pattern (what the mathematics kept saying)

mac-mini-S82's verdict was "structure forgets measure" — every structural facet severed from
the metric residue. The arc's answer: the measure was never analytic fog; it was **exact
finite structure one level down**, reachable by riding the right clock. Coherence
(multiplicative, additive, hierarchical) converts measure to counting; the mirror converts
residuals to zero; the self-dual 7-clock converts the threshold itself into geometry
(tight witnesses = tiling corners). The problem's last inequality lived on the unique clock
whose half-cell IS the origin band — the incompressible frame measuring itself — and fell to
the two oldest perspectives (the origin's cap, the spokes' envelope) spliced at their
crossover.

The owner's counting identity was not a metaphor. It was the load-bearing decomposition,
and the file that closes the arc has its sectors as section headers.
