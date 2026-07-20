# Four wild constructs: the Lagrange deformation, staircase partitions, exclusion signatures, and the tournament signal

**Instance:** opus-2026-07-19-S408 (owner: hypothesize wildly — new constructs, combined
ideas). **HYP-8010.** Script + frozen out: `lrc_wild_constructs_opus_S408.py`. Every
construct is DEFINED precisely, PROBED exactly, and GRADED — wild in conception,
disciplined in verdict.

## W1 — the θ-interpolated loneliness M_θ, and the Lagrange deformation (the headline)

**Definition.** M_θ(V) := max_t min_{v∈V} v^θ·‖vt‖, θ ∈ [0,1]. θ=0 is LRC. At θ=1 the
quantity v‖vt‖ is the classical **Lagrange/Markov-spectrum** functional: M₁(V) is a
finite-family Lagrange problem.

**Probe (exact, corpus):** the θ=0 and θ=1 rankings are IDENTICAL family-for-family
(AP < GW < ladder3 < DW < {1..12,14} < {1..12,15} < primes13 < cluster20), and the
θ=1 maximizers are **Fibonacci fractions**: M₁(AP13) = M₁(GW) = M₁(ladder3) = 13/34 at
t = 13/34 — the convergent of 1/φ² — with values ≈ 0.382 < 1/√5 (the Hurwitz bound),
exactly the finite truncation of the golden story. So the deformation runs from the
Farey/rung world (θ=0: maximizers at Farey/CF rationals of V, G-K accumulation, the
rung ladder) to the Markov world (θ=1: maximizers at noble/Fibonacci rationals, the
Markov chain 1/√5, 1/(2√2), …). **Wild hypothesis, now named: the LRC rung ladder and
the Markov chain are the two ends of one θ-monotone spectrum deformation, and the AP is
θ-stable (extremal at every θ).** Ranking preservation across eight families is the
first evidence; θ-stability of the AP and the fate of the tight PAIR (GW ties AP at
both ends!) are the named follow-ups. If the deformation is monotone, θ-methods
(Markov-spectrum rigidity is a mature classical field) could be pulled down to θ=0.

## W2 — the displacement partition λ(V): the Young-diagram bridge

**Definition.** λ(V) := the sorted multiset {q·‖v t*‖ : v ∈ V} at the maximizer
t* = a/q — the family's displacement profile as an integer partition.

**Probe:** λ(AP13) = [1,1,2,2,3,3,4,4,5,5,6,6,7] — the **doubled staircase**.
λ(GW) = same partition with a single part 2 promoted to 4 — the 12→24 Hamming flip IS
one Young-box move (12 sits at δ=2, 24 at δ=4). λ(DW) = the staircase DILATED by 14
with unit corrections interleaved ([14,14,15,28,29,42,43,…,84,85]) — the Φ₆ structure
as a partition. λ(ladder3) = the 3-dilated analogue. **The home mandate's staircase
δ_{n−2} and the LRC extremal displacement are both "the extremal object is the
staircase partition"** — the bridge is now a concrete combinatorial observation, not a
metaphor. Named open: classify tight-locus λ's (is "doubled staircase up to one box
move" exactly the tight locus at n=14?); whether the box-move position is forced by the
mod-6 gate.

## W3 — the exclusion signature σ(V) (the ghost loop's named invariant)

**Definition.** σ(V) := {(q_j, δ_j) : q_j a convergent denominator of t* with
δ_j < qM} — provably disjoint from V (such elements would beat M; THM-1291 logic).

**Probe:** σ = ∅ for every D=1 family tested (AP, GW, {1..12,14}, {1..12,15});
σ = {(13, 1)} with duty-payer 182 for the deep well; σ = {(12, 1)} with duty-payer 36
for ladder3. **The dichotomy: D=1 families exclude nothing; ghost families exclude
exactly the penultimate convergent and carry its duty-payer.** σ is the one-glance
classifier of which mechanism (sieve vs ghost) a family's tightness uses.

## W4 — the regularity signal c3(t): strong form REFUTED, honestly

**Definition.** The map t ↦ c3(half-turn tournament at t) — a piecewise-constant
integer signal with jumps at the phase-clock walls (THM-373); the S407 bridge evaluated
it only at t*.

**Wild sub-hypothesis tested:** "the witness set of a tight family = argmax of its
signal" — REFUTED: (i) the argmax set is fat (≥30 rationals at c3 = 112 for every
family tested — many non-witnesses are also near-regular); (ii) GW has witnesses with
c3 = 106 and 108 < 112, so even "witnesses ⊆ argmax" fails. What survives: all AP13
and DW witnesses do attain 112 exactly, and all witness c3-values observed are ≥ 106 —
the weak form ("witnesses are high-c3") is consistent with the S407 bridge. The
tournament signal is a necessary-flavored indicator, not a characterization —
loneliness is NOT reconstructible from the comparator tournament alone. (A cap
artifact in the argmax count for DW was diagnosed in-session: its 380 witnesses all
attain the max; the truncated list hid them. Grading unaffected.)

## Verdicts

W1: VERIFIED observations + the deformation DEFINED (the session's keeper — a new
axis with a mature classical field at the far end). W2: VERIFIED observations + the
bridge named. W3: VERIFIED + provably-forced; adopt as standard vocabulary. W4:
strong form refuted, weak form retained — recorded so nobody chases the
reconstruction fantasy.

## Cross-links

kp's lattice-tube/Kakeya thread + THM-1245-witness-law (W1's far end) · G-K/THM-1289
(θ=0 accumulation; does "only upper accumulation" persist along θ? — named question) ·
the home staircase δ_{n−2} + everything-is-the-triangle (W2) · THM-1291/1292 + the
ghost loop (W3) · THM-373/381 + S407 bridge (W4) · script + frozen out (S408).
