---
id: THM-767
title: Deck Kirchhoff rigidity — the balance identity, the zero-variance (unit-resistance) case 7g|c, the EVENT PIERCE closing the r=7 deck stratum (the witness lives at the switching times), and the chamber-locking of exact deck tilings (the strict-set event-crossing barrier)
status: PROVED parts (1)-(3) + corollary + coincidence law; part (4)'s ABSORPTION INEQUALITY WITHDRAWN (mac-mini audit, 2026-07-14 16:00, MISTAKE-146) and REPLACED by the stronger event-crossing barrier
source: opus-2026-07-14-S300 (owner directive: r>=7 deck tilings via squaring-the-square / Smith diagrams / Kirchhoff); correction audit mac-mini-2026-07-14 (routed codex-S3)
depends_on:
  - THM-761   # the sheet frame: fiber-exact core, bad-sheet grids
related: [THM-745, THM-754, THM-760, HYP-6830, HYP-6835]
verification: 04-computation/lrc14_deck_kirchhoff_rigidity_opus_S300.py
  (+ 05-knowledge/results/lrc14_deck_kirchhoff_rigidity_opus_S300.out);
  correction referee lrc14_thm767_correction_referee_opus_S300.py (mac-mini example exact)
---

# THM-767 — Deck Kirchhoff rigidity

**Frame.** V = cP ⊔ W as in THM-761: scale c ≥ 2, core P (|P| = 13 − r), exceptions
W = {w_1, …, w_r} with c ∤ w_a, g_a = gcd(w_a, c), δ = 1/14, sheets t_k = (t0+k)/c,
bad-sheet sets B_a(t0) = {k ∈ Z_c : ||w_a(t0+k)/c|| < 1/14}. The sheets preserve the
core margin exactly (THM-761(i)); everything below concerns the exceptions' deck.

## (1) The balance identity (the conserved current)

> For every exception w and every scale c: **∫₀¹ |B_w(t0)| dt0 = c/7 exactly.**

*Proof.* With g = gcd(w, c), the phase multiset {w(t0+k)/c mod 1 : k ∈ Z_c} is the
translated grid {θ0 + j·(g/c) : j = 0, …, c/g − 1}, each point of multiplicity g, and
θ0 = w·t0/c moves at rate w/c and is periodic in t0 with period (g/c)/(w/c) = g/w.
Over t0 ∈ [0,1) there are exactly w/g ∈ Z full periods. The mean number of grid points
of spacing s in a uniformly sliding arc of length L is exactly L/s; here L = 2δ = 1/7,
s = g/c, multiplicity g, giving g·(c/(7g)) = c/7. ∎

## (2) The two-value lemma and the zero-variance case (unit resistance)

> |B_w(t0)| takes only the two values g·⌊c/(7g)⌋ and g·(⌊c/(7g)⌋+1) — a two-valued
> step function of t0. If moreover **7g | c**, the count is **constant, = c/7, for
> every t0 outside the finite event set** E_w = {t0 : some phase sits exactly on the
> bad-arc boundary ±1/14}, where it drops to c/7 − g.

*Proof.* Grid-in-open-arc counts lie in {⌊L/s⌋, ⌊L/s⌋+1} (THM-761's counting). When
L/s = c/(7g) is an integer, an open arc whose length is an exact multiple of the
spacing contains exactly L/s points unless an endpoint lies on the grid — and then
both endpoints do, giving L/s − 1. Multiplicity g scales both values. ∎

This is the Smith-diagram "unit resistance": all exceptions burn arcs of the SAME
phase-length 2δ = 1/7 (speed changes the wiring — the AP step — never the arc length),
and 7g | c is the Dehn-commensurability condition making every current integral.

## (3) THE EVENT PIERCE — the r = 7 deck stratum closes at the switching times

> Let 7 | c, r = 7, and suppose every exception satisfies **7·g_a | c** (in particular
> gcd(w_a, c) = 1 suffices). Let t0* be an **event moment** of any exception (some
> phase exactly at ±1/14) with t0* in the **closed** core-safe set
> Ḡ_P = {t : min_{p∈P} ||p t|| ≥ 1/14}. Then
> Σ_a |B_a(t0*)| ≤ 6·(c/7) + (c/7 − g) ≤ c − 1, so a **free sheet k\*** exists, and
> at t\* = (t0* + k\*)/c every exception clears ≥ 1/14 (closed), every core element
> clears ||p t0*|| ≥ 1/14 exactly, hence
>
> **M(V) ≥ 1/14.**

*Proof.* By (2) each non-eventing exception contributes at most c/7 and the eventing
one exactly c/7 − g_a ≤ c/7 − 1; the total is ≤ c − 1 < c, so some sheet lies in no
B_a — i.e., every exception's clearance at that sheet is ≥ 1/14 in the closed sense
(the bad condition is the open inequality). The core is fiber-exact. ∎

**Event density (when events must exist).** The event moments of w form two shifted
(1/w)-grids in t0 (exits at w·u ≡ c/14, entries at w·u ≡ −c/14, mod c, u = t0 + k;
for gcd(w, c) = 1 their t0-projections are (1/w)-grids), so consecutive events of w
are ≤ 1/w apart. Ḡ_P has measure ≥ 1 − |P|/7 = 1/7 (six core speeds each burn closed
measure 2δ = 1/7) and at most Σ_P := Σ_{p∈P} p components (complement of ≤ Σ_P arcs),
hence contains a closed interval of length ≥ 1/(7Σ_P). Therefore:

> **Corollary (large exceptions self-destruct).** 7 | c, r = 7, strata 7g_a | c:
> if max_a w_a > 7·Σ_P, then M(V) ≥ 1/14 — no search, no enumeration.
> *Sharpened constant (mac-mini audit): the closed core-safe set contains the
> Lipschitz interval of half-width (M(P) − 1/14)/max(P) around the core optimum —
> for |P| = 6 that is length ≥ 2(1/7 − 1/14)/max(P) = 1/(7·max P) — and the
> stratified event mesh is g_a/w_a, so* **max_a (w_a/g_a) ≥ 7·max(P)** *suffices.*

**The S299 wall falls.** The realized wall instance (c = 7, P = {1..6},
W = {12, 38, 72, 96, 151, 169, 188}: ALL sheets bad at the core optimum t0 = 1/7) is
pierced at **every one of its 203 tested core-safe event moments**, each yielding a
full exact 1/14-witness. THM-761's residual note stands corrected: standing at the
core OPTIMUM was the artifact; **the witness lives at the switching times** — the
optimum is where the deck is most covered; the boundary events are where it tears.

## (4) The coincidence law and the event-crossing barrier
## [CORRECTED 2026-07-14 ~16:15 — the original "KCL absorption inequality" is WITHDRAWN; see MISTAKE-146]

Exits of a live on the u-grid w_a·u ≡ c/14 (mod c) — current w_a per window of length
c; entries of b on w_b·u ≡ −c/14 (mod c). Eliminating u:

> **Coincidence law (stands as proved).** An exit of a and an entry of b share a
> u-value iff **14·d_ab | w_a + w_b** (d_ab = gcd(w_a, w_b)), and then exactly
> **d_ab** coincidences occur per window. (Verified exactly, 300 random triples,
> zero mismatches.) This classifies the deck's double-boundary events.

> **Event-crossing barrier (replaces the withdrawn inequality; STRONGER).** Because
> the bad condition is STRICT (open), a perfect exit/entry handoff still leaves the
> sheet momentarily free: at the coincidence instant both phases sit exactly at
> ±1/14, so the sheet is in neither open set, and maintaining coverage across the
> event would need a third blocker — contradicting multiplicity 1. At a
> non-coincidence exit the multiplicity simply drops to 0. Hence **an exact deck
> tiling cannot cross ANY event, at any c: exact tilings are CHAMBER-LOCKED** —
> they persist only on the open chambers between consecutive events (owner-a event
> mesh g_a/w_a in t0), and every chamber wall inside the closed core-safe set is a
> witness moment.

**What was withdrawn and why (mac-mini audit, routed codex-S3).** The original text
stated: "a maintained exact tiling whose interval contains the full event window
requires Σ_{14·d | w_a+w_b} d_ab ≥ w_a." The barrier above makes that hypothesis
UNSATISFIABLE (no exact tiling survives its first event), so the inequality was
vacuous as scoped — and misleading: the audit's chamber example (c = 7, core {1..6},
W = {1, 4, 5, 6, 8, 9, 10}, t0 = 1/7) is an exact tiling with all-singleton bad
sheets persisting on its whole chamber while the a = 10 absorption capacity is 2
under the plain-14 form and 0 under the 14·gcd form — capacity does not govern
chamber persistence (nothing does; chambers are trivially static). Exact referee of
the example: tiling confirmed; the nearest event (owner w = 10, t0 = 3/20, distance
1/140) pierces the cover exactly as (3) predicts, with full witness t = 43/140 at
clearance exactly 1/14. The mirror arithmetic (14-divisibility of w_a + w_b — the
deck's time-reversal, THM-745 one level down) remains meaningful as the
classification of WHERE double-boundary events occur, feeding the chamber census;
the "1399/1400 generic violation" statistic describes coincidence SCARCITY, not a
persistence criterion.

## (5) Scope, residual, and what this closes

- **Closes:** the r = 7 stratum of THM-761's residual, above the explicit shape bound
  w_max > 7Σ_P, whenever V admits a scale with 7 | c, exactly six c-multiples, and
  7g_a | c strata (coprime automatic for c = 7·squarefree-part tricks). In
  particular every family with exactly six multiples of 7 and one huge non-multiple
  is now closed by an event witness.
- **Residual (honest):** (i) r ≥ 8 — the single-event pierce fails structurally
  (Σ ≥ c persists; a surviving cover at an event moment is realized in the battery:
  P = {5,7,8,13,14}, W = {108,169,143,213,206,197,30,162}, t0 = 19/216); needs
  (r−7)·(c/7)/g + 1-fold simultaneous events. (ii) 7 ∤ c decks — the s-threshold
  alignment problem; there the enemy has slack, covers (not exact tilings) can cross
  events via multiplicity buffers, and the chamber census + s-threshold analysis of
  HYP-6835 governs. (iii) the bounded-shape residue: all w_a/g_a ≤ 7·max(P) with
  7 | c — a normalized finite family per shape, exactly the right feed for the
  band/blocker machinery. (iv) strata with 7g_a ∤ c — descend c → 7·lcm-adjusted
  scales or handle by THM-761's budget.
- **The tournament/engineering thread:** the deck's incidence is r circulant rows;
  the event calculus is its boundary operator; the event arrangement partitions t0
  into chambers (the deck's Smith diagram: nodes = events, chambers = wires), exact
  tilings are chamber-locked, and every chamber wall in the closed core-safe set is
  a witness. The free sheet is the source; the pierce says the source appears the
  moment any current switches.
