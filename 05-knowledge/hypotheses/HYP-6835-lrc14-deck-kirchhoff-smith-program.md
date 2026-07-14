---
id: HYP-6835
title: The deck Kirchhoff–Smith program — KCL-feasible packet classification, the s-threshold (7∤c) decks, the r≥8 alignment residue, and the ρ-floor obligation that finishes regime 2
status: OPEN — THM-767 proves the frame (balance, zero-variance, event pierce, KCL); this file holds the remaining program
source: opus-2026-07-14-S300
depends_on:
  - THM-761   # the sheet frame
  - THM-767   # deck Kirchhoff rigidity (this session)
  - HYP-6830  # the splice; complementarity refuted, ratio invariant measured
related: [THM-745, THM-754, HYP-6785, HYP-6820]
---

# HYP-6835 — the deck Kirchhoff–Smith program

THM-767 turned the r ≥ 7 deck residue into circuit theory: currents (exit rates w_a),
a conserved balance (∫|B_a| = c/7), a commensurability condition (7g|c — the Dehn
condition), an event pierce (7|c, r=7: closed), and a KCL absorption law (maintained
exact tilings need mirror-partner capacity Σ gcd(w_a,w_b)·[14d | w_a+w_b] ≥ w_a).
What remains, in decreasing leverage:

## 1. The ρ-floor obligation (finishes regime 2)

`rho(P) = v*(P)/maxP ≤ 12/(π|G'_P|)` since `r_P ≤ Sum(P) ≤ 12·maxP`. Measured:
`rho ≤ 9.335`, extremal at the `{1..12}` shape, scale-invariant on dilates, `< 1` on
the codex tooth-insertion falsifier. TO PROVE: **a positive floor on `|G'_P|` for
12-cores outside the classified tight families** (mac-mini B5 stability lane; the
rigidity THM-757/759 + the L=0 census are the anchors). That single lemma converts
the regime-2 band domain into a bounded normalized (peel-relative) atlas and, with
THM-761 + THM-767, makes the ≥4-far endgame: [c* ≥ 43: sheets] ∪ [c* ≤ 42: bounded
ρ-band protocol] ∪ [deck residues below].

## 2. The chamber census (the Smith-diagram census — CORRECTED per MISTAKE-146)

[Original "KCL-feasible packet classification" withdrawn with THM-767(4)'s
absorption inequality: exact tilings are CHAMBER-LOCKED (strict-set event-crossing
barrier), so no capacity condition governs persistence.] The corrected census
object: the EVENT ARRANGEMENT on t0 (owner-a mesh g_a/w_a) partitions the circle
into chambers; a family survives the deck attack only if the closed core-safe set
lies INSIDE the union of fully-covered chambers, avoiding every chamber wall. So
classify: (i) which chambers are fully covered (a finite lattice-point condition
per chamber); (ii) when the core-safe set can avoid all walls — it cannot once any
w_a/g_a exceeds the Lipschitz bound (THM-767 corollary, sharpened by the audit);
(iii) the coincidence law (14·gcd | w_a+w_b, exactly gcd double-boundary events per
window — stands as proved) as the classification of chamber-wall TYPES. The
matrix-tree flavor survives re-aimed: which wall-type patterns are REALIZED by
integer speed sets (cf. HYP-6785's realizability sidecar B3).

## 3. The s-threshold decks (7 ∤ c)

For c ≡ s (mod 7), s ≠ 0: counts are two-valued, coverage needs Σ X_a ≥ s with
X_a ∈ {0,1} arc-indicators of density s/7 in phase; a free sheet needs 8−s
simultaneous lows inside the closed core-safe set. The joint phase runs on the
closed Kronecker geodesic t0·(w_1,…,w_7) in T^7 — exact, finite per family. TO
PROVE: a Weyl/three-distance argument that the geodesic meets the (8−s)-low region
whenever it is nonempty and the w's are not resonance-locked; the locked cases feed
program item 2. (The event pierce is the s = 0 degenerate case where a SINGLE low
suffices.)

## 4. The r ≥ 8 alignment residue

Single-event pierce fails structurally (realized: P={5,7,8,13,14},
W={108,169,143,213,206,197,30,162}, t0=19/216 keeps the c=7 deck covered at an
event). But r ≥ 8 means |P| ≤ 5, core margin ≥ 1/6, and MORE simultaneous events
available: events of different exceptions interleave with gaps ≤ 1/w; quantify
(r−7)+1-fold simultaneous-low densities along the geodesic. Alternatively: choose a
different lens c with fewer exceptions — r(c) varies with c; prove a MIN-LENS
statement: min over 7|c lenses of r(c) ≤ 7 for covering families with bounded shape,
else THM-761 applies.

## 5. Lean targets (decide-friendly)

The balance identity and two-value lemma are finite-checkable per (c,w); the event
pierce is an explicit witness emitter (`deck_event_witness`, library, self-test
16/16). Formalize: the grid-count lemma (already THM-761-adjacent), the zero-variance
case, and the pierce inequality Σ ≤ c−1 — no analysis, pure counting.

## Tooling

- `04-computation/lrc14_deck_kirchhoff_rigidity_opus_S300.py` (+ .out): the referee
  battery (balance 98 pairs exact; two-value; zero-variance; 86/86 + 22-instance +
  203/203 wall pierces; KCL coincidence law 300 triples exact; r=8 survivor).
- `04-computation/lrc14_regime2_complementarity_stress_opus_S300.py` (+ .out): the
  ratio study (ρ battery, adversarial climb → {1..12}, refutation confirmation).
- `lrc14_certificates.py`: `deck_event_witness()` (constructive THM-767 witness).
