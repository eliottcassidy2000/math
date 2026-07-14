---
id: THM-767
title: Deck balance, reduced-winding event pierce, and chamber locking of exact sheet tilings
status: PROVED parts (1)-(3), corrected reduced-winding corollary, coincidence law, and event-crossing barrier; original KCL absorption inequality WITHDRAWN (MISTAKE-146); exact seven-owner defect extended by THM-771
source: opus-2026-07-14-S300; correction audits mac-mini-2026-07-14 and codex-2026-07-14-S5
depends_on:
  - THM-761   # the sheet frame: fiber-exact core, bad-sheet grids
corrected_by:
  - THM-771   # exact F=Q+Omega-sigma identity, reduced event mesh, strict endpoint convention
related: [THM-745, THM-754, THM-760, THM-771, HYP-6830, HYP-6835, MISTAKE-146]
verification: 04-computation/lrc14_deck_kirchhoff_rigidity_opus_S300.py;
  correction referee 04-computation/lrc14_thm767_correction_referee_opus_S300.py;
  exact extension 04-computation/lrc14_r7_sheet_endpoint_defect_codex_S5.py
---

# THM-767 — Deck Kirchhoff rigidity

> **Cross-reference (opus-S301):** codex's THM-771 (written concurrently with this file's
> mac-mini correction) is the fuller rewrite of the seven-exception endpoint reduction —
> exact capacity/slack/overlap identity, the reduced-winding pierce, and the corrected
> scale-free cutoff. Its correction list matches the one applied here (two-value display at
> integral C/7; the g/w mesh; strict-boundary safety). For the seven-exception endpoint
> reduction use THM-771; this file remains authoritative for the balance identity, the
> zero-variance case, the fixed-event pierce under 7g_a|c, and the chamber-locking barrier.

> **CORRECTION BANNER (codex-2026-07-14-S5).**  The balance identity and the
> fixed-event pierce under `7*gcd(w,c)|c` are valid.  Three promoted claims were
> not: the integral two-value display used the wrong pair, event density used
> raw `w` instead of reduced winding `w/gcd(w,c)`, and strict safe endpoint
> equality was treated as a KCL handoff.  The raw condition
> `max(w)>7*sum(P)` and the KCL necessity are withdrawn.  **Use THM-771** for
> the proved exact replacement.  The text below retains the valid historical
> derivation but is corrected at every affected point.

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

## (2) The corrected two-value lemma and zero-variance case

Write `C=c/g` and `u=w/g`.  If `C/7` is not an integer, `|B_w(t0)|`
takes the two values

```text
g*floor(C/7),  g*ceil(C/7).
```

If instead **`7|C`** (equivalently `7g|c`), the correct two values are

```text
|B_w(t0)| = c/7-g,  on E_w,
            c/7,    off E_w,                              (2.1)
E_w={t0 : u*t0=C/14 (mod 1)}.                             (2.2)
```

The earlier floor/floor-plus-one display does not apply at the integral
boundary: openness removes both aligned endpoints and produces floor-minus-one.

*Proof.* When `C/7` is nonintegral, the translated `C`-grid meets the open
arc in either the floor or ceiling number of points.  When `C=7m`, the arc is
exactly `m` grid spacings long.  It contains `m` points off alignment and
`m-1` when both endpoints align and are excluded.  Multiplicity `g` gives
(2.1).  The two signed endpoint congruences project to the same set because
their difference is `C/7=m`, proving (2.2). ∎

This is the Smith-diagram "unit resistance": all exceptions burn arcs of the SAME
phase-length 2δ = 1/7 (speed changes the wiring — the AP step — never the arc length),
and 7g | c is the Dehn-commensurability condition making every current integral.

**Necessity guardrail.**  The condition `7g_a|c` cannot be replaced by merely
`7|c`.  At `c=14`, `t0=1/8`, the seven exceptions

```text
W={7,47,3,23,11,53,49}
```

have strict bad-sheet sets

```text
{0,2,4,6,8,10,12}, {10,13}, {0,9}, {6,9}, {5,10},
{2,7}, {1,3,5,7,9,11,13}.
```

They cover all fourteen sheets with total multiplicity `24`, and every phase
has threshold margin at least `1/112`, so the overlapping cover persists on
an open `t0`-neighborhood.  Thus maintained coverage need not be an exact
tiling in incompatible gcd strata.  This row does not satisfy `7g_a|c` for
every owner and therefore does not enter the event-pierce theorem.  THM-771's
identity `F=Q+Omega-sigma` is the exact replacement: ramification surplus can
pay for persistent overlap without producing a free sheet.

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

**Corrected event density.** Formula (2.2) is one coset of the `u=w/g`
grid. It has `u` points and cyclic mesh

```text
                         1/u=g/w,                         (3.1)
```

not `1/w` unless `g=1`. Let `b=max(P)`. A core maximizer, together with
`M(P)>=1/7` for a six-speed core, gives a closed core-safe interval of length
at least `1/(7b)`: the core envelope is `b`-Lipschitz. Therefore the valid,
stronger, scale-free corollary is

> **Corrected large-reduced-winding corollary.** If `r=7`, every owner has
> `7g_a|c`, and `max_a (w_a/g_a) >= 7*max(P)`, then `M(V)>=1/14`.

The older measure/component argument also remains valid after replacing raw
`w_a` by `w_a/g_a`, but gives only the weaker threshold involving `sum(P)`.

**The S299 wall still falls.** The realized wall instance (c = 7, P = {1..6},
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

Mirror coincidence counts may still be recorded as tournament diagnostics,
but they are not a proof certificate for the strict LRC predicate.

## (5) Scope, residual, and what this closes

- **Closes:** the r = 7 stratum of THM-761's residual above the corrected
  reduced-winding shape bound `max(w_a/g_a)>=7*max(P)`, whenever V admits a
  scale with exactly six c-multiples and `7g_a|c` for every owner. In
  particular every family with exactly six multiples of 7 and one huge non-multiple
  is closed only when “huge” survives reduction by its sheet gcd.
- **Residual (honest):** (i) r ≥ 8 — the single-event pierce fails structurally
  (Σ ≥ c persists; a surviving cover at an event moment is realized in the battery:
  P = {5,7,8,13,14}, W = {108,169,143,213,206,197,30,162}, t0 = 19/216); needs
  simultaneous endpoint deficits or a surplus/overlap argument. (ii) effective
  grids `7∤C_a` — ramification surplus can absorb endpoint losses; THM-771's
  primitive `c=21` row realizes this; covers with multiplicity buffers can cross
  events, so the chamber/surplus census in HYP-6835 governs. (iii) the bounded
  reduced-shape residue `w_a/g_a<7*max(P)` — the right feed for band/blocker
  machinery. (iv) `r>=8` and nonintegral effective grids need a separate capacity
  theorem; the withdrawn KCL is not such a theorem.
- **The tournament/engineering thread:** the deck's incidence is r circulant rows;
  the event arrangement partitions `t0` into chambers, and exact tilings are
  chamber-locked. THM-771 identifies the exact preserved quantity
  `F=Q+Omega-sigma`; a tournament on private owner splices is useful diagnostically
  but destroys capacity slack, ramification surplus, and simultaneous overlap.
