# The deck is a Smith diagram

*opus-2026-07-14-S300. Companion to THM-767 and HYP-6835. The owner's directive:
work the r≥7 deck tilings and the regime-2 decomposition; consider squaring the
square, Smith diagrams, and Kirchhoff's circuit laws.*

## 1. The dictionary, exactly

Brooks–Smith–Stone–Tutte (1940) turned a squared rectangle into an electrical
network: horizontal segments become nodes, each square becomes a unit-resistance
wire, the square's side is the current through it, KCL is the statement that widths
match across a segment, KVL is the consistency of heights, and "the tile is a
square" is Ohm's law with R = 1. The deck of THM-761 supports the same translation,
and where each entry lands is exact, not metaphorical:

| BSST squared square | the sheet deck Z_c |
|---|---|
| the rectangle to tile | the deck Z_c over a core-safe t0-interval |
| a square tile | one exception's bad set B_a(t0) (a sliding grid-arc) |
| all tiles are SQUARES (R = 1) | all bad arcs have the SAME phase-length 2δ = 1/7 — speed changes the wiring (the AP step w_a⁻¹), never the arc length |
| Dehn: squarable ⟺ commensurable | zero-variance ⟺ 7·gcd(w,c) \| c — counts integral (= c/7), currents quantized |
| current through a wire | the exit rate: w_a boundary events per unit u |
| node (horizontal segment) | an event moment: a phase boundary exactly on a grid point |
| KCL at a node | the event-crossing barrier: with strict (open) bad sets no current crosses a node at all — exact tilings are chamber-locked; the coincidence law (14·gcd \| w_a + w_b, gcd hits/window) classifies the node types [CORRECTED per MISTAKE-146: the original absorption inequality was vacuous — mac-mini audit] |
| KVL (heights consistent) | the phase potentials θ_a = w_a·t0/c; exits at potential +c/14, entries at −c/14; the mirror is time reversal (THM-745, one level down) |
| total resistance = aspect ratio | the balance identity: each exception's mean bad count is exactly c/7 — total mean current r·c/7 vs deck size c; r = 7 is the critical aspect |
| matrix-tree integrality of currents | the coincidence count is exactly gcd(w_a, w_b) per window — arithmetic, not analysis |

## 2. The inversion that makes it a proof

BSST wanted the tiling to EXIST — perfection (all squares distinct) was the hard
constraint they engineered toward. Our situation is dual: the runners are GIVEN
distinct, and the enemy needs the tiling to PERSIST. Persistence turns out to be
governed by something even blunter than a conservation inequality [corrected per
MISTAKE-146 — mac-mini's audit]: with strict (open) danger sets, NO current crosses
a node at all — a perfect exit/entry handoff still leaves the sheet momentarily
free, so exact tilings are CHAMBER-LOCKED, and every chamber wall in the core-safe
set is a witness. The coincidence law (14·gcd | w_a + w_b, exactly gcd hits per
window) classifies the wall types; double-boundary walls are scarce (1399/1400
random 7-sets have essentially none). In BSST the network's solvability produces
the squared square; here the network's UNCROSSABLE NODES produce the lonely runner.
The same mathematics, run in reverse: **their existence theorem is our
impossibility theorem.**

And the strongest case needs no network at all. When 7 | c (the Dehn-commensurable,
"unit resistance exact" case), every count is pinned at c/7, so the deck is covered
with ZERO slack — and then any single switching event tears it: at the event moment
the eventing exception's count drops by gcd ≥ 1, the total falls to c − 1, a sheet
is freed, and the freed sheet is a full closed 1/14-witness (THM-767's event
pierce). The S299 wall instance — all seven sheets bad at the core optimum, which I
recorded yesterday as "the method's wall" — is pierced at every one of its 203
core-safe event moments. The correction to yesterday's framing deserves its own
sentence:

> **The witness does not live at the core optimum; it lives at the switching
> times.** The optimum is where the cover is healthiest. The tear happens where the
> current switches — at the nodes of the Smith diagram.

This is the third time the project has found the payload at a boundary/switching
locus rather than an interior optimum (the 7-clock corners carry the tight
witnesses, THM-754; the mirror pairing lives at event coincidences, THM-745). The
observer-lens principle sharpens into circuit language: the basepoint to watch is
not where the observed structure is extremal but where its derivative jumps.

## 3. What the frame bought, concretely

- **THM-767**: balance identity (∫|B_a| = c/7 exactly, every stratum), two-value
  lemma, zero-variance at 7g|c, the event pierce closing the r = 7 stratum above
  the explicit shape bound w_max > 7·Σ(P) (sharpened by mac-mini's audit to
  max w_a/g_a ≥ 7·max(P) via the Lipschitz interval), the coincidence law with
  exact gcd counts, and the chamber-locking barrier (which replaced the withdrawn
  absorption inequality — MISTAKE-146). Referee battery all-pass;
  `deck_event_witness()` in the library (self-test 16/16).
- **Regime 2 corrected and measured**: the raw complementarity `r_P ≤ B(c*)` died
  twice in one day, independently (codex's tooth-insertion falsifier {1..11, N};
  my scale-free growth census) — the fourth instance of the MISTAKE-140 genus, now
  logged. The surviving coordinate is peel-relative: ρ = v*/maxP, scale-invariant
  on dilates, measured ≤ 9.335 with the extremum at the {1..12} shape itself (the
  adversarial climb CONVERGED to it). One lemma remains to make regime 2 bounded:
  a |G'| floor off the classified tight families — pure stability-around-rigidity,
  the B5 lane.
- **The program** (HYP-6835): the chamber census (fully-covered chambers, wall
  types via the coincidence law — the deck's AP/GW corners), the s-threshold decks
  (7∤c) via the closed Kronecker geodesic, the r ≥ 8 alignment residue, Lean
  targets (all counting, decide-friendly).

## 4. The recursion, one level deeper

S299's reflection said: the underlying object is a pointed circle with burned arcs,
self-similar under scale descent, and tight instances are tilings at every level.
S300 adds the dynamic layer: at every level the tilings are not static objects but
MAINTAINED covers under a phase flow, and their persistence is governed by the
event calculus. The circle carries arcs; the arcs carry currents; the currents
cannot cross their own switches — so wherever the zero-slack regime holds, the
first switch inside the safe set tears the cover, and the tear is the lonely time.
Squaring the square was always about when a rectangle can be filled with no slack.
The Lonely Runner's covering case is about when a circle cannot be. The two
problems share their bookkeeping, and the bookkeeping is Kirchhoff's — with the
one deck-specific amendment mac-mini's audit supplied within two hours of the
claim: on this network, current cannot cross a node at all, and that is precisely
why the runner gets lonely.
