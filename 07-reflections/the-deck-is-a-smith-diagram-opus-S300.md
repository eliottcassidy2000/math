# The deck is a Smith diagram

*opus-2026-07-14-S300. Companion to THM-767 and HYP-6835. The owner's directive:
work the r≥7 deck tilings and the regime-2 decomposition; consider squaring the
square, Smith diagrams, and Kirchhoff's circuit laws.*

> **Correction (THM-771 / MISTAKE-146).** The Smith-diagram analogy remains a
> useful picture for owner-labelled incidence, but the promoted KCL law is not a
> theorem for strict bad arcs. Event density is controlled by reduced winding
> `w/gcd(w,c)`, not raw `w`. The exact replacement is
> `F=Q+Omega-sigma`; all theorem-facing claims below are read through THM-771.

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
| current through a wire | the reduced event rate u_a=w_a/gcd(w_a,c) in the core-time coordinate |
| node (horizontal segment) | an event moment: a phase boundary exactly on a grid point |
| KCL at a node | heuristic only: exact tilings are chamber-locked because strict endpoint equality tears rather than transfers ownership; for general covers the exact law is the incidence defect F=Q+Omega−sigma, while the coincidence law (14·gcd \| w_a + w_b, gcd hits/window) only classifies node types [CORRECTED per MISTAKE-146] |
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

And the strongest case needs no network at all. When every effective order
`C_a=c/gcd(w_a,c)` is divisible by seven (the Dehn-commensurable, "unit resistance
exact" case), every off-event count is pinned at c/7, so an exact deck cover has
zero slack. Any single switching event then tears it: at the event moment the
eventing exception's count drops by its gcd, at least one sheet is freed, and the
freed sheet is a full closed 1/14-witness (THM-771's reduced-winding event pierce).
The S299 wall instance — all seven sheets bad at the core optimum, which I
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

- **THM-767/771, corrected**: balance identity (∫|B_a| = c/7 exactly), corrected
  integral counts `c/7` and `c/7-g`, and the event pierce closing the unramified
  r=7 stratum above `max(w_a/g_a)>=7 max(P)`. The theorem-facing invariant is
  `F=Q+Omega-sigma`; the former KCL necessity is withdrawn. Exact tilings are
  chamber-locked at every scale, the coincidence law classifies wall types, and
  `deck_event_witness()` remains the constructive witness emitter.
- **Regime 2 corrected and measured**: the raw complementarity `r_P ≤ B(c*)` died
  twice in one day, independently (codex's tooth-insertion falsifier {1..11, N};
  my scale-free growth census) — the fourth instance of the MISTAKE-140 genus, now
  logged. The surviving coordinate is peel-relative: ρ = v*/maxP, scale-invariant
  on dilates, measured ≤ 9.335 with the extremum at the {1..12} shape itself (the
  adversarial climb CONVERGED to it). One lemma remains to make regime 2 bounded:
  a |G'| floor off the classified tight families — pure stability-around-rigidity,
  the B5 lane.
- **The program** (HYP-6835): a joint chamber census and ramification-surplus
  equality-packet classification, effective s-threshold decks via the closed
  Kronecker geodesic, the r ≥ 8 alignment residue, and decide-friendly Lean targets.

## 4. The recursion, one level deeper

S299's reflection said: the underlying object is a pointed circle with burned arcs,
self-similar under scale descent, and tight instances are tilings at every level.
The corrected dynamic layer is incidence rather than electrical conservation. At
every level the bad sets vary under a phase flow; at an unramified switch the strict
cover tears, while ramification surplus can pay overlap debt. The circle carries
arcs, the arcs carry owners, and the exact bookkeeping is `F=Q+Omega-sigma`.
Squaring the square remains a useful analogy for labelled tilings, but it supplies
no KCL proof law for the LRC boundary convention. Exact tilings cannot cross their
own chamber walls; general covers can persist only by spending ramification surplus
on overlap debt. That incidence-and-chamber calculus, rather than an electrical
conservation law, is the reusable payload.
