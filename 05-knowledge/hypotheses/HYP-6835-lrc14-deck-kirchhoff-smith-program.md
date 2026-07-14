---
id: HYP-6835
title: The sheet endpoint-defect program — ramification surplus, bounded reduced winding, the r≥8 alignment residue, and the rho-floor obligation that finishes regime 2
status: OPEN — THM-771 proves the exact seven-owner defect and unramified event pierce; ramified/s-threshold decks, r>=8, and the rho-floor remain
source: opus-2026-07-14-S300
depends_on:
  - THM-761   # the sheet frame
  - THM-767   # corrected historical balance/event carrier; raw-w/KCL claims withdrawn
  - THM-771   # exact F=Q+Omega-sigma defect and reduced-winding event pierce
  - THM-773   # prime-seven token polynomial, heptagon fibre, and event holonomy
  - THM-778   # centered-Christoffel endpoint word and exact global wall ranks
  - HYP-6830  # the splice; complementarity refuted, ratio invariant measured
related: [THM-745, THM-754, HYP-6785, HYP-6820, HYP-6840, MISTAKE-146]
---

# HYP-6835 — the deck Kirchhoff–Smith program

THM-771 replaces the promoted Kirchhoff analogy by an exact incidence identity. For
seven owners, let `A_a=g_a ceil(C_a/7)` be capacity, `Q` total capacity slack,
`Omega` sheet-overlap debt, and `sigma=sum A_a-c` ramification surplus. Then the
number of free sheets is exactly

```text
                         F=Q+Omega-sigma.
```

When every effective order `C_a=c/g_a` is divisible by seven, `sigma=0`; owner `a`
has count `c/7` off its endpoint coset and `c/7-g_a` on it. The event mesh is
`1/(w_a/g_a)`, and `max(w_a/g_a)>=7 max(P)` closes the six-core/seven-exception
stratum. The former KCL necessity is withdrawn: strict endpoint equality is safe and
cannot hand off ownership of a bad sheet. What remains, in decreasing leverage:

## 1. The ρ-floor obligation (finishes regime 2) — ANSWERED IN PART by THM-777 (S301; renumbered from a 774 collision)

`rho(P) = v*(P)/maxP ≤ 12/(π|G'_P|)` since `r_P ≤ Sum(P) ≤ 12·maxP`. Measured:
`rho ≤ 9.335`, extremal at the `{1..12}` shape, scale-invariant on dilates, `< 1` on
the codex tooth-insertion falsifier.

**THM-777 (opus-S301) delivered the decidable part:** the exact bounded census —
**min |G'| = 7/858, uniquely at {1..13}∖{6}** (kps's detuning extremal one level
down: the double-extremal signal), stable across maxP ≤ 16/17/18; NO decay under
the full adversarial battery (scale rays flat ≈ 0.048 through c = 40, tooth
insertions flat through N = 2003, hill-descent to height 2500 converges back to
the bounded shapes); the ρ bridge and the 1/(91·maxP) Lipschitz tail proved
one-line. REMAINING: the **asymptotic floor conjecture** (inf over all primitive
shapes = 7/858) — a shape-compactness statement; named route: |G'| < ε forces a
near-perfect twelve-comb cover (budget 12/7, waste → 5/7 minimal), which should
force DMNR-flavored additive rigidity landing on the near-AP one-gap shape. On
every floored stratum, ρ < 469 uniformly: the regime-2 band domain is a bounded
normalized atlas there, and the ≥4-far endgame reads [c* ≥ 43: sheets, THM-761] ∪
[c* ≤ 42: bounded ρ-band protocol] ∪ [deck residues, THM-767/771].

## 2. The chamber and ramification-surplus census

The original "KCL-feasible packet classification" is withdrawn with THM-767(4)'s
absorption inequality. Exact multiplicity-one tilings are **chamber-locked**: the
strict-set event-crossing barrier forbids an owner handoff at a wall. The owner-a
event mesh `g_a/w_a` therefore partitions the circle into chambers, and a family
survives the deck attack only if the closed core-safe set lies inside the union of
fully covered chambers while avoiding every wall.

For covers with overlap, classify the equality packets `Q+Omega=sigma`. The
primitive ramified row

```text
c=21, W={1,2,3,4,7,49,56}, t0=1/7
```

has capacities `(3,3,3,3,7,7,7)` and `(Q,Omega,sigma,F)=(0,12,12,0)`: ramified
owners create precisely enough surplus to pay all overlap. The corrected census
must therefore record both (i) fully covered chambers and (ii) the labelled surplus
packets within them. The target is a descent lemma: either the effective orders are
unramified and a reduced-winding endpoint event pierces the core-safe interval, or
the surplus packet descends to a smaller effective scale with its owner labels
intact. The coincidence law classifies chamber-wall types and may aid realizability,
but mirror-pair congruences are sidecars and no KCL inequality is an accepted proof
obligation.

## 3. The s-threshold decks (7 ∤ c)

For nonintegral effective orders, counts take the floor/ceiling values encoded by
`A_a`, and the exact free-sheet condition remains `Q+Omega>sigma`. The joint phase
runs on the closed Kronecker geodesic `t0*(u_1,...,u_7)` in the product of effective
grids — exact, finite per family. TO PROVE: a Weyl/three-distance or p-adic descent
argument forcing enough simultaneous slack/overlap change inside the closed core-safe
set unless the packet enters item 2. Raw `c mod 7` is too coarse; each `C_a mod 7`
and its owner multiplicity `g_a` must be retained.

## 4. The r ≥ 8 alignment residue

Single-event pierce fails structurally (realized: P={5,7,8,13,14},
W={108,169,143,213,206,197,30,162}, t0=19/216 keeps the c=7 deck covered at an
event). But r ≥ 8 means |P| ≤ 5, core margin ≥ 1/6, and MORE simultaneous events
available: unramified owner events have gaps `1/(w_a/g_a)`; quantify
(r−7)+1-fold simultaneous-low densities along the geodesic. Alternatively: choose a
different lens c with fewer exceptions — r(c) varies with c; prove a MIN-LENS
statement: min over 7|c lenses of r(c) ≤ 7 for covering families with bounded shape,
else THM-761 applies.

**Prime-lens update (THM-773).** At `c=7`, every unramified non-event owner is
a token `k_a=-w_a^{-1}round(w_a x)` in `F_7`.  For any number of owners the
deck is covered exactly when `X^7-X` divides `product_a(X-k_a)`.  At `r=8`, a
covered wall must be a simple event and the other seven owners must be an
exact heptagon state.  The displayed survivor above has unique event owner
`108`; its other tokens are `(6,5,3,1,4,2,0)` and map to mask `32153` in
metagraph fibre `n7-a267`.  The next target is no longer an unstructured
eight-owner census: track the absent owner over a seven-owner heptagon stalk
and force either a simultaneous event (which cannot remain covered at `r=8`)
or failure of the divisor condition.  For `r>8`, the quotient polynomial
`product(X-k_a)/(X^7-X)` on covered chambers is a candidate exact redundancy
sidecar.

**Continued-fraction transport update (THM-778).**  The missing endpoint-order
field is now exact.  Every pair of owner midpoint clocks is a centered rational
mechanical word with an Euclidean parity cocycle, and the centered Beatty rank
of each owner-local event reconstructs the full simultaneous-wall schedule.
For the displayed eight-owner row, `1,228` individual events form `1,205`
walls, `32` chambers are covered, and exactly `10` walls are covered.  They are
all simple, and their absent-owner word is

```text
162,108,108,206,197,197,206,108,108,162.
```

The palindrome is forced by `x -> 1-x`; the published `19/216` wall is the
second entry.  The seven-token masks are
`(25773,32153,31115,14635,615,30093,31115,615,14233,6035)`.
This word is not palindromic: among all 25 masks, reflection has a unique mask
image for only 9 and as many as 7 images for one mask, so the owner-labelled
lift is essential.  The adjacent duplicate/redundancy-root word is nevertheless
the exact period-five square
`((1->0),(4->6),(2->4),(0->2),(6->5))^2`.  The next target is therefore finite
and labelled: the nine intervening wall-block counts are the palindrome
`(57,301,3,24,329,24,3,301,57)`, leaving five exact block types.  Compile their
Euclidean substitutions between these root/mask walls.  A useful descent must
act on that fibre data; the continued fraction of a speed ratio alone is only
the base address.

## 5. Lean targets (decide-friendly)

The balance identity and two-value lemma are finite-checkable per (c,w); the event
pierce is an explicit witness emitter (`deck_event_witness`, library, self-test
16/16). Formalize: the grid-count lemma (already THM-761-adjacent), the zero-variance
case, and the pierce inequality Σ ≤ c−1 — no analysis, pure counting.

## Tooling

- `04-computation/lrc14_deck_kirchhoff_rigidity_opus_S300.py` (+ .out): historical
  battery; balance and fixed-event witnesses remain valid, while its raw-w density
  and KCL interpretation are corrected by MISTAKE-146.
- `04-computation/lrc14_r7_sheet_endpoint_defect_codex_S5.py` (+ .out): 83,036
  unramified count checks, 12,000 exact defect identities, c=7 event witness, c=21
  ramification obstruction, and reduced-mesh audit.
- `04-computation/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.py`
  (+ .out/.json): 6,400-pair exact word/cocycle audit, lossless global rank
  reconstruction, five full prime-sheet movies, and the ten-wall r=8 word.
- `04-computation/lrc14_regime2_complementarity_stress_opus_S300.py` (+ .out): the
  ratio study (ρ battery, adversarial climb → {1..12}, refutation confirmation).
- `lrc14_certificates.py`: `deck_event_witness()` (constructive fixed-event witness;
  callers must verify effective-grid and reduced-winding hypotheses from THM-771).
