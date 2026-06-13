---
source: oracle-2026-06-01-S535o
status: synthesis + computation (a spectrum of LRC->structure mappings ranked by realizable-class restriction)
tags:
  - lonely-runner
  - tournament-mappings
  - near-graph
  - circular-indifference
  - metric
  - isomorphism-classes
  - multitude
---

# Mapping LRC into Structure Space: a Spectrum Ranked by How Much It Restricts the Realizable Iso-Classes

The S518 mapping (runner system -> closed walk on the circular tournament menu) is
*faithful but non-restrictive*: every circular class is realizable, so "which
classes are exhibitable" has the useless answer "all of them," and the conjecture
hides in the WALK (S519). This note answers the request for **different mappings
where the realizable iso-class set is more meaningfully restricted** — and isolates
the principle that controls restriction.

## The controlling principle: restriction = retained metric

A plain tournament keeps only the **order** (the sign of `(v_i-v_j)t mod 1` in
`(0,1/2)`). Loneliness is **metric** (the `1/n` threshold, S529). So:

> **The realizable iso-class set shrinks monotonically as the mapping retains more
> metric information.** Order-only mappings can't restrict (the conjecture is in the
> walk); metric-bearing mappings restrict, and LRC becomes a *membership* statement.

We build a spectrum and measure `R = |realizable iso-classes| / |all iso-classes of
the target|` (`lrc_mappings_restriction_spectrum_s535.py`).

## The computed rungs

**M0 — order tournament (baseline, S518).** Target: tournaments on `n` vertices.
Realizable = circular menu `2·Fib(n-2)`; but over `t` *all* of them appear, so per
speed-set it is a walk, not a restricted set. `R` looks small but the SET is
unrestrictive. The conjecture is not a membership statement here.

**M1 — NEAR-GRAPH (threshold `1/n`) — the key rung.** Map the runner system at time
`t` to the graph on `{observer}∪{runners}` with `i ~ j` iff circular distance
`< 1/n`. As `t` varies, a speed set realizes a *family* of near-graphs. Computed:

```
 n   realizable near-graph iso-classes / all graphs      R     LRC: observer-isolated exhibited
 4            8 / 11                                     0.73         120/121
 5           20 / 34                                     0.59         120/121
 6           54 / 156                                    0.35         60/61
```

The realizable near-graphs are exactly the **circular indifference graphs** of `n`
points at threshold `1/n` — a genuinely restricted family whose share **shrinks with
`n`** (`0.73 → 0.59 → 0.35`). And the reformulation is clean:

> **LRC@n ⟺ the marked "observer-isolated" near-graph lies in the realizable family
> of *every* speed set** (observer isolated in the near-graph = no runner within
> `1/n` = lonely).

The single misses (`120/121`, `60/61`) are **not** failures: they are the **tight AP
/ regular-polygon** sets, whose lonely instant occurs only at the measure-zero times
`t = k/n` that the finite `t`-grid skips. So the membership statement holds for every
non-tight set and is *exactly tight* (boundary-only) at the regular polygon — the
mapping fingerprints the extremal.

**M2 — `k`-level metric tournament — extreme restriction.** Color each edge by the
distance level `{near (<1/n), mid (1/n..2/n), far (≥2/n)}`. Target: 3-edge-colored
complete graphs. Realizable: `14` (n=4) and `75` (n=5) iso-classes — against
`3^6=729` and `3^10=59049` labeled colorings. The realizable colorings are the
**circular-metric** ones (a circular analogue of a metric space / Robinson matrix);
LRC = "the observer's edges are all at the `far` level" is exhibitable. Restriction
is near-total; the cost is a richer target.

**M3 — resonance / QR (Paley) tournament — arithmetic restriction.** Static map
`i -> j iff (v_i - v_j)` is a QR mod `n` (`n` prime). Realizable collapses to a
*single* Paley-type iso-class in the sample (`R = 0.25` at n=5, `0.018` at n=7).
**Caveat (honest):** for `n ≡ 1 (mod 4)` `-1` is a QR, so the relation is symmetric
(Paley *graph*, not a tournament) — degenerate; the genuine Paley tournament needs
`n ≡ 3 (mod 4)`. So M3 is extremely restrictive but arithmetically delicate, and
loses the dynamical content. It is the "all-metric-is-arithmetic" limit.

## The multitude (further mappings, stated as hypotheses)

Beyond the computed three, a family of more abstract mappings, each retaining a
different slice of metric/arithmetic and each yielding a restricted realizable set:

- **MAP-p (p-adic valuation tournament).** Edge label `v_p(v_i - v_j)` for `p | n*`
  (the prime-power filtration of S534/n=18). Realizable = **ultrametric** label
  patterns (every triangle isosceles) — strongly restricted. LRC reads the deepest
  valuation level; the channels of S533/S534 are its quotients.
- **MAP-nerve (Čech complex of the danger cover).** Map to the simplicial complex of
  which danger-arc families `B_i = {t: ||v_i t||<1/n}` have common `t`. LRC ⟺ the
  nerve does **not** cover `[0,1)` (`H_0` of the complement survives). Realizable
  nerves restricted by a Helly-on-the-circle / periodic-arc constraint.
- **MAP-wire (allowable sequence / wiring diagram).** The cyclic sequence of runner
  orderings is a **closed allowable sequence** (adjacent transpositions at crossings)
  — the wiring diagram of the runner lines. Realizable = **stretchable** sequences
  (line, not pseudoline, arrangements). LRC = the sequence contains a wide-gap
  (observer-isolated) frame. Restriction = stretchability + the linear-flow (common
  point) degeneracy.
- **MAP-matroid (resonance independence).** The speeds define which runner subsets
  are resonance-independent (`Σ m_i v_i = 0` unsolvable with bounded support); the
  realizable **matroids** are restricted (representable over `Q`); LRC couples to the
  matroid's connectivity / the inside-debt support (S533).
- **MAP-CF (continued fraction / Stern–Brocot).** Speed ratios -> CF tree; the
  three-gap (Steinhaus) theorem governs the realizable gap structures; LRC =
  the largest-gap (apex, S530) staying `≥ 2/n` along the CF descent. Restricted by
  the bounded-partial-quotient structure.
- **MAP-Robinson (seriation matrix).** The pairwise time-averaged closeness matrix is
  a **circular Robinson** matrix (a circular seriation); realizable = circular-
  Robinsonian, a classified family. LRC = the observer row is "uniformly far."

## The ranking and the lesson

```
mapping     metric retained         realizable restriction        LRC = membership?
M0 order    sign only               none (walk hides it)          no
M1 near     1/n threshold (1 bit)   circular indifference, R↓ with n   YES (the clean one)
M2 k-level  distance levels         circular-metric (near-total)  yes (richer target)
M3 QR       arithmetic of diffs     single Paley class (delicate) degenerate
MAP-p/nerve/wire/matroid/CF/Robinson  various                    conjectured restricted
```

> **The honest synthesis: the realizable-class set is meaningfully restricted exactly
> to the extent the mapping carries the `1/n` metric.** The near-graph (M1) is the
> sweet spot — minimal extra structure (one threshold bit) over the order tournament,
> yet it turns LRC into a clean membership statement on a shrinking family (circular
> indifference graphs), with the regular polygon as the boundary-tight extremal. The
> order-only menu (M0) cannot restrict, because the metric it discards is the entire
> content of the conjecture (S519/S529 reconfirmed, now as a restriction theorem).

## Verdict / next
- Computed: a spectrum M0–M3 with restriction ratios; the near-graph (M1) gives LRC
  as membership on circular indifference graphs (`R` shrinking with `n`), the AP set
  as the boundary extremal.
- Multitude of further mappings (p-adic, nerve, wiring, matroid, CF, Robinson) posed
  as restricted-realizability hypotheses.
- Concrete next: (1) prove the realizable near-graphs are *exactly* circular
  indifference graphs and characterize the "observer-isolable" subfamily; (2) compute
  MAP-wire (stretchable closed allowable sequences) realizable counts; (3) tie MAP-p
  to the S534 prime-power channels.

## Artifacts
```
04-computation/lrc_mappings_restriction_spectrum_s535.py
05-knowledge/results/lrc_mappings_restriction_spectrum_s535.out
```
Related: S518 (circular menu mapping), S519 (walk non-reducing), S529 (metric/inside
debt), S530 (apex = largest gap), S533/S534 (channels, p-adic), S521o (permutohedron).
