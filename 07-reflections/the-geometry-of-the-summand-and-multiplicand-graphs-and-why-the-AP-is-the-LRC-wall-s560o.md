---
source: oracle-2026-06-02-S560o
status: progress (geometry of the summand vs multiplicand graphs; the AP is the joint extremum = the LRC wall; sumset excess correlates with looseness, r=+0.78)
tags:
  - lonely-runner
  - summand-graph
  - multiplicand-graph
  - freiman
  - sumset
  - geometry
  - pinch
  - sieve
---

# The Geometry of the Summand and Multiplicand Graphs — and Why the AP Is the LRC Wall

Continuing the user's thread (the `A+B=C` summand graph, S559o). Two graphs on ℕ, two
foliations of the plane, related by `log` — and the LRC extremal (the AP) is the point
where *both* foliations degenerate.

## Two graphs, two foliations

| | SUMMAND graph (`a,b → a+b`) | MULTIPLICAND graph (`a,b → ab`) |
|---|---|---|
| level set of node `n` | **antidiagonal** `x+y=n` (slope −1 line) | **hyperbola** `xy=n` |
| in-pair count of `n` | `⌊(n−1)/2⌋` — **linear, dense** | `(τ(n)−[n □])/2` — **tiny, sparse** |
| leaves / sources | only `1,2` | every **prime** (`(1,p)` only; drop the unit ⇒ sources) |
| character | **additive**, complete (its simple shadow is the transitive tournament) | **multiplicative**, the divisibility DAG (sparse incomplete tournament) |

**The log bridge.** `(x,y) ↦ (log x, log y)` carries the hyperbola `xy=n` to the line
`log x + log y = log n`. So **the multiplicand graph *is* the summand graph in
log-coordinates** — the additive antidiagonals are the multiplicative hyperbolas under
`exp`. (This is the `+`/`×` weld of the hyperoperation tower, S548: `∏ = exp(Σ log)`.)
The summand graph's in-degree grows like `n/2`; the multiplicand's stays `O(d(n))` and
spikes only at highly composite `n` — the additive world is dense and uniform, the
multiplicative world sparse and arithmetic.

## The LRC dictionary (S555o / S557)

- A **pinch denominator** `C = v_a + v_b` is a **summand-graph node**: the pinch
  witnesses `t = m/C` live on antidiagonal `C`, and the denominators available to a speed
  set `S` are exactly the **sumset `S+S`** (distinct pair-sums).
- A runner `w` is **cleared** at `t=m/C` (gcd `m,C`=1) iff `C ∤ w` — a **multiplicand
  (divisibility) test**. The sieve THM-369 is the multiplicand shadow.
- oracle-S555o: **the rational pinch IS the sieve.** Geometrically, the two foliations
  **coincide on the integer points** (rational `t`); the **fine pinch (`C>n`) is where
  addition outruns division** — the antidiagonals carry witnesses the divisibility
  hyperbolas cannot see. That is the open core.

## Why the AP is the wall: the joint extremum

The two graphs explain the AP from both sides at once.

**Additive (Freiman).** For distinct-pair sums, `|S \hat+ S| ≥ 2k−3`, with **equality iff
`S` is an arithmetic progression.** So the **AP minimizes the number of distinct pinch
denominators** — the antidiagonals carrying witnesses are as few as possible — giving the
**tightest pinch pigeonhole.** Verified (`lrc_summand_multiplicand_geometry_s560.py`,
n=14): the AP attains `|S\hat+S| = 23 = 2·13−3`; random primitive 13-sets have `52–64`
(excess `27–39`). And **the correlation between sumset excess and the loneliness margin
`M(S)` is `r = +0.78`** — bigger sumset ⇒ looser (more pinch denominators = more room).
The AP is the tightest set (`M = 1/14`) and the unique sumset minimizer.

**Multiplicative (coverage).** The AP `{1,…,13}` contains a multiple of **every** `q ≤ 13`
(`coverage 12/13`, the maximum among bounded sets), so it is the **most sieve-covered**
set — it escapes the divisibility sieve only at the single modulus `q = n = 14` (no
multiple of 14 in `{1,…,13}`). The AP is the small-divisor-coverage maximizer.

> **The AP is the JOINT EXTREMUM:** minimal in the summand graph (fewest pinch
> denominators) and maximal in the multiplicand graph (densest small-divisor coverage).
> Both degeneracies say the same thing — the AP gives the proof the least additive room
> and the most multiplicative obstruction. That is *why* it is the universal hard case
> (the wall, S551; the polynomial-method failure apex, S559/HYP-2063).

## How this simplifies LRC

1. **Index the pinch search by the sumset, test with divisibility.** The two graphs split
   the problem: **addition supplies the candidate witness times** (`t = m/C`, `C ∈ S+S`),
   **multiplication tests clearance** (`C ∤ w`). They coincide on the coarse rationals
   (`C ≤ n`, = the sieve, S555o); the work is the **fine sumset** `S+S ∩ (n, 2·max]`.
2. **Freiman gives a quantitative surplus off the AP.** Any non-AP set has
   `|S\hat+S| > 2k−3` (Freiman `3k−4` theory makes the jump explicit), so it has *strictly
   more* pinch denominators than the AP. The pinch pigeonhole has provable surplus exactly
   where the AP's tightness relaxes — the surplus is the room LRC needs, and it is a
   *combinatorial* (sumset) quantity, not a measure one.
3. **The hard case is forced to be the AP (or near it).** Both extremal characterizations
   are *rigid* (Freiman equality ⇒ AP; max coverage ⇒ `{1..k}`-like). A counterexample
   would have to be simultaneously sumset-minimal and coverage-maximal — i.e. an AP — yet
   the AP is lonely (it lacks a multiple of `n`, so `t=1/n` clears it). This is the
   summand/multiplicand-geometry form of the S556 "tension."

## Verdict / next
- **Geometry:** summand = antidiagonal foliation (dense, additive), multiplicand =
  hyperbola foliation (sparse, multiplicative), `log`-conjugate; LRC pinch lives on `S+S`,
  clearance is divisibility, they meet on the rationals.
- **AP = joint extremum** (min sumset / max coverage) = the wall; sumset excess correlates
  `+0.78` with looseness.
- Concrete next: (1) make "non-AP ⇒ extra pinch denominators clear an extra runner"
  quantitative via Freiman `3k−4` (the first jump in `|S+S|`); (2) the fine sumset
  `S+S ∩ (n,2max]` as the index set for the open fine-regime pinch pigeonhole (S555o
  handoff); (3) phrase HYP-2063's apex (even node's missing midpoint pair) as the unique
  antidiagonal/hyperbola coincidence at `n=2q`.

## Artifacts
```
04-computation/lrc_summand_multiplicand_geometry_s560.py
05-knowledge/results/lrc_summand_multiplicand_geometry_s560.out
```
Related: S559o (summand node = pinch denominator), S555o (rational pinch = sieve),
S557 (pinch radius r/s), S548 (+/× hyperoperation tower / log weld), THM-369 (sieve),
S551 (the wall), HYP-2063 (apex zero-divisor), summand-graph-fermat-zeckendorf.md,
natural-operation-digraphs-and-product-sum-s365.md.
