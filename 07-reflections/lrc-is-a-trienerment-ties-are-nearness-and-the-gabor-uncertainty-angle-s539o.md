---
source: oracle-2026-06-01-S539o
status: synthesis + computation (LRC pairwise structure = a t(r)ienerment; ties = nearness; the Gabor/uncertainty angle)
tags:
  - lonely-runner
  - trienerment
  - a007025
  - ties
  - near-graph
  - gabor
  - uncertainty
  - time-frequency
---

# LRC Is a T(r)ienerment: Ties Are Nearness, and the Gabor Uncertainty Angle

The right pairwise structure for LRC is not a tournament but a **t(r)ienerment**
(the repo's 3-state object: each edge `i→j`, `i←j`, or `i↔j` a dual-headed **tie**;
iso-classes `= A007025 = 7,42,582,21480` for `n=3..6`, refined by `#ties` via
`f(n,k)`, with `f(n,0) = A000568` = tournaments). The Gabor (time-frequency) angle
then controls the ties through an uncertainty trade-off.

## The LRC trienerment: ties = nearness

On `n` vertices (observer `0` + `n-1` runners), define the **LRC trienerment** `T(t)`:
```
i ↔ j  (TIE)   iff   circular distance < 1/n     (the two runners are NEAR);
i → j  / i ← j (directed) by the half-turn order  (far + ahead/behind).
```
So **a tie is exactly a near-pair** — the loneliness-relevant relation. This single
move unifies the whole mapping program:

- the **order tournament** (S518) is the tie-free idealization `f(n,0)`;
- the **near-graph** (S535/S536) is precisely the **tie-subgraph** of `T(t)`;
- the metric mappings are the choice of tie threshold;

and the loneliness target becomes a clean **local** condition:

> **Observer lonely ⟺ the observer has tie-degree `0`** (all its edges directed =
> far from every runner). So **LRC@n ⟺ every speed set reaches a `T(t)` with
> observer tie-degree `0`.**

## A pigeonhole fact: the LRC trienerment is NEVER a pure tournament

Computed realized total-tie-counts start at **1**, never `0`
(`lrc_trienerment_gabor_s539.py`). The reason is pigeonhole: `n` points on the
circle leave `n` gaps summing to `1`, so the minimum gap is `≤ 1/n` — **there is
always at least one near-pair (one tie)**. Hence the tournament slice `f(n,0)` is the
**unrealizable** tie-free limit; the LRC trienerment lives strictly in `k ≥ 1`.
Loneliness is therefore never "the whole graph is a tournament" — it is the *local*
statement that the **observer** is tie-free, even though the graph globally cannot be.

**Boundary correction (codex-2026-06-01-S542 / THM-389).**  The paragraph above
mixes two tie conventions.  Strict ties `dist < 1/n` are the right convention for
`observer tie-degree 0 <=> closed LRC loneliness`, but then globally tie-free
states do occur exactly at regular `n`-gon walls.  Weak ties `dist <= 1/n` are
the convention for the pigeonhole statement that every `n`-point configuration
has a tie.  The durable statement is therefore two-layer: strict observer ties
for local loneliness, weak/compactified ties for the global wall ledger.

## Restriction: a shrinking slice of A007025

```
 n   realizable LRC-trienerment iso-classes / A007025(n)    R       observer-tie-free (LRC)
 4              20 / 42                                     0.48          119/121
 5             100 / 582                                    0.17          80/81
```

`f(4,k) = [4,10,12,10,4,1,1]`, `f(5,k) = [12,50,107,144,131,78,41,13,4,1,1]`; the LRC
trienerment realizes the full tie-range except `k=0`, but only a **restricted,
shrinking** fraction of iso-classes (`R: 0.48 → 0.17`). The misses in the LRC count
are again the **tight AP** sets (observer tie-free only at the measure-zero `t=k/n`).

So the trienerment is a strong "complex vertex/edge object" in the S538 sense: the
third state (tie) carries the metric, and the realizable family is a small, graded
(by ties) slice of `A007025`, with `LRC = local tie-freeness of the observer`.

## The Gabor angle: ties (real clustering) ⟷ character concentration (frequency)

Loneliness is a **sharp real-space feature** — the observer is locally tie-free, a
`2/n` clearance around it (S536's empty observer cells). By the DFT duality
(S536/S537) a sharp space feature is dual to a **frequency spread**: the empty
window forces the character sums `|ĉ_m| = |Σ_j e^{2πi m x_j}|` (S529) to stay
**bounded away from zero**. Measured at the lonely time (`lrc_trienerment_gabor_s539.py`):

```
 n=5: avg max_m |ĉ_m| = 1.81  (of n-1 = 4 runners)
 n=7: avg max_m |ĉ_m| = 2.37  (of n-1 = 6 runners)
```

— a frequency **floor** at the lonely time: you cannot empty the observer's window
*and* make all character sums small. This is a discrete **uncertainty principle**:
the **tie-count** (real-space clustering / discrepancy) and the **character
concentration** (frequency) cannot both be small. The two prior dual mappings are
its two marginals — **sectors = space (S536)**, **harmonics/flows = frequency
(S537)** — and the Gabor picture is their joint phase space.

### The Gabor trienerment (posed)

Take the vertices to be the **(sector, harmonic) cells** `(k, m) ∈ {0..n-1} ×
{1..n-1}` — the discrete time-frequency atoms (short-time Fourier atoms of the runner
indicator). Decorate cell `(k,m)` with its Gabor amplitude `G_{k,m}(t) = Σ_{j∈S_k}
e^{2πi m x_j}`, and put a **trienerment** on the atoms: `→` / `←` by amplitude
order, `↔` (tie) when two atoms are jointly **unresolvable** (the uncertainty
incomparability). Then:

- the realizable Gabor trienerments are restricted by a **discrete uncertainty
  bound** (an atom cannot be sharp in both `k` and `m`);
- loneliness = the observer's space-column (`k ∈ {0, n-1}`) is empty, which by
  uncertainty lights up a spread of frequency-rows — a **specific Gabor signature**;
- it is the joint lift unifying every prior mapping: runners (S518), sectors (S536),
  flows/harmonics (S537), near-graph/ties (S535), and this trienerment.

## Verdict / next
- LRC's pairwise structure is a **t(r)ienerment** with `ties = nearness`; loneliness
  = observer tie-degree `0`; the graph globally always has `≥1` tie (pigeonhole), so
  it is never a pure tournament — the tournament is the unrealizable tie-free limit.
- Realizable LRC-trienerments are a shrinking slice of `A007025` (`R: 0.48→0.17`),
  graded by `#ties` (the repo `f(n,k)`).
- Gabor/uncertainty: loneliness (sharp space hole / local tie-freeness) forces a
  frequency floor on the character sums — a discrete uncertainty principle whose
  marginals are S536 (space) and S537 (frequency).
- Concrete next: (1) build the `(sector,harmonic)` Gabor trienerment and measure its
  uncertainty-restricted realizable classes; (2) prove the frequency floor
  `max_m |ĉ_m| ≥ c(n)` at any observer-tie-free time (a Turán/large-sieve bound);
  (3) the `f(n,k)` grading of LRC trienerments vs the inside-debt order (S533).

## Artifacts
```
04-computation/lrc_trienerment_gabor_s539.py
05-knowledge/results/lrc_trienerment_gabor_s539.out
```
Related: t(r)ienerment machinery (`trienerment_iso_count.py`, A007025/f(n,k)),
S538 (pair-tournament ties), S535 (near-graph), S536 (sectors/DFT), S537
(flows/tension), S529 (characters/inside debt).
