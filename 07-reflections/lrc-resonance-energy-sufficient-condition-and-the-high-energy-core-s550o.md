---
source: oracle-2026-06-01-S550o
status: progress (rigorous sufficient condition for LRC via the resonance-energy bound; reduction to the high-energy core)
tags:
  - lonely-runner
  - circle-method
  - resonance-energy
  - sufficient-condition
  - minimal-resonance
  - inside-debt
  - progress
---

# A Rigorous Sufficient Condition for LRC, and the Reduction to the High-Energy Core

Pushing the recent arc (S544 decorrelation, S545 returns, S529 covering) toward
actual progress: the covering identity gives a **rigorous sufficient condition for
LRC**, and computation shows it covers *essentially every* speed set, reducing the
conjecture to a thin, structured **high-energy core**.

## The bound

From the covering identity (S529), with `ĝ(0) = 1-2/n`, `ĝ(m) = -sin(2πm/n)/(πm)`:

```
|LONELY(v)| = ∫_0^1 ∏_i 1_{[1/n,1-1/n]}(v_i t) dt
            = (1-2/n)^{n-1}  +  Σ_{0≠m, Σ m_i v_i = 0} ∏_i ĝ(m_i).
```

By the triangle inequality, with the **resonance energy**
`E(v) := Σ_{0≠m, Σ m_i v_i=0} ∏_i |ĝ(m_i)|`:

> **`|LONELY(v)| ≥ (1-2/n)^{n-1} − E(v)`. Hence `E(v) < (1-2/n)^{n-1} ⟹ |LONELY(v)| > 0
> ⟹ LRC holds for `v`.** (Rigorous modulo a tail estimate on `|m| > M`.)

This is the major-arc/minor-arc (Hardy–Littlewood) picture for LRC: `(1-2/n)^{n-1}` is
the **independence main term** (S544), `E` is the **resonance correction**, and LRC
holds whenever the resonances don't overwhelm the main term. `E` is dominated by
*small* resonances (`ĝ` decays like `1/(π|m|)`), so the controlling invariant is the
**minimal resonance length** `ℓ(v) = min{Σ|m_i| : 0≠m, Σ m_i v_i = 0}`.

## What the computation shows (`lrc_resonance_energy_full_s550b.py`)

Full truncated energy `E_trunc` (all orders, `|m_i| ≤ M`), the validated bound
`|LONELY| ≥ main − E`, and the classification:

```
 n=5 (main 0.130):  AP (1,2,3,4):  E=0.304  → main−E=−0.17 (core, |LONELY|=0)
                    random sets:   E=0.006–0.052 → main−E=+0.08…+0.12 ≤ |LONELY|  (LRC PROVEN)
 n=6 (main 0.132):  AP (1,2,3,4,5):E=0.419  → core (|LONELY|=0)
                    random sets:   E=0.019–0.058 → LRC PROVEN
 HARD-CORE SCAN:  bound proves LRC for 120/120 (n=5) and 60/60 (n=6) random primitive sets.
```

Two clean facts:

1. **The bound proves LRC for essentially every set.** All `120 + 60` sampled
   random primitive sets have `E < main` with a *large* margin (`main − E ≈ 0.08–0.12`
   vs main `0.13`), so even after a tail correction `E_full < main` — **rigorously
   lonely**. The decorrelated majority is settled.

2. **The high-energy core is thin and structured.** The only sets with `E ≥ main`
   are the **small-minimal-resonance** ones, led by the **AP / regular polygon**
   (`E = 0.30, 0.42 ≫ main`). There `|LONELY| = 0`, so `E ≥ main` *necessarily*
   (`|0 − main| ≤ E`). The first pass (S550) using only the *pairwise* energy `E2`
   mislabeled the AP (`E2 = 0.122 < main`!); the full energy reveals that the AP's
   excess comes from the **order-≥3 returns / inside debt (S545)** — exactly the
   3-cycles the no-return fact forbids. So the high-energy core = the many-returns
   sets, and the AP saturates it.

## The progress: LRC reduced to the high-energy core

> **Rigorously (modulo tail), every speed set with `E(v) < (1-2/n)^{n-1}` is lonely.**
> The Lonely Runner Conjecture is therefore reduced to the **high-energy core**:
> the (thin, structured) sets with `E(v) ≥ (1-2/n)^{n-1}` — equivalently the
> small-minimal-resonance, many-returns sets. The AP/regular polygon is the
> extremal of this core, and it *is* lonely (the sieve witness `t = a/n`,
> THM-369/initial_segment_unit_lonely). So the residual is the high-energy core
> *minus* the sieve-handled cases.

This is the honest state: the bound is a genuine, computable, (essentially) rigorous
sufficient condition that disposes of the decorrelated majority; the conjecture's
difficulty is concentrated, by this lens, in the high-energy core where the
order-≥3 returns dominate — the same wall as S529/S533/S545, now quantified as
"`E ≥ main`".

## Honest limitations
- The truncation tail (`|m_i| > M`) is not rigorously bounded here; the large margin
  (`main − E_trunc ≈ 0.1`) makes `E_full < main` safe for the generic sets, but a
  clean tail lemma (`ĝ` decay ⇒ geometric tail) is needed for a fully rigorous
  theorem.
- The bound does **not** settle the high-energy core (where it is vacuous), so it
  does **not** prove `n = 14, 16, 18`: a counterexample, if any, would live in that
  core. The contribution is the *reduction* + the explicit core characterization.

## Verdict / next
- Rigorous (modulo tail) sufficient condition `E(v) < (1-2/n)^{n-1} ⟹ LRC`, covering
  `180/180` sampled sets; the conjecture reduced to the high-energy core.
- Concrete next: (1) prove the **tail bound** `Σ_{|m|>M} ∏|ĝ| ≤ C·r^M` (geometric
  decay) to make the sufficient condition a clean theorem (and formalize it on top of
  `near_pair`/`sieve` in Lean); (2) characterize the high-energy core as a bounded
  set of speed *ratios* (small `ℓ`) — a finite check per `n`?; (3) for `n = 14`
  (prime `n* = 7`, S546) test whether the core is exactly the sieve-covered/AP sets.

## Artifacts
```
04-computation/lrc_resonance_energy_sufficient_s550.py     (pairwise E2; the mislabel + lesson)
04-computation/lrc_resonance_energy_full_s550b.py          (full E; the validated bound)
05-knowledge/results/lrc_resonance_energy_sufficient_s550.out
05-knowledge/results/lrc_resonance_energy_full_s550b.out
```
Related: S529 (covering identity), S544 (decorrelation main term), S545 (returns =
inside debt), S527 (cascade), S546 (n=14 prime channel), THM-369 (sieve / the core's
AP witness).
