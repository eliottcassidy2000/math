---
source: opus-2026-06-03-S575 (remote-control)
status: NEW exact criterion (Lemma D) + circuit-positivity functional P(S) (P>0 <=> loose, verified 100%) + a PROVED summed corollary; residual = simultaneous-resonance infeasibility
tags: [LRC, endpoint-cover, circuit, positivity, midpoint, congruence, THM-398, C-prime, n14]
---

# Endpoint-cover circuit positivity

**Prompt (user):** work on endpoint cover circuit positivity.

The S574 translator turned "is component `(a,b)` covered by a `v=nw` arc" into
endpoint-owner congruences. Packaging the two endpoint conditions symmetrically gives a
**single exact inequality per component**, a global **circuit** that ties them together,
and a **positivity** functional whose sign *is* loose-vs-tight.

## 1. The exact endpoint-cover criterion (midpoint form)

For a component `C_i=(a_i,b_i)` of `G(S')` (midpoint `m_i`, length `ℓ_i`):
```
C_i ⊆ D_v  ⟺  ∃ j∈ℤ: ‖v a_i − j‖ ≤ 1/n and ‖v b_i − j‖ ≤ 1/n   (both endpoints near ONE integer)
           ⟺  ‖v m_i‖ ≤ 1/n − (v/2) ℓ_i.
```
The left form is the literal **endpoint cover** (each endpoint's `v`-phase within `1/n`
of a common arc index `j`); the right form packages it into the midpoint's `v`-phase
versus a length-discount. **Verified 100%** against direct tight/loose computation
(`2492/2492`…`700/700`, n=6..14) — exact, not heuristic.

## 2. The circuit

`S` is tight ⟺ **every** component satisfies the criterion. These are not independent:
ordering the `M` components around the circle, their arc indices `j_i` are monotone and
**wind once**, `Σ_i (j_{i+1}-j_i) = v`. So a single integer `v` (the multiple of `n`)
must *simultaneously* bring **all** `M` midpoints within `1/n` (phase) of an arc centre
`j/v`. The cover is a **closed circuit** on the `v`-lattice, not `M` separate problems.

## 3. The positivity functional

```
P(S) := max_i ( ‖v m_i‖ + (v/2) ℓ_i − 1/n ).      P(S) > 0  ⟺  S loose;   P(S) ≤ 0 ⟺ S tight.
```
`P` is the per-component cover-deficit, maximised over the circuit. **Verified `P(S)>0`
for every multiple-of-`n` config** sampled (all loose). So C′ ⟺ `P(S) > 0` on the whole
multiple-of-`n` class — the conjecture is now a **single scalar positivity**.

## 4. A PROVED consequence (summing the circuit)

Summing the criterion `‖v m_i‖ + (v/2)ℓ_i ≤ 1/n` over the circuit:
```
S tight  ⟹  Σ_i ‖v m_i‖  ≤  M/n − (v/2)·μ(G(S'))  <  M/n,
```
so **the component midpoints must have average `v`-phase distance `< 1/n`.** This is
unconditional (just adds the `M` inequalities). Measured reality:

| n | mean of avg_i ‖v m_i‖ | tight needs `< 1/n` |
|---|---|---|
| 6 | 0.245 | 0.167 |
| 8 | 0.239 | 0.125 |
| 10 | 0.243 | 0.100 |
| 14 | **0.245** | **0.071** |

The midpoints' average `v`-phase sits at `≈0.245` — essentially the generic value `1/4`
of an unresonant point — while tightness needs it below `1/n`. **The positivity margin
grows toward the frontier** (`0.245` vs `0.071` at n=14). Tightness would demand `v` to
resonate anomalously with *all* component midpoints at once; generically it resonates
with none.

## 5. Why this is the right object (and how it unifies the prior criteria)

`P(S) > 0` is a clean restatement of "the periodic arc family misses some safe
component," now as one scalar per config. It subsumes the earlier proved criteria:
- **Criterion B′** (long component): if `(v/2)ℓ_i > 1/n` then the `i`-th term of `P` is
  `> ‖v m_i‖ ≥ 0`, so `P>0` — the long-interval case.
- **Lemma C** (both small owners): forces the component's endpoints onto exact arc
  centres, which can't coincide — the criterion's common-`j` condition fails.
Both are special ways `P>0`; the residual is exactly the configs where `P>0` holds only
through the **circuit's simultaneous-resonance infeasibility**, not a single component.

## 6. Honest status

- **Lemma D** (exact criterion, both forms): **PROVED**, verified 100% (n=6..14).
- **Circuit-positivity functional** `P(S)`, `P>0 ⟺ loose`: defined; verified `>0` on all
  multiple-of-`n` samples.
- **Summed corollary** (tight ⟹ avg midpoint `v`-phase `< 1/n`): **PROVED**; the real
  average `≈0.245 ≫ 1/n`, margin growing with `n`.
- **C′ ⟺ `P(S) > 0` always**: **OPEN** — the simultaneous-resonance / circuit Diophantine
  infeasibility (the `~1%` large-owner short-component residual of S574, now as one
  scalar).

**Artifacts:** `04-computation/lrc_endpoint_cover_circuit_positivity_s575.py` (+`.out`),
`lrc_circuit_positivity_avg_phase_s575.out`. Folded Lemma D + `P(S)` into THM-398 (§4¾).
Builds on S574/HYP-2105 (translator, Lemma C), THM-398/HYP-2104/2103/2102. New: **HYP-2108**.
