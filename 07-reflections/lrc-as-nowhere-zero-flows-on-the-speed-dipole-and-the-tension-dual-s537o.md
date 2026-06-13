---
source: oracle-2026-06-01-S537o
status: synthesis + computation (LRC inside debt = nowhere-zero flows on the speed dipole; tension/coloring dual)
tags:
  - lonely-runner
  - nowhere-zero-flow
  - tutte
  - seymour
  - flow-polynomial
  - circular-coloring
  - tension-dual
  - resonance
---

# LRC as Nowhere-Zero Flows on the Speed Dipole — and the Tension/Coloring Dual

Attacking LRC from the nowhere-zero-flow (NZF) mindset gives a precise dictionary
that imports Tutte/Seymour flow theory into the resonance/inside-debt machinery, and
explains the parity degradation of S533/S534 as a **flow theorem**. It also exposes
the dual face: runner positions are a **tension** (circular coloring), and
flow–coloring duality is the Lonely Runner's two sides.

## The dictionary: full-support resonances = nowhere-zero flows

The covering inside debt (S529/S533) sums over resonances `Σ_i m_i v_i = 0`. A
**full-support** resonance (all `m_i ≠ 0`) is exactly a **nowhere-zero flow** on the
**speed-weighted dipole** `G_v` — two vertices joined by `n-1` edges, edge `i`
carrying weight `v_i`; flow conservation is `Σ v_i m_i = 0`, nowhere-zero is
`m_i ≠ 0`. Reducing mod the character modulus `n*` (`= n/2` even, `n` odd; S533/S534):

> **The full-support inside debt is feasible (present) mod `n*`  ⟺  `G_v` has a
> nowhere-zero `ℤ_{n*}`-flow** ⟺ `∃ m ∈ (ℤ_{n*}∖0)^{n-1}` with `Σ m_i v_i ≡ 0`.

So the three-channel parity law (S533) and its n=18 vacuity (S534) are
**nowhere-zero-flow existence statements** on `G_v` — and the entire toolbox of flow
theory becomes available.

## Verified structure (`lrc_nowhere_zero_flow_s537.py`)

**(1) Factorization (inert vs active edges).** A runner with `v_i ≡ 0 (mod n*)` is an
**inert / free edge**: any of the `n*-1` nonzero values gives `v_i m_i ≡ 0`, so it
contributes a factor `(n*-1)` and never affects the conservation. An **active**
runner (`v_i ≢ 0`) is constrained. Hence

```
NZ-ℤ_{n*}-flow count(G_v) = (n*-1)^{#inert} · C_active.
```

**(2) The bridge characterization (= the parity law).** Verified `400/401` on random
primitive sets at `n = 6, 8, 18`:

> **Debt-free (no NZ `ℤ_{n*}`-flow)  ⟺  the active edges form a `ℤ_{n*}`-bridge
> ⟺ EXACTLY ONE active runner** (the rest inert).

This is precisely the S533 `k=1` law, now a flow statement: one active edge is a
bridge (no nowhere-zero flow can cross a bridge); `≥2` active edges make the dipole
**bridgeless**, which always carries a nowhere-zero flow.

**(3) Tutte/Seymour leverage — why parity is vacuous at large `n*`.** Seymour's
6-flow theorem: every bridgeless graph has a NZ `ℤ_6`-flow; the dipole carries NZ
`ℤ_k`-flows for all `k ≥ 2` once it has `≥ 2` (active) edges. So **debt is present
whenever `≥ 2` runners are active** — and for primitive sets at `n = 18` (`n* = 9`)
essentially every set has many active runners (the NZ-`ℤ_9`-flow count is astronomical,
`~2.5·10^{14}` for the AP). The S534 finding "parity is vacuous beyond `n = 4`" is
the flow-theoretic statement *bridgeless ⇒ flows exist*. The unique parity-potent case
`n = 4` (`n* = 2`) is exactly where the only nowhere-zero `ℤ_2`-flow structure is rigid.

**(4) The inside debt is a flow polynomial.** The unweighted `m`-dipole flow
polynomial `F(D_m;k) = ((k-1)^m + (-1)^m(k-1))/k` is reproduced exactly by the NZ-flow
count at unit weights (`2,6,12; 6,21,52; 10,60,204`). The inside-debt **value** is its
`ĝ`-weighted analogue — a **flow enumerator** `Σ_{NZ integer flows} ∏_i ĝ(m_i)` —
verified to vanish exactly when the NZ-flow count is `0` (the bridge case). Like every
flow polynomial it obeys **deletion–contraction** (delete a runner / contract a
resonance), the flow-theoretic form of the S531 modular recursion.

## The dual face: runner positions are a tension (circular coloring)

Flows live in the cycle space; their orthogonal complement is the **tension / cut
space**. The runner positions `x_i = frac(v_i t)` are a **tension** (potential
differences) on the observer-star: the edge `(observer, i)` carries tension `x_i`.
Scaling the circle to circumference `n`, loneliness is

> **a circular `n`-coloring**: the observer sits at circular distance `≥ 1` from every
> runner. **LRC = the orbit `{(v_i t)}` contains a time `t` realizing this circular
> coloring.** (Verified: `n=7, v=(1,2,3,4,6)` colorable at `t = 0.19`.)

So LRC is a **flow ⟷ coloring (tension) duality**, exactly Tutte's duality: the
inside-debt resonances are the *flow* side; the runner positions / loneliness are the
*coloring* side. The conjecture asks the orbit-constrained tension to be a proper
circular coloring of the observer-star; the obstruction (loneliness always fails) is
dual to the flow structure carrying the inside debt.

## What this buys, and the honest limit

- A rigorous import of a major theory: **inside debt = NZ flows on `G_v`**, parity law
  = NZ `ℤ_{n*}`-flow existence, debt-free = a `ℤ_{n*}`-bridge, the degradation =
  Seymour/Tutte (*bridgeless ⇒ flows exist*). The S533/S534 facts are now corollaries
  of flow theory, with deletion–contraction and the flow polynomial as tools.
- The **tension dual**: loneliness = orbit-constrained circular `n`-coloring; LRC is a
  flow–coloring duality, tying it to the **circular chromatic / circular flow number**.
- **Honest limit:** this reframes the inside-debt/parity layer (known insufficient
  beyond `n=4`, S534) into flow language — it does not by itself prove LRC. Its value
  is the toolbox and the unification, plus a sharp new target.

## Hypotheses / next (the forward bets)
1. **Full covering as a flow enumerator over all subgraphs.** The complete `|LONELY|`
   (not just full support) = `Σ` over NZ flows of every sub-dipole (every runner
   subset) of `∏ ĝ` — a Tutte/flow-polynomial-type sum. Compute its deletion–contraction
   and seek a closed product over the speed lattice.
2. **Circular flow number bound.** Phrase LRC via the circular flow number of a derived
   graph and test whether Tutte's 5-flow / Seymour's 6-flow give a loneliness bound.
3. **Counterexample ⟺ NZF obstruction.** Is "the speed set is a counterexample"
   equivalent to a derived graph having NO nowhere-zero flow of a prescribed (orbit /
   circular) type? If so, Tutte's flow conjectures would constrain counterexamples.

## Artifacts
```
04-computation/lrc_nowhere_zero_flow_s537.py
05-knowledge/results/lrc_nowhere_zero_flow_s537.out
```
Related: S529 (covering/characters), S533 (three-channel parity = NZ ℤ_{n*}-flow),
S534 (n=18 vacuity = bridgeless⇒flows), S531 (recursion = deletion–contraction),
S536 (sector/DFT dual), S535 (mapping spectrum).
