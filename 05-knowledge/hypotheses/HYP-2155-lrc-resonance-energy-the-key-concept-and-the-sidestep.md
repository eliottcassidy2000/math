---
id: HYP-2155
status: SYNTHESIS + VERIFICATION — the resonance energy E(v) is the key LRC concept; it is the
  dynamical-face (measure) obstruction, large & blind on the core; the SIDESTEP is the additive-face
  construction (sieve), partitioned by the Vitali boundary. Verified; unifies HYP-2053/2054 with the
  two-faces frame.
source: claudebox-2026-06-03-S605
related:
  - HYP-2053
  - HYP-2054
  - HYP-2150
  - HYP-2120
  - HYP-2145
  - THM-369
---

# HYP-2155: the resonance energy is the key concept — and you sidestep it

## The key concept (HYP-2053, oracle-S550)

The covering identity `|LONELY(v)| = main + Σ_{0≠m: Σ m_i v_i=0} ∏ ĝ(m_i)`, `main = (1−2/n)^{n−1}`,
gives the **resonance energy**

```
E(v) = Σ_{0≠m, resonance} ∏ |ĝ(m_i)|,        |LONELY(v)| ≥ main − E(v),
```

so `E(v) < main ⇒ LRC(v)`. The "resonances" are the integer relations `Σ m_i v_i = 0` — the
runner frequencies summing to zero — and `E(v)` is their `ĝ`-weighted mass. **This `E(v)` is exactly
the S585 relation-lattice theta in L¹ form**: the same `Σ_{m∈Λ} ∏ ĝ(m_i)` whose `m=0` term is the
independence baseline `main`, here taken in absolute value as the worst-case correction.

## What it sees and what it is blind to (verified, `lrc_resonance_energy_sidestep_s605.py`)

- **Bulk:** generic / circuit-free / translated sets have `E(v) < main` (e.g. `0.07, 0.10, 0.12`
  vs `main ≈ 0.13`) — **LRC proven by the energy bound**. (Circuit-free sets have no length-3 term;
  their `E` starts at length 4 — the S585/HYP-2120 grading.)
- **High-energy core:** the AP / regular polygon has `E(v) ≫ main` (`{1..5}: 0.36`, `{1..6}: 0.59`),
  carried by the **length-3 (3-term fusion `v_a+v_b=v_c`) resonances** (the largest single block).
  The bound fails, and (HYP-2054, the Vitali wall) `measure(LONELY)=0` exactly — the measure/energy
  is intrinsically **blind** to the core.

## The sidestep (verified)

On the core you cannot push the resonance energy below `main` — so you **sidestep it**. Abandon
measure; exhibit the witness by **construction**. Verified: every high-energy core is lonely at the
rational sieve `t = 1/n`, with gap exactly `1/n` (the n-gon vertex, a measure-0 witness invisible to
the energy). `{1..5}→t=1/6`, `{1..6}→t=1/7`, Fibonacci→`t=1/6`, all tight.

## The synthesis: the resonance energy is the dynamical face; the sidestep is the additive face

The resonance energy unifies the whole recent arc as **the dynamical-face obstruction**:

- `E(v) = Σ` over the **relation lattice = the resonances = the multiplicative/doubling structure**
  (HYP-2150's dynamical face). It is a **measure** quantity; it is dominated by the **length-3
  fusions**, which are the 3-term relations whose 2-block carries the apex (HYP-2145/HYP-2063). It is
  large exactly on the worry-set, and blind to it.
- The **sidestep is the additive face** (HYP-2150): the rational sieve `t=a/n` and the mod-`2n−1`
  transversal rigidity (the 64 self-converse classes, all lonely). It is **construction**, not
  measure; it is apex-free; it handles the core that the energy cannot.

> **The LRC architecture, one sentence.** LRC = bound the resonance energy on the bulk (the
> dynamical face / measure, where `E < main`) ⊔ **sidestep** it on the core (the additive face /
> construction, the sieve `t=a/n`), partitioned by the **Vitali boundary** (HYP-2054) where measure
> goes blind. The resonance energy is the key concept; the sidestep is the only way past it.

This says precisely why every prior attempt that tried to *bound* the resonance energy on the core
stalled (S550's wall): the core's `E` is irreducibly large because its length-3 fusions are
maximal, and no measure bound sees the measure-0 witness. The work is not to bound `E` but to
**route around it** — and the route is the additive-face construction the human has been pointing at
(the 64 classes mod `2n−1`, the sieve).

## Open / next

- Quantify the sidestep: is the core (where `E ≥ main`) exactly the set with a short minimal
  resonance length `l(v)=min Σ|m_i|` (the AP-led small-`l` family)? If so, the core is a *finite*
  check per `n` (bounded `l`), and the sieve `t=a/n` (THM-369) handles the AP extremal.
- Bound the resonance-energy tail `Σ_{l>L} ∏|ĝ| ≤ C r^L` rigorously (HYP-2053's open tail) — clean
  on the bulk; combine with the construction on the core.
- Map the dominant length-3 fusions to the 2-block (HYP-2145): the resonance energy's hardest part
  is the apex 2-block, the same point the additive face dissolves.

**Artifacts:** `04-computation/lrc_resonance_energy_sidestep_s605.py` (+`.out`),
`07-reflections/lrc-resonance-energy-and-the-sidestep-s605.md`. Unifies HYP-2053 (resonance energy),
HYP-2054 (Vitali), HYP-2150 (two faces), HYP-2120/HYP-2145 (relation-lattice/blocks), THM-369 (sieve).
