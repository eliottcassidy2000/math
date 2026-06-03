---
source: claudebox-2026-06-03-S614
status: REFLECTION — Collatz seen through the repo's 2-adic/doubling/resonance machinery; it is the
  multiplicative resonance twin of the Lonely Runner. First Collatz content; two Lean modules.
tags: [collatz, syracuse, LRC, resonance, 2-adic, doubling, HYP-2117, parity-vector, lemma-A-B,
  binary-signature, formalization]
---

# Collatz is the Lonely Runner wearing a multiplicative mask

**Prompt (human):** consider Collatz and how it is a similar question to work we have covered.

It is the same question. The repo spent months on the Lonely Runner and built a machine for one
shape of problem: a multiplicative structure pitted against the 2-adic, with a *resonance* as the
obstruction and a *conjecture* that the resonance is trivial except at the base. Collatz is that
shape exactly, and the machine transfers without bending.

## The same two sides

In the Lonely Runner the multiplicative side is the speed set `{v_i}` and the 2-adic side is the
doubling map `x↦2x mod n` (HYP-2117), degenerate at the even `n=2q` seam. In Collatz the
multiplicative side is the `×3` and the 2-adic side is the `÷2`; the shortcut map `T(n)=n/2` or
`(3n+1)/2` folds one halving into every step, and `v₂(3n+1)` is the per-step 2-adic height — the
dynamical version of S596's rigidity-height `v₂(n)`. That height is geometric, `P(k)=2^{−k}`,
`E[k]=2`, which is why a Collatz step drifts by `3/4 < 1`: the 2-adic side wins on average. The
Lonely Runner's "even `n` is the hard frontier" and Collatz's "the `÷2` beats the `×3`" are the same
statement about the 2-adic coordinate.

## The same resonance

The Lonely Runner obstruction is the additive resonance `Σ m_i v_i = 0` — the speeds conspiring —
and its weighted mass is the resonance energy (HYP-2155). Multiply a Syracuse cycle around its loop
and you get `2^K · ∏ n_i = ∏(3n_i+1)`, i.e. `2^K = 3^L · ∏(1+1/3n_i)` — a near-equality `2^K ≈ 3^L`,
the powers of 2 and 3 conspiring. That is the additive resonance's multiplicative twin: a cycle *is*
a resonance, and the conjecture (no nontrivial cycle) *is* circuit-freeness for the 2-adic/3-adic
problem. Only the trivial `L=1` resonance is feasible with small elements, exactly as the Lonely
Runner's only "tight" base is the arithmetic progression. I formalized the cycle equation
(`cycle_resonance`) and `2^K ≥ 3^L` — the first Collatz content in the repo's Lean, sitting on the
same `Fin`/`Finset` machinery as the LRC modules.

## The same Lemma A / Lemma B split, and the same signatures

The Lonely Runner splits into circuit-free (randomness, equidistribution, Lemma A — the bulk) and
3-term/AP (structure, the resonant core, Lemma B). Collatz splits the same way: a *balanced* parity
signature (odd-density below `log₃2 = 0.631`) contracts — that is Tao's "almost all orbits attain
almost bounded values," the Lemma-A randomness side — while the resonant cycles are the structured
core. And the bridge object is the same kind of thing the tournament work has used throughout: a
**binary signature**. The Lagarias parity vector is a bijection `[0,2^K) ↔ {0,1}^K`; the first `K`
parities are a 2-adic signature of `n mod 2^K`, and the shortcut map is a shift on them — precisely
HYP-2117's doubling dynamics, with the `+1` of `3n+1` playing the role of the IFS `T`-branch's `+1`.
I formalized the inductive heart of that bijection (`shortcut_mod_pow`: a step advances the 2-adic
scale by one). Even THM-004's binary-signature/descent calculus is the same genre.

## Why the masks differ

The Lonely Runner's resonance is *additive* (`Σ m_i v_i=0`) and its 2-adic seam is *static* (the
modulus `n`); Collatz's resonance is *multiplicative* (`2^K≈3^L`) and its 2-adic structure is
*dynamical* (the per-step `v₂`). That is the whole difference: Collatz is the Lonely Runner with the
additive group replaced by the multiplicative, and the static clock replaced by the iterated map.
The transcendence of `log₂3` (so `2^K≈3^L` only at deep convergents) is the multiplicative echo of
the Lonely Runner's "the resonant speeds are a thin core." Both conjectures say: the resonance is
trivial, and you reach the base — the lonely time, or the integer 1.

**Artifacts:** `04-computation/collatz_lrc_resonance_s614.py`; math-lean `Math/Collatz/Resonance.lean`,
`Parity.lean`; new **HYP-2175**. Builds on HYP-2117, HYP-2155, HYP-2140, HYP-2145, THM-004.
