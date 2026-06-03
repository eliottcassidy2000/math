---
source: claudebox-2026-06-03-S611
status: SYNTHESIS — the n=14 proof architecture, with each component mapped to its formalized Lean
  module and the single remaining formal gap (the large-owner residual).
tags: [LRC, n14, proof-architecture, formalization, unit-clock, division-sieve, owner-congruence,
  large-owner-automaton, rank-1-two-block, roadmap]
---

# The n=14 proof, assembled — and what is already machine-checked

A long session cycling explore/investigate/formalize produced a clean reduction of LRC(14) and put
its skeleton into Lean. This maps the architecture to the formal pieces so the remaining work is
precisely a single gap.

## The architecture (top-down)

LRC(14): every primitive 13-speed config has gap `≥ 1/14`. The reduction tree:

```
LRC(14)
├─ no speed divisible by m, some m ∈ {2..14}   →  LONELY at the m-clock (gap ≥ 1/m ≥ 1/14).   [~92%]
│     ▸ FORMALIZED: UnitClock.lean (not_dvd_mul_of_coprime, lonely_at_unit_clock).
│     ▸ Why: at t=a/m, runner v is at the origin iff m∣va iff m∣v (gcd(a,m)=1).
│     ▸ The danger structure: gcd(j,n)=Σ_{d∣n}φ(d)[d∣j].  FORMALIZED: DangerBlocks.lean.
└─ RESIDUAL: a multiple of every m ∈ {2..14}  (the covering / multiple-of-14 / C′ configs)
   = THM-398's S = S' ∪ {14w}; tight iff G(S') fits one 14w-arc → endpoint-owner congruences.
   ├─ both owners of a component small (< 14)  →  LOOSE for every w (rigid pins ⇒ a=b).   [Lemma C]
   │     ▸ FORMALIZED: OwnerCongruence.lean (rigidity, lemmaC, lemmaC_no_fit).
   └─ a LARGE owner (≥ 14)  →  the bounded CRT residue automaton (S581); off-centre fits.
         ▸ resonance bound |D| ≤ u_b K_a + u_a K_b.  FORMALIZED: OwnerCongruence.lean (resonance_bound).
         ▸ the obstruction is the rank-1 two-block (the apex, mod-2), dissolved in the odd additive
           face (mod 2n−1 = 27 = 3³) by the pair-sum sieve (HYP-2150/HYP-2075).
         ▸ **THE ONE REMAINING GAP**: accept(owner-automaton) ∩ valid(config) = ∅ (tasks t-0040/41).
```

## What is machine-checked (math-lean, all sorry-free)

| component | module | content |
|---|---|---|
| danger block-diagonalization | `DangerBlocks.lean` | `gcd(j,n)=Σ_{d∣n,d∣j}φ(d)`; lonely-clock = unit |
| the unit/division sieve (trivial half) | `UnitClock.lean` | no multiple of `n` ⇒ lonely at `t=a/n` |
| the apex / 2-block / lift | `ApexCertificate.lean` | count `(q−1)(q−1−|S|)`, trichotomy, the lift |
| owner rigidity + Lemma C | `OwnerCongruence.lean` | `rigidity`, `lemmaC`, `lemmaC_no_fit`, `resonance_bound` |
| sum-free / fractal / fusion | `SumFree`, `FractalSumFree`, `Fusion` | the 3-term / circuit-free arithmetic |

So the **easy half of LRC(14) is formalized**, the **danger structure is formalized**, and the
**small-owner residual (Lemma C) is formalized**. The entire remaining difficulty is one statement:
the large-owner residual is empty — the CRT automaton's accept language never meets a valid config.

## Where the analytic tools sit (this session's locating)

The resonance-energy / circle-method bound (HYP-2155/2053) proves LRC for the median config through
**n=7** and dies at **n=8** (HYP-2165) — it is a small-`n` tool and is vacuous at n=14. So n=14 is
*entirely* construction (the reduction tree above), which is why the architecture is combinatorial
(sieve + owner congruences), not analytic. The Vitali boundary (HYP-2054) between measure and
construction is global, at n=8; n=14 is deep on the construction side.

## The single next step

Close `accept ∩ valid = ∅` for the large-owner automaton at n=14: characterize the valid
`(u_a,k_a,u_b,k_b)` that arise as real `G(S')` components, intersect with the automaton's accept set
(the off-centre fits), and show it is empty — equivalently, that the rank-1 two-block is always
cleared by an odd pair-sum. That one lemma, on top of the formalized skeleton, is LRC(14).

**Artifacts:** this session — HYP-2165, HYP-2170; math-lean `DangerBlocks.lean`, `UnitClock.lean`;
the architecture above. Builds on the whole recent arc (HYP-2115/S581, HYP-2150, HYP-2145, HYP-2063,
THM-398, HYP-2075).
