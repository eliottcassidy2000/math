# Consolidated structure analysis: one object, four faces (S614–S617)

*claudebox, 2026-06-03. A consolidation of the recent arc — Collatz↔LRC (HYP-2175, S614), the iterated-log
altitude principle (HYP-2180, S615), the apex-lift certificate sheaf (HYP-2185, S616), and the covering-depth
master object with the additive-chain collapse family (HYP-2195, S617). The claim: these are four readings of a
single structure.*

## The one object

Across all four sessions the recurring object is the same: **a multiplicative/arithmetic structure tested against
the 2-adic clock, with a resonance as the only obstruction, and an iterated-log altitude counting how deep the
obstruction is nested.** Spelled out as a single template:

```
  STRUCTURE          ×  CLOCK            obstruction = RESONANCE         conjecture
  -------------------------------------------------------------------------------------
  speeds {v_i}          t ∈ ℝ/ℤ          Σ m_i v_i = 0  (small m)        gap ≥ 1/(n+1)   (LRC)
  ×3                    ÷2 (2-adic)       2^K = 3^L·∏(1+1/3n_i)           no nontrivial cycle (Collatz)
```

LRC and Collatz are the *same* template (S614, HYP-2175): a multiplicative object vs the 2-adic, resonance the
obstruction, the conjecture asserting the resonance is trivial except at the base (the lonely time / the integer 1).
Everything below is the internal anatomy of that template, assembled from the four sessions.

## The four faces of the obstruction

The obstruction — "is there a nontrivial resonance" — has four equivalent presentations, one per session. They are
the same statement; the table is the consolidation.

| Face | Object | "no obstruction" reads as | Session |
|---|---|---|---|
| **Measure** | covering-depth `p_k`, `p₀ = meas{lonely}` | `p₀ > 0` | S617 / HYP-2195 |
| **Sheaf** | apex-lift certificate sheaf, `H⁰` | `H⁰ ≠ ∅` | S616 / HYP-2185 |
| **Resonance** | relation lattice `Σ m_i v_i = 0` | only trivial relations active | S614 / HYP-2175, HYP-2155 |
| **Altitude** | iterated-log depth of the certificate | depth `= 1` (single averaging) | S615 / HYP-2180 |

The bridges, each verified or formalized this arc:

- **Measure ⟺ Sheaf.** `p₀ > 0 ⟺ H⁰ ≠ ∅`. The depth distribution is the *measure-theoretic refinement* of the
  sheaf's global sections: `H⁰` sees nonempty-or-not, `p₀` sees the measure. (S617 §open; S616 built `H⁰` as
  `CertLocus`, the intersection-of-arc-complements — literally `{depth = 0}`.)
- **Measure ⟺ Resonance.** The collapse family `{p₀ = 0}` is *exactly* the additive chains — sets closed under
  `a+b=c`, i.e. rich in `(1,1,−1)` resonances (S617, exhaustive small search; AP **and** sporadic `(1,3,4,7)`,
  `(1,3,4,5,9)`). Mechanism: clock-distance subadditivity `‖(a+b)t‖ ≤ ‖at‖+‖bt‖` pins the third arc (formalized
  `dZ_chain_le`). So "collapse" = "resonance-rich" is one sentence.
- **Resonance ⟺ Sheaf.** A resonance is a cross-lane dependence; in the sheaf it is the failure of the local
  certificates to glue. The apex (the multiple-of-n, σ-fixed lane) is the maximal resonance — it forbids the whole
  plane, emptying `H⁰`; the `r/p` lift restores codimension 1 (formalized `certLocus_eq_empty_of_apex`,
  `certLocus_apex_lift_nonempty`).
- **Altitude over everything.** "How nested is the obstruction" is the iterated-log depth (S615). #logs = #averaging
  levels = #independent relations broken = wall codimension. Verified in the covering picture: the singleton wall
  `{1,2}` opens as `p₀ ~ ε^{0.99}` — `(loglog)¹`, one level; multi-relation walls bend sub-linear.

## The 2-adic seam is the load-bearing wall

Every face localizes its obstruction to the **same place**: even `n = 2q`.

- Dynamical (doubling `x→2x mod n`): collapses to a rank-1 2-block, the apex/`⟨−1⟩`-fixed point `n/2` (HYP-2150).
- Sheaf: the apex is the unique **σ-fixed lane** of the antipodal involution `v↦−v` — the ramification point of the
  double cover (S616). The obstruction sits exactly at the fixed point.
- Resonance: the danger block-diagonalizes; rigidity "leaks" only through the `gcd(·,n)>1` characters — the 2-block
  (HYP-2145).
- Collatz: the entire problem is `2^a` vs `3^b`; `v₂(3n+1)` is the per-step rigidity-height (S614/S596).

This is the perspective key in its consolidated form: **loneliness/contraction fails to glue precisely at the fixed
locus of the antipodal automorphism.** Not chirality — automorphism rigidity, localized to the ramification lane.
The two faces of even `n` (the dynamical apex and the additive collapse) are one seam seen twice.

## What is now formalized (math-lean, all sorry-free)

```
  Collatz/Resonance.lean      cycle resonance 2^K∏n=∏(3n+1), 2^K ≥ 3^L           (S614)
  Collatz/Parity.lean         shortcut map, 2-adic shift (Lagarias bijection)     (S614)
  IteratedLog/Altitude.lean   geometric-descent engine a_i ≤ ρ^i a_0 + C/(1−ρ)    (S615)
  LonelyRunner/ApexCertificate.lean §Sheaf
                              CertLocus, gluing, apex obstruction, lift           (S616)
  LonelyRunner/CoveringDepth.lean
                              dZ subadditivity, additive-chain link, 1st-moment    (S617)
```

## The consolidated open problem

One target now subsumes the loose ends of all four sessions:

> **Characterize the collapse locus `{p₀ = 0}` and prove it equals the additive-chain / resonance family; show its
> only points are tight (gap exactly `δ`), never sub-`δ`.**

Because the four faces are bridged, a proof on any face transfers: the additive-chain ⟺ collapse characterization
(measure face) is the LRC tightness statement (resonance face) is the `H⁰` non-vanishing off the apex (sheaf face),
and its difficulty is graded by the altitude (singleton walls — rank 1, `(loglog)¹` — first; the apex/`2q` seam —
the genuinely nested level — last). For `n = 14` this is concretely two tasks (S616): the transverse lanes glue
nonempty, and the single lifted apex stalk stays nonempty. The covering-depth picture adds the quantitative target:
`p₀(14) > 0` with the additive-chain configs as the measure-zero boundary.
