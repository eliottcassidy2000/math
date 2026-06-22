# Creative swings at all three LRC(14) nodes (post-main-sync) — including the dead ends

*kind-mendel-2026-06-22-S4. The owner asked: pull the latest from main, use it as inspiration to
attempt all 3 open nodes, be creative/outside-the-box, and report ALL ideas tried — even the ones
that didn't pan out. This is that honest ledger. Standout: Swing 3 (bounded-denominator witness).*

## What main looked like after the sync (the inspiration)

Rebasing onto `origin/main` pulled in a burst of work that built directly on my S1–S3 (decorrelation
unification, sheet-counting). Key state:
- **HYP-2851 is the team's own retraction** of the kps-S32 "VERIFIED-CLOSED" overclaim ("HYP-2844/2845
  closure OVERCLAIMED; the proved-δ bound = 0 at the worst P={2,3,4,5,6}"). **My S1 audit was vindicated
  and acted on.**
- **Node 1** boundary core: closed for the pure dilation (my sheet-counting + mac-mini's discretization
  lemma `ρ_K ≥ ρ*−arcCount/V` + codex's mod-13 certificate). Three-gap lemma (HYP-2853) is the converter.
- **Node 3**: a lovely claimed **ζ(2) floor** `c_q ≥ 3/π² = 1/(2ζ(2)) ≈ 0.304` (HYP-2856, via rate-V
  Farey neighborhoods + Mertens), plus a spectrum-intersection floor `R' ≥ 0.53` via L2 Cauchy–Schwarz /
  Abel signed tail (HYP-2857/2861).
- **HYP-2862**: Nodes 2 & 3 are **one** backbone — L2-Parseval **√-cancellation** with a shared
  generic(Weyl/Erdős–Turán) + resonant(Farey) split.

Genuine remaining gaps after all that: the tail bounds are **verified numerically at fixed H, not
uniformly**; the ζ(2) floor is **un-restricted** (no `G_P`); Node-1 general (non-AP) clusters.

---

## Swing 1 — verify/probe the ζ(2) floor and the G_P-restriction gap  *(PARTIAL — de-risked, not closed)*

Verified the un-restricted floor `c_q = meas{maxgap{frac(jx):j≤2q-2} > 1/q}` ≈ 0.44 for q=3,5,7,9
(my grid; the team's "~0.55" differs slightly — convention/grid — but both sit comfortably above
`3/π²=0.304`). The Mertens/Farey lower bound to `3/π²` is sound in spirit.

**The probe that mattered:** does the floor survive intersecting with `G_P` (the gap the team flagged
as "modulo the general-cluster G_P version")? **Yes, robustly.** `meas(GOOD_E ∩ G_P)` came out
`0.26–0.65` across configs — all **4.7–11.5× m_P** — *including the very `P={2,3,4,5,6}` that zeroed
the team's conservative δ-bound* (HYP-2851). So the floor's **truth** is not in danger from `G_P`; the
only gap is the **uniform proof**. Outcome: didn't close it, but removed the scariest doubt and
confirmed the ζ(2)/Farey mechanism is attacking a real, robust quantity. Files:
`lrc14_creative_swings_kindmendel.py`.

## Swing 2 — Beurling–Selberg extremal functions for the floor  *(TOOL, not a closure)*

Idea (outside the team's L2-CS/Abel): replace `1_{G_P}` by its **Selberg minorant** `m` — a trig
polynomial of degree ≤ N with `m ≤ 1_{G_P}` and *exact* defect `∫(1_{G_P}−m) = r_P/(N+1)`
(`r_P`=#arcs). Then
> `meas(G_P ∩ cover^c) ≥ ∫ m·1_{cover^c} = Σ_{|n|≤N} m̂(n)·conj(\hat{1_{cover^c}}(n))`,
a **finite, rigorous** lower bound (band-limited, one-sided). Cleaner than the Abel/Cauchy–Schwarz
tail control because the sign is built in and the error is exactly `r_P/(N+1)`.

**Why it didn't close anything:** it yields a rigorous *per-config* finite certificate, but not a
*uniform* one — you still must bound the finite sum below over all configs, which is the same
uniformity wall. It is essentially power-equivalent to the team's tail bounds (all are rigorous
Fourier truncations), just packaged more cleanly. Verdict: a clean rigorous packaging worth keeping
for the eventual Lean certificate, but not a new closure.

## Swing 3 — the BOUNDED-DENOMINATOR WITNESS  *(STANDOUT — a possible elementary route)*

**Reframe the hard case constructively.** By THM-523 a counterexample must be a *covering* 13-set
(contains a multiple of every `q∈{2,…,14}`). Ask: does every covering set have an **explicit
small-denominator** lonely point `τ = a/D`, i.e. `‖s·a/D‖ ≥ 1/14` for all `s∈S`?

**Empirically, yes, with a small absolute bound:**
- The **loosest known** covering set `{1,…,11,13,84}` (M=7/89, the binding case of uniform looseness):
  witness `τ = 17/41` (**D=41**).
- **Random covering sets stay bounded regardless of speed size:** max `D` was 25 / 25 / 23 / 26 for
  speeds up to `10³, 10⁴, 10⁵, 10⁶` (mean ≈ 18). **D does not grow with magnitude.**

**Why magnitude is irrelevant (the structural insight):** a witness `τ=a/D` depends *only* on the
residues `{s mod D}` — a finite amount of information. So the existence of a good `a/D` is a question
about residue patterns mod `D`, completely independent of how large the speeds are; the covering
condition constrains those residues. This is exactly why `D` is bounded even at `10⁶`.

**Status / reduction.** This is the conjecture
> **(HYP-2864) every primitive covering 13-set has `τ=a/D` with `D ≤ D₀` (empirically `D₀=41`).**
It would prove LRC(14)'s hard case **with zero analysis** (no Part-A, no cap, no decorrelation). It is
a sharpening of the repo's HYP-2566 (uniform looseness: covering sets have `M ≥ 7/89 > 1/14`): bounded
`D` ⟸ uniform looseness, but **not** via naive Dirichlet — Lipschitz/Dirichlet only gives
`D ≤ max_speed/(2ε)` (unbounded), whereas the truth is bounded. So bounded-`D` is a genuine
**arithmetic** fact (residue-pattern finiteness), stronger than the metric statement, and the right
target is: *classify the residue patterns mod `D≤D₀` of covering sets and show one always admits a good
multiplier `a`.* Files: `lrc14_bounded_witness_conjecture_kindmendel.py`, `..._crt_witness_stress_...py`.

---

## Ideas considered but not developed (honest)

- **Character sums / Fano-plane over F₇** for the decorrelation (kps-S34 already saw `q=7` is
  QR(7)/Fano-balanced). Plausibly the cleanest route to the *apex-prime* decorrelation, but I did not
  develop a Gauss-sum bound.
- **Second-moment / Paley–Zygmund** for nonemptiness of the lonely set — the first moment *is* the
  floor; the variance would need the same correlation control, so likely circular.
- **Direct LRC(14)→LRC(7) via the `14=2·7` doubling** (Cayley–Dickson / double-round-robin THM-378):
  appealing but I found no actual reduction; speculative.
- **Geometry-of-numbers / Minkowski** lonely-lattice-point: overlaps the existing covolume/HYP-2606
  singular-series line; no new traction found.

## Honest meta-conclusion

The team's recent ζ(2)/√-cancellation/three-gap work is real progress and Swing 1 confirms the witness
floor is **robustly true** (G_P does not kill it). But the **uniform rigor** of the analytic route
(Nodes 2 & 3) is still open — the tail bounds are per-config/numerical, and Beurling–Selberg (Swing 2)
packages but does not transcend that. The most promising *new* direction is **Swing 3**: the hard
(covering) case may admit an **elementary, magnitude-independent, bounded-denominator witness**,
sidestepping the entire analytic apparatus. It reduces to a finite residue-pattern statement and
refines HYP-2566. If I were to spend the next session anywhere, it is here: try to *prove* `D ≤ D₀` by
classifying covering residue patterns mod the small `D`'s that actually occur (14, 16, 41, …).

→ HYP-2864 (new), HYP-2847/2848 (mine), HYP-2851/2856/2857/2861/2862/2566 (team), THM-523, OPEN-Q-108.
