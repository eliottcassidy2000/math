# The next step: the residual is the narrow-class ℤ/2 of ℚ(√21), and exhibiting √21 means building the covering-min certificate as a left–right Cayley-complex cocycle (Annals LTC connection)

*opus-2026-07-01-S27. The owner asked to compute the next step in the ι-odd certificate in ℚ(√−3,√−7),
contemplate exhibiting √21 in the covering-min certificate, and weigh three external links. The next step is
arithmetically sharp — the residual is a **narrow-class ℤ/2** — and one of the links (Annals, good LTCs via
left–right Cayley complexes) is the exact cohomological machinery for exhibiting it.*

## The next step (verified): the residual is the narrow class group ℤ/2 of ℚ(√21)
Going one level past the Gauss sums (HYP-3818), compute the **units and genus theory** of the bridge field:
- `ℚ(√21)`: fundamental unit `(5+√21)/2`, **norm = +1** ⇒ the negative Pell `x²−21y²=−1` is unsolvable ⇒ the
  **narrow class number is `2·(wide)`**: a residual **ℤ/2** in the narrow class group.
- Genus theory confirms it: `disc = 21 = 3·7`, two ramified primes (both ≡3 mod 4), so the narrow 2-rank is
  `t−1 = 1` ⇒ `Cl⁺(ℚ(√21)) ⊇ ℤ/2`.
- The two ι-odd primes each carry the same: `ℚ(√3)` and `ℚ(√7)` have norm-**+1** units (ℤ/2 each); the
  covering's other prime **`ℚ(√61)` has norm −1** — narrow = wide, **no residual**. So `61` is the *clean*
  (ι-even, ≡1 mod 4) split factor, and the entire obstruction concentrates on the ι-odd `3·7 = 21`.

So the "next step" is not another Gauss sum but the **narrow class group**: the residual is the **ℤ/2 that genus
theory forces on `ℚ(√21)` from its two ι-odd ramified primes 3 and 7**. This is exactly the shape of the
OPEN-Q-108 odd index — a **ℤ/2 Borsuk–Ulam degree** (THM-582) — now realized as an ideal-class obstruction. The
narrow class group *is* the sign/orientation (ι-graded) refinement of the wide one, so **narrow-ℤ/2 = the ι-odd
degree**, and it lives on 3·7 and vanishes on 61 — matching precisely which primes are "hard."

Two more pieces of the ladder:
- **Hasse–Davenport lift:** over `𝔽_{p²}`, `g₂(χ_p) = −g(χ_p)² = −(i√p)² = +p`. The ι-odd `i√p` becomes the
  ι-even prime `p` one level up — the odd certificate is a "square root" of the even prime, and √21 is the cross
  square root `√(3·7)`.
- **Compositum:** `√−7 ∈ ℚ(ζ₁₄)` (the tight-AP heptagon), `√−3 ∈ ℚ(ζ₁₈₃)` (the covering), so `√21 ∈
  ℚ(ζ_{lcm(14,183)}) = ℚ(ζ_{2562})`. **√21 surfaces only when both the heptagon (7) and the Eisenstein covering
  (3) are present** — never in either regime alone.

## Exhibiting √21 in the covering-min certificate — the Cayley-complex route
The covering-min binds at `Φ₆ = 3·61`, and its atomic certificate (2 atoms `{t*,1−t*}`, `t*=14/183`) has moments
in `ℚ(ζ₁₈₃)` — carrying `√−3` (and the clean `√61`) but **not** `√−7`. The heptagon `√−7` enters from the LRC-14
constraint (`N=14=2·7`, the tight-AP structure). So **exhibiting √21 = writing a certificate that is
simultaneously valid for the covering (√−3) and the tight-AP heptagon (√−7)** — a certificate over the
biquadratic `ℚ(√−3,√−7)`, whose *residual obstruction* is the narrow-class `ℤ/2` of `ℚ(√21)`.

The concrete target: the certificate is a **2-cocycle** — the cup product of the two ι-odd 1-cocycles (`i√3`
Eisenstein, `i√7` heptagon) — and the obstruction to trivializing it is the narrow-`ℤ/2`. **This is a
left–right Cayley complex.** (Link 2 below.) Build the square complex on the group `ℤ/2×ℤ/2 = Gal(ℚ(√−3,√−7))`
with the two involutions (antipode `ι`, QR-negation) as the two generating sides; the covering-min certificate
is a **cocycle** on it, and `√21` (the narrow-ℤ/2) is the nontrivial cohomology class that must be exhibited.
The far-element/`χ(lonely set)` sharpening (S25) is where this cocycle should be read off.

## The three links (assessed)
1. **Cornell CS 6840 (Algorithmic Game Theory), Sept 12** — early AGT, i.e. **equilibrium existence via
   Brouwer's fixed-point theorem** (Nash), plus the SOS-for-games / spectral-certificate literature that
   neighbors it. *Relevance: the fixed-point + certificate toolkit I've been using* (Brouwer/Borsuk–Ulam/SOS,
   HYP-3814). Tangential but on-theme: Nash = Brouwer = my easy side; the hard side needs the topological degree.
2. **Annals 2026 (203-2, p.03): "Good Locally Testable Codes" (Dinur–Evra–Livne–Lubotzky–Mozes)** — resolves
   the c³-problem via **left–right Cayley complexes**. *Relevance: DEEP.* This is exactly the machinery of my
   certifying cohomology (S25) and the Cayley transform (S22): a Cayley square-complex whose **cohomology is a
   locally-testable code**, i.e. a *locally certifiable* certificate. The √21 residual = a cohomology class of
   such a complex; "locally testable" = the covering-min certificate checkable on local data. This is the
   right framework to **construct** the √21 cocycle — the strongest of the three leads.
3. **GitHub `Pengbinghui/pipeline-math`** — an AI prover–verifier pipeline (GPT + **Lean**) resolving open
   problems, including "tiling complements." *Relevance: methodological* (mirrors this repo's multi-agent
   prover/court process and codex's Lean formalizations, e.g. THM-346), with a thematic "tiling" tie. A model
   for formally verifying the covering-min certificate once constructed.

## Status
- **Verified (the next step):** `ℚ(√21)` fundamental unit `(5+√21)/2` norm **+1** ⇒ narrow-class **ℤ/2**
  residual; genus 2-rank 1 (disc 21=3·7); `ℚ(√3),ℚ(√7)` norm +1 (ℤ/2), `ℚ(√61)` norm −1 (clean); Hasse–Davenport
  `i√p ↦ p`; `√21 ∈ ℚ(ζ_2562)` (needs both 7 and 3).
- **Grounded:** the OPEN-Q-108 residual = the **narrow-class ℤ/2 of ℚ(√21)** (an ideal-class realization of the
  ℤ/2 Borsuk–Ulam odd index); it lives on 3·7, vanishes on 61.
- **Route to exhibit √21:** build the covering-min certificate as a **cocycle on the left–right Cayley complex**
  of `ℤ/2×ℤ/2` (the two involutions), with `√21` the nontrivial (locally-testable) class — the Annals LTC
  framework is the tool; `χ(lonely set)`/far-element is where to read it.
- **Honest:** the class-field arithmetic is exact and classical; "the covering-min certificate is a Cayley-complex
  cocycle with √21 the narrow-ℤ/2 class" is the sharpened conjecture/route, not yet a construction; the AGT link
  is tangential, the LTC link is the substantive one.

Related: HYP-3818 (the biquadratic bridge / √21), OPEN-Q-108/THM-582 (the ℤ/2 odd index), HYP-3817 (certifying
cohomology / chain complex), HYP-3814 (Brouwer/Borsuk–Ulam), HYP-3812/klein (Φ₆ metric-irreducible), kps-S21
(the √−7 vs √−3 split). External: Annals 203-2 p.03 (LTCs / left–right Cayley complexes — the cohomology tool);
CS6840 (Brouwer/Nash + SOS certificates); pipeline-math (Lean prover-verifier). HYP-3819 (this). Script:
04-computation/biquadratic_narrow_class_group_residual_opus_20260701.py.
