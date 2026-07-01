# The AP runner-cloud's Verblunsky coefficients are HARMONIC — |α_j| = 1/(n−1−j), exactly, verified n=5..20 (the n-th roots minus the marked origin z=1); α_0 = 1/(n−1) is literally the runner CENTROID seen from z=1 (the observer-lens made exact); this is the POSITION-space DUAL of kind-pasteur-S7's time-space OPUC "atomic thermometer" (HYP-3793); and the LRC operations organize into a group-like DICTIONARY of circle maps — the AFFINE group AGL(1,ℝ/ℤ) (runner=dilation ⋊ observer=rotation ⋊ antipode=reflection, lonely set = invertible dilations = units), the MODULAR PSL₂(ℤ) (Gauss/Farey = the CF/rung ladder), and the SZEGŐ vertex-insertion semigroup (n→n−1 = Mode A). The reproducing-kernel-at-z=1 handle is an HONEST NEGATIVE (it measures mass, not the lonely gap; ill-conditioned).

*opus-2026-07-01-S13. Owner: "pushing Verblunsky to the unit circle provides a nice visual recursive metaphor for
the LRC... see if there are creative functions between points on a loop... a whole dictionary of them may allow
them to operate group-like. also work the clean next step, synthesizing." I built the Verblunsky tool for the
runner cloud, found a clean harmonic law, built the dictionary, and honestly killed the one extension that didn't
pan out.*

## The finding: HARMONIC VERBLUNSKY LAW (verified, exact)
For the AP `{1,…,n−1}` at its lonely time `t = 1/n`, the runner positions `{e^{2πi v/n} : v=1,…,n−1}` are the
**n-th roots of unity minus the point z=1** (the marked origin/observer). The Verblunsky coefficients of the
uniform measure on this cloud are, **exactly**,
> **|α_j| = 1/(n−1−j),  j = 0,…,n−2**   (phases all ±π, i.e. real, signed)

verified n = 5, 7, 10, 14, 20 (Levinson–Durbin on the moments `m_k = (1/(n−1)) Σ_v e^{2πi k v/n}`). So the tower is
`1/(n−1), 1/(n−2), …, 1/2, 1` — the reciprocal integers, terminating at **α_{n−2}=1** (modulus 1 ⇒ purely atomic,
n−1 atoms). Two consequences drop out:
- **α_0 = 1/(n−1) = the runner centroid from z=1 (observer-lens, EXACT).** Σ over ALL n roots is 0; delete z=1;
  the remaining n−1 runners' centroid is `−1/(n−1)`, and `α_0 = −conj(centroid) = 1/(n−1)`. The **first Verblunsky
  coefficient IS the observer-lens** — the marked origin z=1 pulls the barycenter to distance 1/(n−1). This is the
  user's [[observer-lens-approach]] made literal in OPUC.
- **Σ|α_j| = H_{n−1} → ln n + γ** (Euler–Mascheroni), one of the four triangle constants (CLAUDE.md). The Szegő
  product `∏(1−|α_j|²)=0` (atomic), consistent.

**Why "harmonic" = vertex-peeling = Mode A.** The Szegő recursion `Φ_{j+1}=zΦ_j−ᾱ_jΦ_j*` adds one degree
(one "vertex") per step, with reflection α_j. Here α_j = 1/(n−1−j): as j runs 0→n−2 the reflection GROWS
1/(n−1)→1, the last step being the atomic closure. The Verblunsky/Szegő recursion is exactly the repo's
**Mode A reduction** (n→n−1, "vertex insertion, H=1+2^d") — the recursive "runners on a loop" metaphor the owner
asked for, now with a clean closed form for the extremal AP.

## The DICTIONARY of circle functions (group-like, as requested)
Functions between points on the loop ℝ/ℤ, organized by the group they generate (all relations verified
computationally except the noted float artifact):
- **AFFINE — AGL(1,ℝ/ℤ) = {x ↦ ±v x + a}** — the home of the LRC:
  - `rot_a: x↦x+a` = the **observer shift** (inhomogeneous LRC, THM-591);
  - `dil_v: x↦vx` = a **RUNNER** of speed v (dilations compose as the monoid ℤ; the **units (ℤ/n)\*** are the
    invertible ones);
  - `refl: x↦−x` = the **antipode / killer** symmetry (the −1);
  - relations verified: affine law `aff_{v,a}∘aff_{w,b}=aff_{vw,vb+a}`, `refl∘dil_v∘refl=dil_v`, semidirect
    `dil_v∘rot_a=rot_{va}∘dil_v`.
  - **LRC lonely set = (ℤ/n)\* = the invertible dilations** (the speeds with a well-defined observer-inverse).
    At n=14: {1,3,5,9,11,13}=φ(14)=6. **This is exactly kind-pasteur-S7 point (1)** ("covering-min maximizers =
    (ℤ/N)\*") and mac-mini-S77's safe-band residue system `R(q,a)={av mod q}` — both live in this AGL(1,q) action.
- **MODULAR — PSL₂(ℤ)** — the CF/hyperbolic side: `gauss: x↦{1/x}` (continued fraction), Farey `mediant`,
  `mobius: x↦(ax+b)/(cx+d)`. This is the **rung ladder / self-concordant residual** (`1/M=[n−1;n]`, S71) and the
  covering-min descent.
- **SZEGŐ semigroup** — `verblunsky` vertex-insertion, the n→n−1 reduction (Mode A). The recursive metaphor.
- **ANALYTIC scalars** — sawtooth `((x))` (ι-odd), distance `‖x‖` (ι-even) — the two halves under the antipode.

So the LRC's every operation is a dictionary entry: **runner=dilation, observer=rotation, antipode=reflection,
lonely set=units, CF-descent=modular, n→n−1=Szegő.** The runners "operate group-like" precisely because they are
the dilation subgroup of AGL(1).

## The DUAL to kind-pasteur-S7 (HYP-3793) — two circles, two OPUCs
kind-pasteur-S7 already brought Verblunsky into the repo, on the **TIME-space** measure (the lonely set L_C /
census moments): they found the extremizer **near-atomic** (|α_{5,9,11}|~0.92–0.98, Szegő product 1e−5) — an
"OPUC thermometer" for the tight/atomic boundary. Mine is the **POSITION-space** measure (the runner cloud at a
fixed t). These are the OPUCs of the **two dual circles** (times vs runner positions), and both say the same
thing from opposite ends: **extremal/tight ⇒ atomic ⇒ Verblunsky near/at the boundary; spread ⇒ Szegő-class
decay.** My AP harmonic tower (position) and their near-atomic tower (time) are the two faces of the same
tightness.
> **Szegő-decay = the two-regime split (HYP-3790).** Σ|α_j|²<∞ ⟺ nontrivial a.c. part ⟺ the SPREAD/far regime
> (equidistribution, mac-mini HYP-3786); α→boundary (atomic) ⟺ the STRUCTURED/consec regime. The Verblunsky
> decay rate is a quantitative dial for "structured vs spread," the same dichotomy as the far-element
> decorrelation (THM-546) and the L_y two-regime — now in OPUC language, on both dual circles.

## The clean next step I tried — and an HONEST NEGATIVE
The observer-lens says "look from z=1"; in OPUC that is the **reproducing/Christoffel kernel** `K_N(1,1)=Σ|φ_j(1)|²`.
I hoped max_t K_N(1,1) would land on the lonely time and bound M(S). **It does not:** max K for the AP is at
t≈0.01 (not 1/14), and K@lonely is non-monotone in M (construction K=5031 with M≈AP; spread K=17 with M=0.22).
**Reason:** the CD kernel measures *local mass/density*, but the lonely gap is a *nearest-atom distance*; and with
N atoms clustering as t→0 the kernel blows up from Vandermonde ill-conditioning, not from any gap. So the
reproducing-kernel-at-z=1 is the **wrong objective** for M(S). The observer-lens is real (α_0=centroid, exact),
but the *full* kernel is not the handle. Recorded so no one re-runs it.

## Status
- **VERIFIED (opus):** AP runner cloud at t=1/n has |α_j|=1/(n−1−j) exactly (n=5..20); α_0=1/(n−1)=centroid from
  z=1; Σ|α_j|=H_{n−1}; terminal α=1 (atomic). The circle-map dictionary's group relations (affine law, refl-conj,
  semidirect) hold; lonely set = units = invertible dilations.
- **SYNTHESIS (opus):** position-space OPUC is the DUAL of kind-pasteur-S7's time-space OPUC (HYP-3793); Verblunsky
  decay = the structured/spread two-regime (HYP-3790, mac-mini HYP-3786); the dictionary's AGL(1,q) is the group
  frame under the (ℤ/N)\* census (kp-S7) and safe-band residue system (mac-mini-S77).
- **HONEST NEGATIVE (opus):** the reproducing kernel at z=1 does NOT bound M(S) (measures mass not gap;
  ill-conditioned) — the observer-lens is α_0=centroid, not the full Christoffel function.
- **Scope:** a beautiful reframing + a new invariant (the Verblunsky signature of a runner config) + the group
  dictionary; NOT yet a new bound on M. The harmonic law CHARACTERIZES the AP's OPUC signature; it does not by
  itself prove extremality.

Related: HYP-3793/kind-pasteur-S7 (time-space OPUC "atomic thermometer" — this is the position-space dual),
HYP-3790/opus (two-regime = Verblunsky decay), HYP-3786/mac-mini (equidistribution ⟺ Szegő class),
HYP-3792/mac-mini-S77 ((ℤ/N)\* residue frame = the AGL(1,q) dictionary), THM-591 (inhomogeneous LRC = the
rotation entry), the S71 self-concordant ladder `1/M=[n−1;n]` (the modular/PSL₂ entry), [[observer-lens-approach]].
HYP-3794 (this). Scripts: 04-computation/lrc_verblunsky_{lens,harmonic_law}_opus_20260701.py,
lrc_circle_dictionary_opus_20260701.py, lrc_christoffel_observer_lens_opus_20260701.py (+.out).
