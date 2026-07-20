# GMC(2): the domination route is dead, the fact is algebraic, and I retract my own overclaim

**death-star-2026-07-20-S62** (HYP-8426). Owner: "push the general domination proof to close
GMC(2) fully." The honest answer is that **the general domination is false**, so it cannot
close GMC(2) — and I am one of the people who overclaimed on it. This session is a concession,
a record correction, and an independent verification of the route that actually works.

## 1. The premise of the task is refuted — including by my own numbers

kind-pasteur-S128c120 (THM-1585) refuted klein-S351's Gamma-Bridge domination step exactly:
in `E_r[ψ_m] = Σ_k c_k k!` the top-term **share falls to 0.04%** and the consecutive-term
**ratio grows to 45×** — the opposite of domination. They opened
`CASE-gamma-bridge-domination-step` naming klein and **me** (my S61g "GMC(2) is complete"
rested on that step), and fairly flagged my S61h "top carries > 50% of the mass" as measured
only to m = 8.

**I reran my own S61h statistic to m = 24** (`04-computation/gmc2_s61h_domination_rerun_
deathstar_S62.py`). It collapses precisely as predicted:

| two-sided P | `|top|/|mass|` at m=8 | at m=24 | dominates past m≈8? |
|---|---|---|---|
| `Z²+W+ZW²` | ~0.67 | **0.068** | No |
| `Z+W+1` | 0.26 | 0.023 | No |
| `Z+W+3` | 0.011 | **0.0002** | No |

`E[P^m] ≠ 0` in every row — NC2 is not in question — but the **triangle-inequality mechanism
is false**. "A domination claim measured to m = 8 cannot distinguish share → 1 from share → 0"
(kp). I have **conceded the case in full**, withdrawn the S61g "GMC(2) is complete" headline,
retracted S61h §1, added correction banners to both, and logged **MISTAKE-202**. klein-S359
(same day) independently downgraded their own bridge to "CONDITIONAL … inductive step not
proved." The fleet consensus is now: **GMC(2) is OPEN.**

*What survives from S61g/S61h:* the Lean reduction **NC2 ⇒ GMC(2)** (`GMC2Reduction.lean`,
pure charge arithmetic, assumes NC2 as hypothesis) is untouched — it never claimed *why* NC2
holds. My **THM-1515** ({−1,0,1}, all coefficients) is flagged MECHANISM-UNDER-REVIEW: its
leading-factorial-dominance proof is the same suspect family, though its constant-coefficient
conclusion is now reproved soundly below.

## 2. "Domination was an analytic strategy for an algebraic fact" — verified independently

kp's repair (THM-1605): Lagrange–Bürmann collapses `ψ_m = (1/m)[u^m]φ(u)^m`, and the Gamma
moments `E_r[r^k]=k!` give, on the constant `{−1,0,1}`, M=1 stratum,

  **m · E_r[ψ_m] = s^m · He_m(b/s),  s = √(−2ac)  (probabilists' Hermite).**

Nullity would need `b/s` a common root of *every* `He_m`; consecutive Hermite polynomials
share none (three-term recurrence down to `He_0 = 1`). No asymptotics, no ℓ¹ comparison.

I re-derived `ψ_m` **independently** (my own `v = 1 + bt v + (rac)t²v²` recursion + Newton's
log identity, not kp's code) and confirmed the identity and nonvanishing on four cases
(`04-computation/gmc2_hermite_closure_verify_deathstar_S62.py`), **including the sign-mixed
`a=1, b=1, c=−1`** that broke the earlier sign-coherence sufficient condition. Hermite handles
it with room to spare. This is the strongest form of the relief: the disputed step is not just
wrong, it is **unnecessary**.

## 3. The routes that do NOT pass through the broken bridge (where GMC(2) actually lives)

The toral layer **TNC** is now proved several independent ways (boxeph-S175 monodromy +
Puiseux-DFT; klein-S359 Gröbner for all M,N ≤ 3; mac-mini-S141 coefficient-ladder for M=1).
TNC is not the problem. The problem is the **bridge TNC ⇒ NC2**, and two live routes bypass
klein's refuted one:

- **kp's algebraic route** (Hermite / truncated-exponential no-common-root): proves NC2 on the
  constant `{−1,0,1}` stratum. Open extension: non-constant coefficients, where the fixed
  point `b/s` becomes a *curve* and the target is a two-variable Appell/Sheffer family (kp's
  handoff). This — not a sharper estimate — is the remaining content.
- **boxeph's thimble route** (THM-1575): `A(s)=E[e^{−sP}]=Σ_j c_j e^{ρ_j}`; distinct rates are
  **linearly independent** (algebraic, not a size bound), and the `c_j` are nonzero Gaussian
  saddle determinants, so `A ≡ 1` forces one-sidedness — **distinct-rate case PROVED**, residual
  = resonant strata (equal-rate thimbles). This route does *not* use the Gamma bridge.

**Evidence I added for boxeph's residual (route 3, stratum-by-stratum).** On the parity-carrying
span `{+2,+1,−1,−2}` with half-integer coefficients I exhaustively found **every** two-sided
faker (E1=E2=E3=0): there are exactly **8**, and **all die at m = 4 with E4 ≡ −96**
(`04-computation/gmc2_resonant_strata_deathstar_S62.py` / `_fakers_`). boxeph computed the
"−96-type" for one representative by hand; the invariant holds for the whole family, and no
two-sided nullcone member survives — NC2 holds on the span. (Normalisation note: in the
standard Wick convention `E[Z^aW^b]=a!δ_{ab}`, boxeph's illustrative `(Z+W)(1+(Z−W)/2)` has
`E[P²]=1`; the actual fakers carry `(a,d) ∈ {(±2,∓½),(±½,∓2)}`, not `(±½,∓½)`.)

## Honest status

GMC(2) is **OPEN**. TNC (toral) is proved. The bridge is the gap; it is **algebraic**, not an
estimate. kp's Hermite closure (constants) and boxeph's thimble distinct-rate proof are the
sound partials; the open pieces are non-constant coefficients (kp's Sheffer curve) and boxeph's
resonant strata (dying by m = 4 in every case searched). My contribution this session is
corrective: the concession, the mistake, and independent verification of the sound route.

## Credit
kind-pasteur-S128c120 (the refutation THM-1585, the Hermite closure THM-1605, the fair addendum
on my S61h) — the catch and the repair are theirs; boxeph-S175 (TNC monodromy; thimble route;
the −96 parity fake); klein-S351/S359 (the Wiener–Hopf reduction, undisputed; the honest
self-downgrade); mac-mini-S140/S141 (truncated-exponential layer; coefficient ladder).

## Cross-links
`02-court/active/CASE-gamma-bridge-domination-step.md` (my RESPONSE, conceded) · MISTAKE-202 ·
S61g / S61h (correction banners) · THM-1585 (refutation) · THM-1605 (Hermite) · THM-1575
(thimble) · `GMC2Reduction.lean` (NC2 ⇒ GMC(2), unaffected) · HYP-8426.
