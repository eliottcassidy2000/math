# The dyadic resolvent lens: 120, 320, the OCF weights, and a Galois view of the LRC(14) crux

*kind-mendel-2026-06-27-S9. Owner: integrate incoming + past work, understand what remains, and let the
solvable-quintic resolvent (roots 2,−4,8,−16; coefficients 120, 320; root = ⁵√2−⁵√4+⁵√8−⁵√16) show a
creative synthesis. Result: 120 and 320 are repo invariants that are exactly the symmetric functions of the
dyadic OCF weights, giving a Galois/solvability lens that unifies the three live attacks on the one remaining
crux. Pulled main throughout (kps-S254: Move A reduced to a finite constant-chase; mac-mini-S69: binding = bounded core).*

## 1. What remains in LRC(14) — the integrated picture

The proof is the **THM-079 (H=21) template** (my S8): reduce-to-atom + a Moon/forcing step.
- **Move A (reduce to the bounded core = the tournament Mode-A peel).** R1 remove-large / R2 omit-prime / R3
  dilation. *Essentially done* per the team: mac-mini-S69 "far helps, binding = bounded core" (adding a far
  speed pushes the miss-PGF Lee-Yang zeros **outward**, ρ 1.49→2.0; multi-far is **not** the obstruction);
  kps-S254 "multi-far floor reduced to a **finite constant-chase**" (R'≥0.642, EH not needed). My S7 `lcm`
  family shows this half is irreducibly analytic but its content is now just constants.
- **Move B (bound the atom = the Moon/forcing step) — THE REMAINING CRUX (★).** The bounded covering core has
  `M > 1/14`, because the tight locus is `{AP, GW}` (both non-covering, optima at `t=a/14`), and covering
  forces a multiple of 14 onto the observer at every `t=a/14` (the apex-7 floor) ⇒ `M>1/14`. Verified (S8).
  The open content is the **bounded-core extremality**: "consec/AP is the unique tight minimizer," equivalently
  the tight-locus characterization, equivalently (kps-S31y) "over-cover ⟺ exact K₃ (= I(K₃,2)=7, forbidden)."

So **everything reduces to one statement (★)**, and it now has **three equivalent live forms**:
| view | object | status |
|---|---|---|
| analytic (my S7/S8) | three-gap (Steinhaus) rigidity: only APs are 3-gap-rigid ⇒ max bimodal coverage | finite/algebraic handle |
| Lee-Yang (team, S69/S254) | miss-PGF `P(z)=Σp_t z^t` (deg 6): floor via `min|z|` of its zeros (consec `min|z|=1.489`, verified) | finite constant-chase |
| Galois (this session) | the dyadic resolvent of the OCF; consec as a solvability branch-point | new lens, below |

## 2. Where 120 and 320 live in this repo (and the dyadic resolvent)

The owner's resolvent quartic `x⁴+10x³−120x²−320x+1024` has roots **{2,−4,8,−16} = {±2^d}** — the **OCF
dyadic weights** (since `H = I(Ω,2) = Σ_k α_k 2^k`, the conflict-graph independence polynomial at fugacity 2).
Its symmetric functions are exactly our invariants:
- `e₁ = −10`, `e₂ = −120`, `e₃ = +320`, `e₄ = 1024 = 2¹⁰`.
- **120 = max |ΔH| at n=7 = 5!** (the arc-flip metagraph; sequence `2,4,12,32,120`; verified n=5→12, n=6→32,
  n=8 ≥ 360 — so 320 is **not** max-ΔH(8)). 120 is also `|S₅|` = the icosahedral/quintic symmetry order.
- **320 = #distinct H-values at n=8** (320 of 330 odd values in [1,661]) **= #{β₃=1 tournaments at n=6}**
  (THM-110/226) **= the THM-155 denominator** `−3n(n²−1)(n²−9)/320`.
- **1024 = 2¹⁰ = 2^{1+2+3+4}** = product of the apex energies; it is the apex coefficient scale.

So the owner's "120 and 320" are two *genuinely different* tournament invariants (an arc-flip energy and a
homology/spectrum count) that the dyadic resolvent ties together as `e₂, e₃` of the OCF weights. That is the
coincidence worth following: **the OCF's dyadic weights `2^d` behave like resolvent roots, and our headline
small-n invariants are their symmetric functions.**

## 3. The Galois/solvability synthesis (the improved argument)

The quintic `x⁵+20x³+20x²+30x+10` is *solvable* precisely because its resolvent has the structured dyadic
roots, giving the closed form `⁵√2−⁵√4+⁵√8−⁵√16` (alternating dyadic radicals). The tournament analogue:

> **The OCF `H = Σα_k 2^k` is "solvable" in the same dyadic sense — every H is a finite alternating-sign
> combination of the dyadic weights — and the LRC(14) obstruction is the one *forbidden* dyadic value,
> `7 = I(K₃,2) = 1+2·3`.** (`14 = 2·7`, and the apex-7 Moon step says covering forces this forbidden
> K₃/H=7 structure.)

This reframes the crux (★) as a **solvability/branch-point** statement, exactly as codex-S271 hinted
("solvable stationary normal forms mark branch-point events" in the miss-PGF's quintic derivative `P'(z)`).
The concrete improvement it suggests, joining the three views:

**FTA-duality bridge (coefficients ↔ roots).** The team bounds the Lee-Yang floor `min|z|` of `P(z)`
*numerically* per config (the open constant-chase). But `P(z)`'s **coefficients** `p_t` are governed
*algebraically* by three-gap rigidity (my S7: the AP's `N_E` is exactly computable from its ≤3 gap-lengths).
So: **use the AP's three-gap-exact coefficients to pin its resolvent/symmetric functions, and bound `min|z|`
from the symmetric functions (Newton's inequalities / the resolvent), turning the numerical Lee-Yang floor
into an algebraic one** — and the extremality (★) becomes "the AP uniquely makes the resolvent solvable
(branch-point), maximizing the cover." This unifies analytic (three-gap), spectral (Lee-Yang), and algebraic
(resolvent) into one extremality, which is the cleanest the crux has looked.

Why this could finish (or sharply advance) (★): the constant-chase the team is stuck on is "uniformly bound
`min|z|` of the bounded-core miss-PGFs." Newton's inequalities relate `min|z|` to the coefficient ratios
`p_t`, which three-gap controls *exactly* at the AP and *monotonically* off it (the bimodality majorization,
S7). So the algebraic route replaces "chase the constant over all bounded configs" with "prove a coefficient
(Newton/Maclaurin) inequality, extremal at the three-gap AP" — a finite, structural target.

## 4. Verified this session (exact / sampled)
- consec miss-PGFs (exact `p_t`): `consec_8 = [481,359,250,156,119,75,30]/1470`, Lee-Yang `min|z| = 1.489`
  (matches mac-mini's 1.49); `consec_13` (n=14 AP) `min|z| = 1.654`. (`04-computation/lrc14_misspgf_resolvent_kindmendel.py`)
- resolvent `{2,−4,8,−16}`: `e₁..e₄ = −10,−120,320,1024` (so 120, 320 are its `e₂,e₃`).
- max |ΔH|: n=5→12, n=6→32, n=7→120 (=5!); n=8 ≥ 360 (so 320 ≠ max-ΔH(8); 320 = distinct-H(8)/β₃(6)).

## 5. Honest status & leads
The Galois/resolvent material is a **lens**, not yet a proof step — but it gives the one remaining crux (★) a
third, *algebraic* face and a concrete finish-route (the FTA/Newton-inequality bridge above). Net: LRC(14) =
[Move A: done mod finite constants, team] + [Move B (★): the bounded-core extremality, now attackable via
three-gap **and** Lee-Yang **and** the dyadic resolvent, all the same statement].
- **Lead 1 (the bridge):** prove `min|z|(P) ≥ ρ*` via Newton/Maclaurin inequalities on `p_t`, extremal at the
  three-gap AP — converts the team's numerical Lee-Yang constant-chase into an algebraic extremality (★).
- **Lead 2 (branch-point):** make codex-S271's "solvable `P'(z)` normal form marks branch points" precise:
  the AP should be the unique config where `P'` (the quintic) hits the solvable dyadic resolvent.
- **Lead 3 (icosahedral):** 120 = `|S₅|`; the quintic↔icosahedron (Klein) and the repo's icosahedral threads
  (five-as-bridge, icosahedral-fifteen) — is the bounded-core extremality an icosahedral/`A₅` rigidity?

→ HYP-3133 (new), HYP-2906/2898/2880 (mine), kps-S254/S31y, mac-mini-S69, THM-079/110/155/226, OPEN-Q-108.
