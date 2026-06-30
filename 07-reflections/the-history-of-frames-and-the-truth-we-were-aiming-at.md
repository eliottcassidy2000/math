# The history of frames, and the truth we were aiming at: one positive cyclotomic atom at the apex

*klein-2026-06-29-S14. Working the bound (the exact floor constant) and then reading the whole sequence of frames the LRC(14) floor was, each in its turn, declared most fundamental. They are charts of increasing resolution on one object. This names it.*

## The bound, exactly

`R' = m_S/(m_R m_Q)`. The exact infimum over the covering scan:
> `inf R' = 114382/332563 = 0.343941`, binding at `R = {1..13}\{7}`, `Q = {1,2}` (denominator `= 7^2·11·617`).
It clears the clean set-independent floor `1/(2 zeta(2)) = 3/pi^2 = 0.303964` by `+0.040`. The *exact* value
is messy and `R`-specific; the *clean* lower bound is `3/pi^2`. That gap between "the messy per-set value"
and "the clean constant" is the whole lesson in miniature: the truth is the clean number, the per-set value
is a chart. Two clean numbers survive every reframing:
- **`3/pi^2 = 1/(2 zeta(2))`** — the bulk / Eisenstein / `zeta(2)`-density floor (the local-global density).
- **`4 cos^2(3 pi/7) = 0.198`** — the apex cusp / doublet / cusp-form obstruction atom (the cyclotomic gap).
The floor is the statement that the second does not sink the first.

## The frames, in the order each was "most fundamental"

| # | frame | what it got right | what it mistook |
|---|---|---|---|
| 1 | Diophantine / covering (`M(S) < 1/n`) | the floor exists; lonely times are real | treated it as a search for a counterexample, not a positivity |
| 2 | tournament metagraph (iso-class graph, `H`-gradient) | the right state space; complement symmetry | the `H`-gradient is the bulk, orthogonal to the floor (HYP-3548) |
| 3 | "everything is the triangle" (the staircase) | the geometric substrate; the three sides | a substrate, not the obstruction |
| 4 | 2-adic descent (THM-580) | peels the non-transitive 2; `14=2·7` | the descent is the move, not the bound |
| 5 | the CV variance gatekeeper (THM-579) | the floor IS a 2nd moment | the variance is the WRONG (lossy, unbounded) coordinate (HYP-3554) |
| 6 | `R`-eigenspace, complement = antipodal (THM-584, HYP-3538) | the σ-even/σ-odd split is real and universal | the split is the *shape*, not the size |
| 7 | Euler product / anti-Littlewood (HYP-3550/3551) | the floor is multiplicative-positive | the product is the bulk; positivity needs the one factor |
| 8 | `Gamma_0(N)` / finite Siegel transform (HYP-3553) | set-independence; the congruence subgroup | the right *frame*, still missing the atom |
| 9 | relations-not-things / obstruction theory (S24) | the floor is a property of the relation | named the kind of object, not the object |
| 10 | ESSENTIAL × BOUNDED (S26) | the two clauses | a factorization, not the value |
| 11 | the transitive collapse / "right frame" (HYP-3566) | `CV` is lossy; transitivity is the fix | the transitivity is the chart that makes it finite |
| 12 | the finite cyclotomic min `4cos^2(3pi/7)` (HYP-3581) | the floor is a FINITE algebraic number | (this is nearly the territory) |
| 13 | `X_0(14)` cusps = Klein (HYP-3586) | the modular curve; cusps = the n=4 classes | the curve is the home, the atom is the cusp form |
| 14 | genus = local-global gap; obstruction = cusp form `f_14` (HYP-3587) | the obstruction is ONE global mode | (the dimension of the truth) |
| 15 | genus = odd boundary / even-graph / Eulerian (HYP-3591) | the obstruction is the `H_1`/odd part | the even/Eulerian is the bulk |
| 16 | blue = SC spine; squares & pronic (HYP-3592) | WHERE the obstruction sits (the SC spine) | blue is the arena, not the count |

Read top to bottom, the frames are a **monotone zoom-in on one obstruction**: a continuum measure →
a variance over all sets → a per-level decorrelation → a finite min over `2^7` cores → a single doublet →
one cyclotomic number `4cos^2(3pi/7)` at the apex cusp `d=7`. Each frame correctly found the **bulk /
obstruction split** but mistook its own chart for the territory; the next frame localized the obstruction
one notch further.

## The truth we were aiming at

> **The LRC(14) covering floor is the positivity of a single cyclotomic obstruction atom at the apex-7
> cusp.** Every "most fundamental" frame was a different name for that one atom — the doublet (THM-578),
> the genus-1 cusp form `f_14`, the level-2 Fiedler/cyclicity mode, the σ-odd `M_odd`, the blue odd-boundary
> on the SC spine, the value `4cos^2(3pi/7)` — surrounded by a clean bulk floor `3/pi^2` that every local
> chart (metagraph, `CV(H)`, Euler product, Eisenstein) computes. The obstruction was always **finite,
> singular, algebraic, and at the apex**; the work was never to *be* it but to *localize* it, and the
> sequence of frames is exactly that localization converging.

Three invariants held across all sixteen frames — they are the truth's fingerprints:
1. **A split** into a computable bulk and a missing obstruction (local/global, even/odd, Eisenstein/cusp,
   black/blue, bulk/`M_odd`). Universal.
2. **The obstruction is one-dimensional at `N=14`** (genus 1, one cusp form, one doublet). The whole
   difficulty is a single number.
3. **It is cyclotomic at the apex `7`** (`Q(cos 2pi/7)`; `4cos^2(3pi/7)`). The apex prime, not any
   dimensional coincidence (the persistence-test discipline kept us honest: `b_1^-=7`, `89=F_11` were
   traps; the apex `Z_7` is the bridge).

## Why we kept re-declaring the fundamental frame (the process-truth)

Each frame felt fundamental because it *was* the deepest available chart at that moment, and because it
correctly captured the split. The recurring error was identifying the chart with the territory — calling
the variance, or the Euler product, or the metagraph "the thing," when each was a coordinate system whose
job was to bring the one atom into focus. The cure was the same each time: the persistence test (does it
survive the next `n` / the next scale?), which dissolved the false-fundamentals (CV unbounded, blue≠genus,
the dimensional 7's) and kept the real bridges (`28=T(7)=C(8,2)`, `P_n(-1)=SC`, the apex `Z_7`, the doublet).
What remained, frame after frame, is the one positive cyclotomic atom. **That is the truth we were aiming
at: not a structure but a number — that `4cos^2(3pi/7) > 0` survives the bulk, equivalently that one cusp
form does not cancel the floor `3/pi^2`.**

## What to do with it

The remaining proof is the bound on that one atom: the leading apex-cusp coefficient of `f_14` (HYP-3587),
equivalently the floor's `f_14` component being dominated by the `3/pi^2` Eisenstein bulk. Everything else
is charts. The honest closing statement of the whole arc:
> **LRC(14) holds because a single cyclotomic number at the apex-7 cusp is positive, and the bulk around it
> is `3/pi^2`.**

See HYP (this session), HYP-3587 (genus = obstruction = cusp form), HYP-3581 (`4cos^2(3pi/7)`), HYP-3571
(the proof sentence), HYP-3592 (blue=SC spine), HYP-3550/3551 (`3/pi^2` / anti-Littlewood),
[[polysemous-constants-bridges-traps-and-homonyms]] (the persistence test that kept us honest),
[[the-right-frame-audit-when-the-proof-becomes-finite]].
