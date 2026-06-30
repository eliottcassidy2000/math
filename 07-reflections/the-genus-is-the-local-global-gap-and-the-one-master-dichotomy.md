# The genus is the local–global gap: the bulk is constant, the obstruction is the cusp form, and every dichotomy is one

*klein-2026-06-29-S11. Working the concrete next step (decompose the Gamma_0(14) 2nd moment into Eisenstein + f_14) and asking what the genus shift fundamentally MEANS. The answer reorganizes the whole project into a single dichotomy.*

## What the genus fundamentally represents

Weight-2 modular forms on `Gamma_0(N)` split `M_2 = Eisenstein ⊕ S_2(cusp forms)`, with
```
dim Eisenstein_2 = (#cusps) − 1      [the BULK / boundary space]
dim S_2          = genus(X_0(N))     [the OBSTRUCTION / global space].
```
Across the LRC(2p) family the **bulk is CONSTANT**: every `X_0(2p)` has 4 cusps, so `dim Eisenstein = 3`,
always. Only the **obstruction grows**: `dim S_2 = genus = 0,0,1,2,2` for `p=3,5,7,11,13`. So:

> **The genus is the dimension of the space of global modes that the boundary (cusp) data does not
> determine.** A weight-2 form is fixed by its values at the cusps *up to* the genus-dimensional freedom of
> cusp forms (which vanish at every cusp). The genus is exactly the kernel of "form ↦ its cusp values" —
> the **local–global gap**.

This is what the `0→1` shift at `N=14` means, made precise. Below the shift (genus 0, LRC(6)): the floor is
determined by boundary/cusp data alone — a local computation, the Hasse principle, the Euler-product floor
(HYP-3550), the metagraph/`CV(H)` bulk rehearsal SUFFICES. At the shift (genus 1, LRC(14)): there is exactly
**one** global mode — the cusp form `f_14` (the rank-0 elliptic curve `14a`) — that no boundary computation
sees. That single mode is the obstruction. "The testbed models the bulk, not the cusp" (klein-S4) is exactly
"the rehearsal is Eisenstein; the missing piece is the 1-dimensional cusp form."

## The concrete next step, worked

"Decompose the `Gamma_0(14)` 2nd moment into Eisenstein + `f_14` and bound the cusp-form piece at the apex
cusp." Three concrete facts:

1. **The decomposition is `3 + 1`:** three Eisenstein series (one per cusp, minus the constant) plus one
   cusp form `f_14`. The Eisenstein part is the bulk the rehearsal already controls; the work is the single
   `f_14`.
2. **Cusp forms VANISH at the cusps.** So at the apex cusp `d=7` itself, the `f_14` contribution is `0` —
   the *value* of the 2nd moment at the cusp is pure Eisenstein (locally computable). The obstruction is the
   **rate of vanishing**: the leading coefficient of `f_14`'s `q`-expansion at the `d=7` cusp. "Bound the
   cusp-form piece at the apex cusp" = bound that one leading coefficient — a finite, specific number, not a
   function. (This sharpens mac-mini S29's "the proof lives at the cusp": the binding is the *approach* to
   the cusp, where Eisenstein gives the leading term and `f_14` the subleading obstruction.)
3. **Rank 0 makes it non-degenerate:** `L(f_14,1) ≠ 0`, so the obstruction is a fixed nonzero thing — the
   floor is bounded away from `0`, not marginally.

So the remaining work has collapsed to: bound the leading apex-cusp coefficient of a single, explicit,
rank-0 weight-2 newform. That is the entire `f_14` content of the floor.

## The one master dichotomy

Every "two-index" split this project has found is the SAME split, now named on the automorphic side:

| LOCAL / bulk / determined | GLOBAL / obstruction / missing |
|---|---|
| Eisenstein series | cusp form `f_14` |
| boundary / cusp data | interior / period data |
| Hasse principle / local-global holds | the local-global gap (genus) |
| `sigma`-even / R-even (THM-584) | `sigma`-odd / R-odd |
| Brouwer / SOS-provable bulk | Borsuk–Ulam / the signed obstruction |
| metagraph / `CV(H)` / transitive rehearsal | the cusp-form mode the rehearsal misses |
| Euler product floor (HYP-3550) | anti-Littlewood positivity (HYP-3551) |
| the `2^7`-core finite cyclotomic min (HYP-3581) | the genus-1 global correction |

`dim(global side) = genus`. For genus 0 the right column is empty and the problem is the left column (solved).
For genus 1 the right column is one-dimensional — one cusp form — and that is the whole difficulty. This is
the synthesis: **the project has been studying the local–global gap of `X_0(2p)`, and its dimension is the
genus.** Everything computable (metagraph, signed cycle index, CV(H), the transitive collapse, the finite
cyclotomic floor) is the LOCAL column; the one missing thing is the single GLOBAL mode.

## Other synthesizing points (procedurally generated)

- **Genus 1 ⟺ the doublet binds; genus ≥ 2 ⟺ larger cores bind.** Computed: at `p=7` the binding core is
  the doublet (gap `0.198 = 4cos^2(3pi/7)`); at `p=11,13` (genus 2) cores like `{0,1,2,3,7}` bind *below*
  the doublet (gaps `0.0078, 0.0049`). So genus 1 is special: the obstruction is the **simplest possible**
  configuration (a 2-element doublet, THM-578). Beyond genus 1 the obstruction is irreducibly complex. This
  is the sharp reason **`N=14` is the last bespoke-tractable case** — not just "Mersenne∩Heegner∩3mod4," but
  "the last `N` whose obstruction is a doublet."
- **`nu_2 = 0 ⟺ Paley`, and it keeps the genus clean.** For apex `≡3 mod 4` (3,7,11) there are no order-2
  elliptic points, so `genus = 1 + psi/12 − (#cusps)/2` with no fractional corrections — the genus is a
  clean integer driven by `psi(2p)/12`. The Paley/Borsuk-Ulam pillar (THM-581) is what makes the genus
  well-behaved.
- **The floor constant is plausibly a period / special L-value of `f_14`** (BSD-flavored, speculative): if
  the floor `= ` an integral against `f_14`, its value relates to `L(f_14,1)` or the real period `Omega` of
  `14a`. Worth checking whether `inf R' = 0.344` or `4cos^2(3pi/7)` carries an `Omega_{14a}` factor. (Run
  the persistence test before believing any numeric match — the polysemy discipline.)
- **anti-Littlewood = the global obstruction to local vanishing.** Littlewood (`inf q‖qa‖‖qb‖ = 0`) is
  local-at-every-scale; the LRC floor's positivity (HYP-3551) is the GLOBAL cusp-form obstruction to that
  vanishing. Anti-Littlewood IS "the genus-1 mode forbids the local product from reaching 0."
- **`14a`'s bad reduction at `2,7` = the local cusp factors.** Conductor `14 = 2·7`; the reduction types at
  `p=2,7` (the local `L`-factors / Tamagawa numbers) are the local data at the doubling and apex cusps. The
  global object (the curve) is assembled from these local pieces plus the one global mode — the local-global
  picture again.

## The shift in seeing

The deepest readjustment: the project's many beautiful structures are not parallel discoveries — they are
the **single local–global dichotomy of `X_0(2p)`**, viewed in different charts, and the **genus measures the
global side**. LRC(2p) is solved exactly when genus `= 0` (everything is local/boundary); it is first hard
at genus `1` (`N=14`), where the obstruction is one cusp form `f_14`, concentrated as a vanishing-rate at the
doublet apex cusp; it is intractably complex at genus `≥ 2`. The whole proof is now one finite local
computation (the `2^7`-core cyclotomic min) plus one global number (the apex-cusp leading coefficient of
`f_14`).

See HYP (this session), HYP-3586 (cusps=Klein, hardness=genus), HYP-3581 (the finite cyclotomic floor),
[[the-cusps-are-the-klein-group-hardness-is-the-genus]], [[the-right-frame-audit-when-the-proof-becomes-finite]],
[[polysemous-constants-bridges-traps-and-homonyms]], THM-584, THM-578, HYP-3547, HYP-3550/3551.
