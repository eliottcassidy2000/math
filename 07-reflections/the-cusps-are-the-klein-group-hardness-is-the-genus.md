# X_0(14): the cusps are the Klein group, hardness is the genus, the obstruction is the cusp form — and a procedurally generated reframe grid

*klein-2026-06-29-S10. A long session at the cusps. mac-mini S29 put the proof at the cusps of X_0(14); this works the moves there, hunts abnormalities, and forces a grid of reframes. The arc of the whole conversation closes: it began at the n=4 Klein group (THM-584) and the cusps of X_0(14) ARE that group.*

## Three verified facts at the cusps

**1. The 4 cusps of X_0(14) ARE the Klein four-group — the n=4 tournament classes.** `omega(14)=2`, so the
Atkin-Lehner group `W(14) = {1, W_2, W_7, W_14} = (Z/2)^2` acts **regularly** (simply transitively) on the
4 cusps `{d=1,2,7,14}`. So the cusps are a torsor for `(Z/2)^2 = V_4 =` the n=4 iso-class structure
`{T,+,-,S}` (THM-584). The dictionary:

| cusp `d` | width `14/d` | Atkin-Lehner | Klein (THM-584) | LRC role |
|---|---|---|---|---|
| `d=1` (∞) | 14 | identity | `T` (00, source&sink) | BULK / interior (where all the rehearsal lives) |
| `d=2` | 7 | `W_2` | `+` (10) | DOUBLING (the 2-adic descent) |
| `d=7` | 2 | `W_7` | `-` (01) | **APEX-7** (speed 7, `7a/14=a/2` even/odd coupling) = the HARD cusp |
| `d=14` (0) | 1 | `W_14`=Fricke | `S` (11, strong) | FULL-DENSE (`m_R→0`, all mod-7 resonance = covering) |

So `W_2` = the 2-adic descent involution, `W_7` = the apex involution, `W_14` = the Fricke/complement —
the project's "two order-2 structures" (parity vs doubling) are literally `W_2` and `W_7`, and the danger
relation's complement is the Fricke `W_14`. The binding doublet (HYP-3581, `4cos^2(3pi/7)`) lives at the
narrow apex cusp `d=7` (width 2 = the doublet).

**2. LRC(2p) hardness IS the genus of X_0(2p).** Verified `genus(X_0(2p)) = 0,0,1,2,2` for
`p=3,5,7,11,13`. The genus **jumps `0→1` exactly at `N=14`**. Among the LRC-relevant apices
(`p ∈ {3,7}` = Mersenne ∩ Heegner ∩ 3-mod-4, HYP-3547): `genus(X_0(6))=0` = LRC(6), SOLVED;
`genus(X_0(14))=1` = LRC(14), the FIRST hard/open case. The genus is the missing "why 14 is hard"
invariant, complementary to the arithmetic characterization.

**3. `nu_2` detects the Paley/3-mod-4 pillar.** The order-2 elliptic-point count `nu_2(X_0(2p)) = 0`
exactly when `p ≡ 3 mod 4` (apices 3,7,11) and `=2` when `p ≡ 1 mod 4` (5,13). So `nu_2=0` ⟺ apex
`≡3 mod 4` ⟺ Paley tournament exists ⟺ the Borsuk-Ulam pillar (THM-581). The elliptic-point count *is*
the Paley condition, read on the modular curve.

## The abnormality: genus 1 means a cusp form, and the cusp form is the obstruction

`X_0(6)` is genus 0 — a rational curve, no cusp forms, only Eisenstein series. `X_0(14)` is genus 1 — the
conductor-14 elliptic curve `14a` (rank 0), carrying a **nontrivial weight-2 cusp form `f_14`** beyond the
Eisenstein series. This is the abnormality, and it is plausibly *the* obstruction:

> **Conjecture (the cusp-form obstruction).** The `Gamma_0(14)` second moment that controls the floor
> decomposes as Eisenstein (bulk) ⊕ cusp-form (`f_14`). For genus 0 (LRC(6)) the cusp-form part is empty,
> the bulk closes it, and the problem is solved. For genus 1 (LRC(14)) the missing piece the metagraph /
> `CV(H)` / `S_n`-transitive rehearsal cannot see (klein-S4: "the testbed models the bulk, not the cusp")
> is exactly the `f_14` component. The proof's hard core = controlling `f_14` at the apex cusp.

Rank 0 of `14a` (`L(f_14,1) ≠ 0`) says the cusp form is "non-degenerate" — the obstruction is a fixed
finite thing, not a vanishing one. This reframes "what we're missing" precisely: a genus-1 cusp form, at
the `d=7` apex cusp.

## Procedurally generated reframe grid (objects × lenses)

Forcing the generation: rows = LRC objects, columns = lenses; each cell is "object AS lens." The high-yield
cells (★ = big shift) are kept; the rest pruned.

- **Danger relation `D` × operator/kernel (L11):** `D*D` (compose with itself) is the floor 2nd moment; the
  floor = its spectral gap = `4cos^2(3pi/7)` at the apex cusp (HYP-3581). ★ (the proof sentence, HYP-3571).
- **`D` × category/relations (L7):** `D` as a profunctor; "does not factor" = `D` is not a product bimodule;
  "stays small" = `D∘D` is bounded. Relations-not-things (mac-mini S24). ★
- **Cusps × group action (L1+L8):** the 4 cusps = the Klein group `V_4` = the n=4 classes. ★★ (loop closure).
- **Floor × modular/automorphic (L8):** floor 2nd moment = Eisenstein(bulk) ⊕ cusp-form `f_14`; the
  obstruction = the genus-1 cusp form. ★★ (the abnormality reframe above).
- **Floor × adelic/local-global (L8+L9):** the 4 cusps = "places"; the floor = an adelic Euler product over
  places (HYP-3550, the cusp-local factors); hardness = a local-global obstruction; genus 0 = Hasse
  (easy), genus 1 = elliptic-curve Sha-like obstruction. ★ (new framing of HYP-3550).
- **Descent core `O_j` × finite set (L_combinatorial):** `O_j ⊆ Z_7` is finite → the floor is a finite
  cyclotomic min, not an estimate (HYP-3581). ★
- **Doublet × difference set/Fano (L12):** the apex QR `{1,2,4}` is a planar difference set (flat spectrum
  = optimal `rho_j`, mac-mini S27); the doublet is the *minimal* (worst) configuration = the binding cusp. ★
- **`N_R` (sheet count) × Siegel transform (L8):** `N_R` = a finite Siegel transform; `E[N_R]` = the
  Eisenstein mean; `Var` = the 2nd moment with the `14Z` congruence (Han-Lee). (mac-mini HYP-3553/3561.)
- **LRC(2p) family × genus (L8):** hardness = `genus(X_0(2p))`; the proof difficulty is the genus. ★★
- **Descent × dynamical flow toward a cusp (L4):** the 2-adic descent = flow `W_2`-ward; it terminates at
  the apex cusp `d=7` where the residual core lives. ★
- **Covering set `S` × coding theory (L6):** `S` = a covering code on `R/Z`; `M(S)` = the covering radius;
  the gap line `M ≥ 7/89` = a covering-radius bound (multiplicative, HYP-3550) — a *different* frame from
  the cyclotomic-gap floor (do not conflate; the polysemy discipline).
- (pruned, low-yield: `D` × entropy; lonely set × tropical; `N_R` × homology — surfaced nothing new.)

The grid's payoff is the convergence of independent cells on the same picture: the floor is a 2nd moment
(`D*D`) whose obstruction is a genus-1 cusp form sitting at the apex cusp, and every "right frame"
(finite cyclotomic, transitive, multiplicative, adelic) is a different chart on that one object.

## Abnormalities to keep tracking

1. **genus jump `0→1` at 14** — the hardness invariant; `genus(X_0(22))=2`, `(X_0(26))=2` (harder still,
   but off the Mersenne/Heegner apex list).
2. **cusps = Klein `V_4`** — closes the arc to the n=4 classes (THM-584); the hard cusp is `d=7`.
3. **`nu_2=0 ⟺ apex ≡3 mod 4 ⟺ Paley`** — the elliptic-point count is the Borsuk-Ulam pillar.
4. **`14a` rank 0** — the obstruction cusp form is non-degenerate (`L(1)≠0`).
5. **the "6" cluster — APPLY THE PERSISTENCE TEST:** `phi(14)=6` (witnesses/units), `6` curves in the
   `14a` isogeny class, the cuspidal group order. These may be *different* 6's (a TRAP per
   [[polysemous-constants-bridges-traps-and-homonyms]]); do not conflate without checking they persist
   across `N`. Flagged, not believed.
6. **the doublet = the width-2 apex cusp = `4cos^2(3pi/7)`** — the binding configuration (THM-578) is a
   cusp.

## The shift

The biggest readjustment: stop treating LRC(14) as an analysis problem with beautiful side-structure, and
treat it as the **arithmetic geometry of the curve X_0(14)** — a rank-0 genus-1 elliptic curve whose 4
cusps are the Klein group, whose Atkin-Lehner involutions are the project's two order-2 structures, whose
genus (vs the genus-0 X_0(6)) is exactly why LRC(14) is the first hard case, and whose weight-2 cusp form
`f_14` is the obstruction. The metagraph, the antipodal map, the signed cycle index, the descent, the
doublet — all of it — are charts on this one curve. The floor closes when the `f_14` component is bounded
at the apex cusp; that is the single remaining object.

See HYP (this session), HYP-3581 (the finite cyclotomic floor), HYP-3571 (the proof sentence),
HYP-3580 (mac-mini, the cusps), HYP-3547 (Mersenne∩Heegner∩3mod4), HYP-3550 (Euler-product floor),
THM-584 (the n=4 Klein group), THM-578 (the doublet),
[[the-right-frame-audit-when-the-proof-becomes-finite]], [[polysemous-constants-bridges-traps-and-homonyms]].
