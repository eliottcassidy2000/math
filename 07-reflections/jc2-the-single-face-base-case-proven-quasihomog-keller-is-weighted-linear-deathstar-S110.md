# JC(2), the single-face base case proven: a quasi-homogeneous Keller component is a weighted-linear coordinate

**death-star-2026-07-22-S110** (HYP-8955). Owner: work the next tool-matched sub-target (S109) — prove the
local⟹global step for a Keller component with a **single resonant face**, combining the DvdK-face nonvanishing
(S106/S231) with a descent bound (S225). **JC(2) remains open.** This session settles the *base case* of that
sub-target — a single face with **no** lower-order terms, i.e. `f` quasi-homogeneous — completely and verifiably,
and in doing so pins down exactly which of the S109 conditions are *free* (provable now) versus *JC(2)-hard*.

## The base-case theorem (proven, verified)

> **Let `f ∈ ℂ[x,y]` be `w`-quasi-homogeneous (weights `w=(w₁,w₂)`, `gcd=1`, `w`-degree `δ`). Then `f` has a
> Jacobian mate `⟺` `f = a·y + b·xᵖ` or `a·x + b·yᵠ` (a weighted-linear binomial, `a≠0`) `⟺` `f` is a coordinate.**

Verified (`jc2_quasihomog_base_case_deathstar_S110.py`, `mate_exists` by exact linear algebra): mates exist for
`y+x²`, `y+x³`, `x+y²`; and fail for every non-weighted-linear quasi-homogeneous form — `x²+y²`, `x²+xy+y²`,
`(x+y)²`, `xy`, `x²y+xy²`, `y²+x³`.

**Proof (⟹), as a ladder that realizes the three repo tools:**
1. **`f` is primitive** — the FREE lemma (no JC): if `f=P(h)` with `deg P ≥ 2` then `{f,g}=P'(h){h,g}`, so `=1`
   forces `P'(h)` a unit, hence `P` linear — contradiction. (Verified: `x²`, `(x+y)²`, `x²+2x` all no-mate.)
   So a mate forces **connected** generic fibers.
2. **`δ < w₁+w₂`** — the S107 **weight/valuation** obstruction, exactly. `{f,g}=1` (which has `v_w=0`) forces
   `v_w(f) ≤ w₁+w₂`; for `w`-homogeneous `f`, `v_w(f)=δ`, so `δ ≤ w₁+w₂`, and `δ=w₁+w₂` needs `g` of `w`-degree 0
   (a constant, giving `{f,·}=0`), so **`δ < w₁+w₂`**. Geometrically this is *genus 0* (it excludes `y²+x³`, the
   cusp, `δ=6 > 5`, whose fibers are genus 1).
3. **`δ < w₁+w₂` forces all monomials onto the axes** (`i,j ≥ 1 ⟹ iw₁+jw₂ ≥ w₁+w₂ > δ`), and `w`-homogeneity of
   a single degree leaves at most one `xᵖ` and one `yᵠ`: `f = a xᵖ + b yᵠ` (a binomial). Its **face polynomial**
   `φ` — the 1-variable collapse `f = xᵅyᵝ·φ(τ)`, `τ = x^{w₂}/y^{w₁}` — is the S106/S101 **DvdK face object**, and
   for a binomial `φ` has degree `= p/w₂ = q/w₁`.
4. **The DvdK-face is single-branch:** primitive + a mate forces `min(p,q)=1` (`φ` linear ⟺ one branch at
   infinity ⟺ one place ⟺ `≅ℂ`); the alternatives are exactly the *fiber ≇ ℂ* failures — `φ` with `≥2` distinct
   roots is `≥2` places at infinity (`x²+y²`), a repeated root or a monomial factor is non-primitive
   (`(x+y)²`, `xy`), matching (1). Hence `f = a y + b xᵖ` = a coordinate. ∎

So for **one face**, the S109 slogan is a *theorem*: `mate ⟺ fiber ≅ ℂ ⟺ coordinate`, with the three tool-conditions
(primitive / `δ<w₁+w₂` / `φ` linear) each supplying one of connected / genus-0 / one-place-at-infinity. The
`ax^p+b y^q` genus/place count is the classical binomial-curve fact; here it is the DvdK single-branch condition.

## The honest sharpening of S109

Working the base case forced a correction to how S109 read. Of the S109 conditions on a Keller component:
- **FREE (provable now, no JC):** `f` primitive ⟹ connected generic fibers (§1). This is genuine but weak —
  `xy` is primitive with fiber `ℂ*` and no mate, so connectedness alone misses the punctures.
- **JC(2)-HARD (not free):** "`mate ⟹ generic fiber ≅ ℂ`" and its corollary "leading form is a power of one
  linear form." The natural proof (flow the Hamiltonian `X_f`, `X_f(g)=1`, to trivialize the fiber to `ℂ`)
  **has a real gap — completeness of the flow**, equivalently *properness of `g` on the fibers*, which is exactly
  the JC(2) crux. So S109's "`mate ⟹ fiber ≅ ℂ`" is the *reformulation's open direction* (Kaliman supplies only
  the converse `fiber ≅ ℂ ⟹ coordinate`), not a free lemma.

**The base case is precisely where that hard direction becomes provable** — because for quasi-homogeneous `f` the
fiber topology is *explicit* (a binomial/`μ_δ`-cover), so properness is not needed; the weight bound `δ<w₁+w₂`
(§2) replaces flow-completeness.

## The descent (the remaining sub-target) and where the tools go

For a general **single principal face** `f = Φ_w + (lower w-degree)`, with `Φ_w` the `w`-leading face:
- If `Φ_w`'s DvdK-face `φ` is **not** linear (`≥2` distinct roots), the principal branch structure already gives
  `≥2` places at infinity — the multi-root/*resonant* pole of the dictionary; the counterexample-excluding claim
  here is a DvdK-face nonvanishing (S106 orbit-product / boxeph S231 certificate), applied per face.
- If `φ` **is** linear, `Φ_w` is weighted-linear and the classical **Abhyankar–Moh** descent (a weighted-
  triangular automorphism killing the top face, then recursing) reduces `f` toward the base case; the number of
  steps is boxeph S225's **coprime-interval / Lamé** descent-termination bound. Termination ⟹ coordinate.

So the sub-target reduces to two tool-matched pieces, both now *named and located*: (a) **multi-root face ⟹ no
mate** (a per-face DvdK nonvanishing — the honest content, since "`≥2` places ⟹ no mate" is itself the JC-hard
direction and must be earned face-by-face, not assumed), and (b) **linear-face descent terminates** (AM + S225).
The base case is (a)+(b) with zero descent steps.

## Honest scope

- **Proven + verified:** the base case (quasi-homogeneous Keller `⟹` coordinate). This is a small/known-flavored
  stratum (weighted-linear coordinates are classical), but it is the first case where S109's `mate ⟹ fiber ≅ ℂ`
  is *earned* rather than assumed, and it realizes all three tools (primitivity, S107 weight, S106 DvdK-face) as
  the three factors of `fiber ≅ ℂ`.
- **Corrected:** S109's "`mate ⟹ fiber ≅ ℂ` / single-root leading form" is the JC(2)-hard direction (flow-
  completeness / properness gap), not free; the free lemma is primitivity only.
- **Open:** the descent for a nontrivial principal face (≥1 lower-order step) — (a) the per-face DvdK nonvanishing
  beyond binomials, and (b) an unconditional AM/S225 termination bound. JC(2) remains open.

Cross-links: S109 (the Hamiltonian reformulation this sharpens; §"honest ladder" corrects its free-vs-hard split),
S107 (the weight obstruction = §2; the resonance dictionary's single-root pole = §4), S106 (the DvdK face `φ` =
§3–4), boxeph S225 (descent-termination = the linear-face recursion), boxeph S231 (monomial certificate = the
per-face nonvanishing), THM-2063 (one-fiber-linear tame class), external: Abhyankar–Moh–Suzuki (the descent /
one-place-at-infinity), Kaliman (fiber `≅ℂ ⟹` coordinate). Script
`04-computation/jc2_quasihomog_base_case_deathstar_S110.py` (+ `.out`). HYP-8955.
