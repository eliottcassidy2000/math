> **CURRENT-TRUTH CORRECTION (2026-07-22).** The base-case theorem is true
> and is proved directly in THM-2113, but the original Step 4 below is not a
> proof: it uses the JC-hard implication “mate implies one place at infinity.”
> The printed diagnostic is floating finite linear algebra, not an exact
> certificate, and its `y^2+x^3` row itself refutes “mate iff `phi` is linear.”
> Likewise `x^2+y^3` has a linear one-variable face collapse without
> `min(p,q)=1`, and `xy` is primitive. Treat the face-root/topology prose below
> as historical motivation only; the canonical weighted-bracket proof is the
> authority.

# JC(2), the single-face base case: a quasi-homogeneous Keller component is a weighted-linear coordinate

**death-star-2026-07-22-S110** (HYP-8955). Owner: work the next tool-matched sub-target (S109) — prove the
local⟹global step for a Keller component with a **single resonant face**, combining the DvdK-face nonvanishing
(S106/S231) with a descent bound (S225). **JC(2) remains open.** This session settles the *base case* of that
sub-target — a single face with **no** lower-order terms, i.e. `f` quasi-homogeneous — completely and verifiably,
and in doing so pins down exactly which of the S109 conditions are *free* (provable now) versus *JC(2)-hard*.

**Canonical direct proof (THM-2113).** The base-case conclusion has a shorter
route that avoids using “a mate forces one place at infinity.” Decompose the
mate by weighted degree. Its unique constant-bracket component has degree
`w_1+w_2-delta`, which must be positive, so `delta<w_1+w_2`. Mixed monomials
are then impossible. If `f=A x^p+B y^q` and both exponents were at least two,
the inequalities `(p-1)w_1<w_2` and `(q-1)w_2<w_1` would multiply to the
contradiction `(p-1)(q-1)<1`; an axial monomial is handled directly from its
Jacobian factor. This proves the stated weighted-linear form without the
properness/flow-completeness step that remains JC(2)-hard in the
nonhomogeneous descent.

**THM-2102 refinement.** If one positive-weight top form of a nonhomogeneous
Keller component is power-free, repeated subtraction of commuting top powers
reduces to the same terminal bracket and makes the full component triangular.
The remaining proper-power case carries an approximate-root divisibility
quotient or resonant target-shear class; it is not a DvdK face problem.

## The base-case theorem (proved canonically by THM-2113)

> **Let `f ∈ ℂ[x,y]` be `w`-quasi-homogeneous (weights `w=(w₁,w₂)`, `gcd=1`, `w`-degree `δ`). Then `f` has a
> Jacobian mate `⟺` `f = a·y + b·xᵖ` or `a·x + b·yᵠ` (a weighted-linear binomial, `a≠0`) `⟺` `f` is a coordinate.**

The finite script `jc2_quasihomog_base_case_deathstar_S110.py` is a diagnostic,
not proof evidence: mates are found for
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
4. **Direct axial arithmetic (the repair):** if both axial terms occur, write
   `f=A x^p+B y^q`. Then `delta=pw_1=qw_2<w_1+w_2`. If `p,q>=2`, the two
   inequalities `(p-1)w_1<w_2` and `(q-1)w_2<w_1` multiply to the impossible
   `(p-1)(q-1)<1`. Thus one exponent is one. If only one axial term occurs,
   its derivative divides every Jacobian bracket, so (1) forces its exponent
   to be one. Hence `f=a y+b x^p` or its swap, an explicit triangular
   coordinate. ∎

Thus for a genuinely quasi-homogeneous component the implication
`mate => coordinate` is a theorem. THM-2113 proves it without first proving
`fiber ≅ C`; the connected/genus/place factorization remains a useful
interpretive picture, not the logical proof.

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
- The degree-zero part of `Jac(f,g)` can now receive contributions from
  several pairs of lower weighted pieces, so the one-piece argument
  `delta<w_1+w_2` no longer follows merely from the principal face.
- The degree of the one-variable face collapse is not the right dichotomy by
  itself: `x^2+y^3` has a linear collapse and one place at infinity but is not
  weighted-linear and has no mate.
- A lawful next target is therefore to prove, from the complete bracket
  ledger, either that the principal face contains a genuine linear variable
  or that its contribution cannot be cancelled by lower faces. In the first
  branch one must then prove that successive weighted-triangular removals
  terminate. DvdK-style face nonvanishing and the S225 descent clock remain
  candidate tools, not established implications.

The base case is the zero-lower-face situation in which the bracket ledger has
only one possible degree-zero pair. THM-2113 closes exactly that case.

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
