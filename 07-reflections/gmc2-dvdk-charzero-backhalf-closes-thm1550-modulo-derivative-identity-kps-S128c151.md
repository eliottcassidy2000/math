# The char-0 back half: THM-1550 closed exp/log-free, modulo exactly the derivative identity

*kind-pasteur-2026-07-22-S128c151. Owner: keep working collaboratively to finish the one analytic
lemma. This session I claimed a specific non-colliding sub-piece (broadcast first), formalized it
kernel-pure, and composed it onto mac-mini's S165 reduction — turning the sole survivor `h(0,t)=1`
into the strictly simpler **exp/log-free** derivative statement `d_t(h(0,t))=0`.*

## The collaborative split (broadcast, then built)

mac-mini (S165) had just reduced **both** DvdK routes to the single scalar identity `h(0,t)=1`
(`coeff_zero_smallRootFactor_mul_unit : P.coeff 0 · h(0,t) = −t·r₀`, so `Π = c·t ⟺ h(0,t)=1`), and
flagged the remaining work as the `[x⁰]`-Laurent / Fredholm-determinant content. death-star (S114) had
taken the additive TOP-companion Lagrange (`=1`) and the dihedral synthesis, explicitly **staying off
the [x⁰] crux** and holding Lean pending a split.

I broadcast a claim on the **char-0 back half of death-star's exp/log-free route** — the piece nobody
was on — then built it. This is exactly the division mac-mini welcomed: they own the Laurent front
half, I own the char-0 closing.

## Delivered (both kernel-pure `[propext, Classical.choice, Quot.sound]`, wired into root)

**`GMC2DvdKCharZeroClosing.lean`** (self-contained, imports only Mathlib — Mathlib-PR-quality):
- `coeff_eq_zero_of_derivativeFun_eq_zero` — over a char-0 field, `derivativeFun f = 0 ⟹` every
  positive-index coefficient vanishes (`coeff_derivativeFun` + `(n+1)` a nonzero factor in char 0).
- `eq_C_of_derivativeFun_eq_zero` — hence `f = C (constantCoeff f)`. **This converse is *not* in
  Mathlib** (Mathlib has only the forward `derivativeFun_C`/`derivativeFun_one`); it is a clean,
  reusable contribution in its own right.
- `eq_one_of_derivativeFun_eq_zero` — with `constantCoeff g = 1`, `derivativeFun g = 0 ⟹ g = 1`.
- `factorCoeff0_eq_of_unit_eq_one` — trivial closing: `g=1 ∧ P₀·g = rhs ⟹ P₀ = rhs`.

**`GMC2DvdKMultiplicativeClosing.lean`** (composition onto mac-mini's Weierstrass objects):
- `smallRootFactor_coeff0_eq_of_derivative_vanishes` — given `h(0,0)=1` and `d_t(h(0,t))=0`,
  `(smallRootFactor R M).coeff 0 = −t·r₀`. A three-line term-mode composition of my char-0 closing
  with mac-mini's `coeff_zero_smallRootFactor_mul_unit`. Hence `Π = (−1)ᴹ P.coeff 0 = c·t`.

## The mathematical point — the survivor got strictly simpler

mac-mini/death-star both scoped the last identity as the `[x⁰]`-in-annulus **`h(0,t)=exp(−∑D_m tᵐ/m)`**
(Fredholm determinant, `log det(1−tf) = −∑ tr(fᵐ)tᵐ/m`). death-star's insight was that the exp/log is
**avoidable**: differentiate instead of integrate. This session formalizes exactly that. Because
"zero formal derivative ⟹ constant" holds in char 0, we never form `exp` or `log`:

> the sole survivor `h(0,t)=1` ⟸ `d_t(h(0,t)) = 0` (the log-derivative identity under `D_m=0`) `∧ h(0,0)=1`.

So the remaining crux is no longer a transcendental series identity but the **derivative-form**
`d_t(h(0,t)) = 0`, i.e. `d_t(h(0,t))/h(0,t) = −∑_{m≥1} D_m t^{m-1} = 0` under `D_m=0` — the root-free
`[x⁰]`-Laurent log-derivative identity, mac-mini/death-star's frame lane, with **no exp, log, Puiseux,
or Fredholm determinant** needed to finish from it.

## Honest overlaps and the infra fix

- `factorCoeff0_eq_of_unit_eq_one` (abstract) overlaps mac-mini's concrete
  `coeff_zero_smallRootFactor_mul_unit`; theirs is the load-bearing one. My **net-new** is the char-0
  converse (`eq_C_of_derivativeFun_eq_zero` and friends) — not in Mathlib, not done by the fleet — and
  the composition that makes the exp/log-free closing explicit on the real objects.
- **Fleet-wide infra fix:** the root `TournamentH7.lean` on `origin/main` carried unresolved
  git-conflict markers (`<<<<<<< / ======= / >>>>>>>`, from death-star-S114's push), so
  `lake build TournamentH7` (full library) was un-parseable for everyone. I resolved it (kept both
  valid imports `GMC2Thm2067Reduced` + `GMC2FullRootPhi`) and wired in the three DvdK-closing modules.

## Remaining (named)
1. `d_t(h(0,t)) = 0` under `D_m=0` — the `[x⁰]`-Laurent log-derivative identity in derivative form
   (mac-mini/death-star frame lane). Exp/log-free now.
2. `h(0,0)=1` — the normalization, a direct consequence of the distinguished factor `P ≡ Xᴹ mod t`
   (so `h ≡ 1 mod t`); dischargeable from mac-mini's Weierstrass structure.

Both are hypotheses of `smallRootFactor_coeff0_eq_of_derivative_vanishes`; discharging (1) closes the
multiplicative route (and, via the shared Abel-duality content, feeds the additive `b=1` wrapper).

## Cross-links
`GMC2DvdKCharZeroClosing` / `GMC2DvdKMultiplicativeClosing` (this session) ·
`GMC2DvdKWeierstrass.coeff_zero_smallRootFactor_mul_unit` (mac-mini-S165) ·
`GMC2GeneratingFunction.generatingFunction_eq_one` (boxeph-S240: `F(t)=1`) · `GMC2TopLagrange`
(death-star-S114) · THM-1550 · THM-2101 §6 (Abel duality) · HYP-9000 (my consolidation) · HYP-9005.
