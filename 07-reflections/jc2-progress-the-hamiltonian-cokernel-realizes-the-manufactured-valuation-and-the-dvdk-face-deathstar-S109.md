# JC(2) progress: the Hamiltonian cokernel realizes the S107 manufactured-valuation and the S106 DvdK-face — and grounds them in fiber theory

**death-star-2026-07-22-S109** (HYP-8950). Owner: work the S107/S108 ideas into progress toward a planar
Jacobian-conjecture proof. **JC(2) is open; this is not a proof.** It is concrete progress: an operator
reformulation (verified) that (i) turns S107's speculative "manufactured-valuation orbit-product" into a *verified*
weight/valuation obstruction, (ii) makes S106's DvdK transfer concrete (the resonant faces are literal DvdK
constant-term conditions), (iii) pins the genuine residual to a single global invariant, and — the honest payoff —
(iv) shows the whole S106/S107 route is the *same obstruction* as the classical fiber-`≅ℂ` approach
(Abhyankar–Moh–Suzuki, Kaliman), de-speculating S107 by grounding it in established theory.

## 0. The reformulation (known, verified in-repo)

For `f ∈ ℂ[x,y]`, `f` is a **Keller component** (has a Jacobian mate `g`, `jac(f,g)=1`) `⟺` `{f,g}=1` is solvable
`⟺` **`1 ∈ im(D_f)`**, where `D_f := {f,·} = f_x∂_y − f_y∂_x` is the Hamiltonian derivation. The obstruction module
is `coker(D_f) = ℂ[x,y]/{f,ℂ[x,y]}` — the Brieskorn/relative-de-Rham module, whose generic rank is
`dim H¹(generic fiber)`. `1 ∈ im(D_f) ⟺ [dx∧dy]=0` there. Standard facts assemble to:

> **JC(2) `⟺` every Keller component is a coordinate `⟺` [ `f` has a mate `⟹` generic fiber of `f` is `≅ ℂ` ].**

(A mate forces `df≠0` everywhere, so fibers are smooth; the conjecture is that they are `ℂ`, not a higher-genus or
punctured curve. "Fiber `≅ℂ ⟹` coordinate" is Kaliman 2002 / Abhyankar–Moh–Suzuki.) **Verified** (`mate_exists`
by exact linear algebra, script): coordinates `x, x²+y, y+x²` have mates; `x²` (fiber = 2 lines), `xy` and `x+x²y`
(fiber `≅ℂ*`), `x²+y³` (genus > 0), `x³` (3 lines), and homogeneous deg-2 forms have **none** — exactly matching
"generic fiber `≅ℂ`."

## 1. S107's manufactured valuation, realized and verified: the weight obstruction

S107 proposed transferring S106's *manufactured valuation* to JC(2) but left the crux open. Here it is, concrete:
each positive weight `w=(w₁,w₂)` is a valuation at infinity (`v_w(x^iy^j)=iw₁+jw₂`). From the bracket-degree count,
`{f,g}` has `v_w ≥ v_w(f)+v_w(g)−w₁−w₂`, so `{f,g}=1` (which has `v_w=0`) forces `v_w(f) ≤ w₁+w₂`. Contrapositive:

> **Valuation obstruction.** If some positive weight `w` has `v_w(f) > w₁+w₂` (every monomial of `f` strictly above
> the line `iw₁+jw₂=w₁+w₂`), then `f` has **no mate** — so it is not a Keller component.

**Verified**: catches `x²` (`w≈(2,1)`), `x²+y³` (`w≈(3,2)`), `x³` — all correctly no-mate; and correctly finds *no*
obstructing weight for the coordinates `x, x²+y`. This is the exact JC-analogue of S106/THM-2067's valuation
contradiction (`v(Π)=1 ≠ 0`): there the manufactured `t`-adic valuation opposed a norm; here a weight/valuation at
infinity opposes the `jac=1` normalization. It generalizes/reframes codex's THM-2045 grading obstruction to *every*
positive weight.

## 2. S106's DvdK object, realized: the resonant faces

The weight obstruction is *sufficient, not necessary*: `xy` and homogeneous forms sit **on** a resonant face
(`v_w(f)=w₁+w₂` exactly), dodging §1. There the question "`1 ∈ im(D_{f^{w\text{-face}}})`" for the `w`-quasi-
homogeneous face is a **1-variable constant-term condition** — the `ℂ^*`-weight action collapses a `w`-quasi-
homogeneous polynomial to a Laurent polynomial in one variable `τ`, and the face-Hamiltonian's surjectivity onto
the constant is exactly a **DvdK-type `CT ≠ 0` nonvanishing** (the S106/S101 object). Concretely for a homogeneous
face of degree `n≥2`: `{f_n,g}` is homogeneous of degree `n+deg g−2`, so hitting the constant needs `deg g = 2−n<0`
— impossible; the face always obstructs at leading order, and any mate must be supplied by lower-order terms.
**This makes S106's transfer precise**: the counterexample-excluding condition on each resonant face is a DvdK
constant-term nonvanishing, to which the S106 orbit-product / monomial-certificate (boxeph S231) apply, and the
sweep over weighted faces is exactly boxeph S225's **descent-termination** (Newton-slope / coprime-interval Euclid).

## 3. The genuine residual, pinned

`x+x²y = x(1+xy)` has generic fiber `≅ ℂ^*` (parametrize by `x∈ℂ^*`), so **no mate** — yet it **dodges every
weight** (Newton polygon `{(1,0),(2,1)}` is low) *and* is not a single resonant face. So neither §1 nor a single §2
face sees it: the obstruction is the **global** invariant `coker(D_f) =` fiber cohomology (`H¹(ℂ^*)=ℂ≠0`), which no
one valuation detects. This is the precise JC(2) core: a Keller component must have `1∈im(D_f)` *globally*, i.e.
generic fiber `≅ℂ` — and the valuation/face obstructions (§1–2) are the *local* shadows of that global condition.
The descent (§2) terminates `⟺` the local conditions force the global one — the open content.

## 4. What this buys (honest scope)

- **Not a proof, and largely a synthesis of known objects** (Hamiltonian reformulation; Brieskorn module;
  fiber-`≅ℂ ⟹` coordinate = Kaliman/AMS; the quasi-homogeneous degree obstruction = THM-2045-flavored). All §0–3
  facts are verified in-repo.
- **The genuinely useful step:** S107's manufactured-valuation route and S106's DvdK transfer were *speculative*;
  here they become a **verified weight obstruction** (§1) and a **literal DvdK-face condition** (§2), and both are
  shown to be *local shadows of one global invariant* `coker(D_f) = H¹(\text{fiber})` (§3) — i.e. the S106/S107
  route and the classical **fiber-`≅ℂ`** route are the *same obstruction*. That grounds the repo's algebraic route
  in established JC theory and says precisely where it can and cannot bite: valuations/faces catch the "high" and
  "resonant" strata (the resonance dictionary's tame + face-resonant poles), but the `ℂ^*`-fiber residual is
  global and needs the descent to terminate — the real frontier, matching boxeph S225.
- **Concrete next target:** prove the local→global step for the first nontrivial class — Keller components whose
  Newton polygon forces a *single* resonant face — by combining the §2 DvdK-face nonvanishing (S106 orbit-product /
  S231 certificate) with a descent-termination bound (boxeph S225's coprime-interval Frobenius number). This is a
  well-posed, tool-matched sub-problem, not the whole conjecture.

Cross-links: S107 (resonance dictionary + manufactured-valuation route — §1 realizes its crux, §3 locates its
residual), S106 (orbit-product/valuation — §2 is its literal transfer), S108 (unit distances = same cancellation
side; the DvdK-face here is the cancellation object), boxeph S225 (descent-termination = the §2 face sweep),
boxeph S231 (monomial certificate — the effective face nonvanishing), codex THM-2045 (grading obstruction = §1
special case), THM-2063 (one-fiber-linear tame class = the fiber-`≅ℂ` pole), codex DC2/Poisson thread (the
Hamiltonian `{f,·}` framing), external: Kaliman 2002 / Abhyankar–Moh–Suzuki (fiber-`≅ℂ ⟹` coordinate), Brieskorn
module. Script `04-computation/jc2_hamiltonian_derivation_obstruction_deathstar_S109.py` (+ `.out`). HYP-8950.
