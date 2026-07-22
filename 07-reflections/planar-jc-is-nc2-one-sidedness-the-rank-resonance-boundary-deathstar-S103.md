# The binary symmetric-Hessian subcase is one-sided; the rank/resonance comparison is heuristic

**death-star-2026-07-21-S103** (HYP-8905). Owner: work to prove the planar Jacobian Conjecture, pulling in past
threads creatively. Target is genuinely **JC(2)** — THM-1300 supplies the
verified dimension-three counterexample. The proved result here is narrower:
the binary homogeneous symmetric-Hessian stratum is one-sided and tame. The
rank/resonance and GMC2 cycle comparisons are proposed analogies only.

> **Scope correction (MISTAKE-229).** There is no proved
> `NC2 => GMC(2) => JC(2)` chain. The symmetric reduction of a general planar
> Keller map lands in four variables, so the exact two-variable calculation
> below is a classified symmetric subcase, not the missing reduction target.
> The rank/cycle comparison is a heuristic and does not classify Keller
> collisions.

## 1. The 2-D nilpotent Hessian is NC2 one-sidedness (proved, verified)
De Bondt–van den Essen reduce JC to the symmetric case `F = x + ∇P`, `P` homogeneous, **Hessian nilpotent**;
Zhao's VC relates Hessian nilpotence to `Δ^m(P^m) = 0 ∀m`. This is formally
moment-like, but it is not the Gaussian functional `E[P^m]` without an
additional transform. In two
variables the `2×2` symmetric Hessian is nilpotent ⟺ `trace = ΔP = 0` (harmonic) **and** `det = 0`. For harmonic
homogeneous `P = A z^d + B z̄^d` (`z=x+iy`) the exact computation gives
```
    det(Hess P) = −4·d²(d−1)²·A·B·|z|^{2(d−2)}  =  0   ⟺   A·B = 0.
```
So **`Hess P` nilpotent ⟺ `P` is one-sided** — `P = A z^d` (holomorphic) or `B z̄^d` (antiholomorphic), a power of
the isotropic linear form `x±iy`. This is *verbatim* NC2's conclusion (a nullcone member is charge-one-sided),
realized inside the Jacobian problem. Verified numerically (d=3,4,5): two-sided `A=B=1` has `det ≈ 15.6, 6.7, 2.0
≠ 0`; one-sided has `trace=det=0`.

**These maps are tame.** For `P = A z^d`, `∇P = Ad·z^{d-1}·(1,i)`, so `F₂ − iF₁ = −iz` is *linear* — the map is
one-fiber-linear (codex THM-2063), and `z = F₁ + iF₂` is recovered linearly, then `x,y` — **invertible**. And
one-sided `P` is holomorphic hence harmonic (`Δ = 0`), so `Δ^m(P^m)=0` holds trivially: **Zhao's VC is automatic
on the one-sided locus, exactly as GMC2's moments vanish there.**

## 2. Rank threshold analogy (not a boundary theorem)
A `2×2` nilpotent matrix (`trace=det=0`) has **rank ≤ 1** — a *single* isotropic direction. A nilpotent `n×n`
can reach rank `n−1`, so **dim ≥ 3 admits rank ≥ 2 = multiple isotropic directions**. This is precisely where the
theory breaks:
- **dim 2:** nilpotent Jacobian is forced rank ≤ 1 ⟹ rows `∇H₁, ∇H₂` are **parallel** (verified for both the
  symmetric `v=(1,i)` isotropic case and the non-symmetric `v=(0,1)` case) ⟹ `H₁,H₂` functionally dependent ⟹
  one-fiber pencil (codex THM-2063) ⟹ tame. One direction = **one-sided**.
- **dim ≥ 3:** nilpotent matrices can have rank at least two. The verified
  Keller counterexample of THM-1300 also lives in dimension three, but no
  theorem here derives its collision from that rank fact.

**The analogy.** S101's unique-versus-coincident cycle split suggests comparing
rank one with multi-direction resonance. Both programs have an easy one-sided
pole and a harder multi-object region. They are not one nullcone: Zhao's
Laplacian iterates and GMC2's Gaussian expectation are different functionals,
and no transfer theorem is supplied here. The comparison may schedule searches
for a rank-sensitive certificate; it proves no Jacobian statement beyond
Section 1.

## 3. Honest scope and what it buys
- **Proved (verified):** the symmetric planar Jacobian / 2-D Zhao-VC / Hessian-nilpotent case — `Hess P` nilpotent
  ⟹ `P` one-sided ⟹ `x+∇P` one-fiber-linear tame. This is a genuine result, and it *is* NC2's conclusion.
- **NOT proved:** full JC(2). The relevant general symmetric reduction lands
  in four variables. The rank-one observation and THM-2063 close special
  planar strata; they are not a reduction of every planar Keller map to a
  remaining non-symmetric rank-one case.
- **What it buys:** one more exact empty cell in the planar atlas and a
  rank-sensitive analogy worth testing. A genuine transfer of GMC2
  noncancellation would still require an explicit map from Gaussian channels
  to Laplacian or Keller data and a proof that it preserves the relevant
  vanishing predicate.

Cross-links: THM-1830 (NC2/GMC2 conditional map, not a JC implication), THM-1435 (Zhao VC, Hessian-nilpotent), THM-2063 (codex
one-fiber-linear planar Keller — = my one-sided maps), THM-1300 (Alpöge counterexample, dim 3), S101/HYP-8878
(unique vs coincident cycle), S102/HYP-8879 (LRC=GMC2 resonance), S217/218 (arithmetic entropy), memory
`gmc2-domination-dead-fact-is-algebraic`, `nc2-gmc2-lean-formalization-state`. Script
`04-computation/planar_jc_is_nc2_one_sidedness_deathstar_S103.py` (+ `.out`). HYP-8905.
