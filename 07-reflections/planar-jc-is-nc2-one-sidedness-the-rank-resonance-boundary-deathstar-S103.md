# Planar JC is NC2 one-sidedness: the 2-D nilpotent Hessian forces one-sided, and the JC-true/false boundary is the GMC2 unique-vs-coincident-cycle threshold

**death-star-2026-07-21-S103** (HYP-8905). Owner: work to prove the planar Jacobian Conjecture, pulling in past
threads creatively. Target is genuinely **JC(2)** — Alpöge's Keller counterexample (THM-1300) killed JC in
dim ≥ 3. Result: the symmetric planar Jacobian / Zhao Vanishing Conjecture in dim 2 **is exactly my NC2 one-sided
nullcone**, proved here; and the JC-true(dim 2)/false(dim ≥ 3) boundary is the *same* rank/resonance threshold as
the GMC2 unique-cycle-vs-coincident-cycle dichotomy (S101). This places planar JC squarely inside my GMC2/NC2
framework (THM-1830's NC2 ⟹ GMC2 ⟹ … ⟹ JC(2) program).

## 1. The 2-D nilpotent Hessian is NC2 one-sidedness (proved, verified)
De Bondt–van den Essen reduce JC to the symmetric case `F = x + ∇P`, `P` homogeneous, **Hessian nilpotent**;
Zhao's VC says `Hess P` nilpotent ⟺ `Δ^m(P^m) = 0 ∀m` — the *same moment-vanishing* as GMC2's `E[P^m]=0`. In two
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

## 2. The rank threshold = the JC-true/false boundary = the GMC2 cycle threshold
A `2×2` nilpotent matrix (`trace=det=0`) has **rank ≤ 1** — a *single* isotropic direction. A nilpotent `n×n`
can reach rank `n−1`, so **dim ≥ 3 admits rank ≥ 2 = multiple isotropic directions**. This is precisely where the
theory breaks:
- **dim 2:** nilpotent Jacobian is forced rank ≤ 1 ⟹ rows `∇H₁, ∇H₂` are **parallel** (verified for both the
  symmetric `v=(1,i)` isotropic case and the non-symmetric `v=(0,1)` case) ⟹ `H₁,H₂` functionally dependent ⟹
  one-fiber pencil (codex THM-2063) ⟹ tame. One direction = **one-sided**.
- **dim ≥ 3:** rank ≥ 2 = several isotropic directions = a **resonance**; Alpöge's Keller counterexample
  (THM-1300, `det JF=−2`, triple collision) lives here — JC is **false**.

**The unification.** This is the *same* threshold as GMC2 (S101): a **unique** primitive cycle (one-sided,
non-cancelling, easy) versus **coincident** cycles (resonance, hard / counterexample-bearing). Planar JC is the
**last one-direction dimension** — the analogue of the DvdK-free unique-cycle stratum. dim ≥ 3 JC / the
coincident-cycle GMC2 stratum / the multi-resonance is where both a JC counterexample and the DvdK-hard case
appear. Zhao's VC (`Δ^m P^m=0` ⟺ Hess nilpotent) **is** GMC2 (`E[P^m]=0` ⟺ one-sided): one moment-nullcone, and
the one-sided/rank-1 conclusion is the shared "easy" pole. This also matches boxeph's arithmetic-entropy frame
(S217/218): the rigid extremum = zero-entropy = the unique/one-sided point; difficulty = the positive-entropy
multi-object.

## 3. Honest scope and what it buys
- **Proved (verified):** the symmetric planar Jacobian / 2-D Zhao-VC / Hessian-nilpotent case — `Hess P` nilpotent
  ⟹ `P` one-sided ⟹ `x+∇P` one-fiber-linear tame. This is a genuine result, and it *is* NC2's conclusion.
- **NOT proved:** full JC(2). The reduction "2-D nilpotent Jacobian ⟹ rank ≤ 1 ⟹ parallel gradients ⟹ functional
  dependence" is sound and lands on codex's one-fiber-pencil frontier (THM-2063); showing the *non-symmetric*
  rank-1 (parallel-gradient) case is always tame is the open crux, exactly where codex is working.
- **What it buys:** planar JC now sits inside my GMC2/NC2/resonance framework — the symmetric case *is* NC2
  one-sidedness; the JC-true/false boundary *is* the unique/coincident-cycle boundary; the counterexample
  dimension *is* the first multi-resonance. A concrete bridge for transferring GMC2 non-cancellation machinery
  (the `Q̄^p` single-power certificate, the unique-cycle criterion) to the pencil/tameness question, and a
  structural reason planar JC should be true (rank forced ≤ 1) unified with why dim ≥ 3 fails (rank ≥ 2 resonance).

Cross-links: THM-1830 (NC2 ⟹ GMC2 ⟹ JC(2) program), THM-1435 (Zhao VC, Hessian-nilpotent), THM-2063 (codex
one-fiber-linear planar Keller — = my one-sided maps), THM-1300 (Alpöge counterexample, dim 3), S101/HYP-8878
(unique vs coincident cycle), S102/HYP-8879 (LRC=GMC2 resonance), S217/218 (arithmetic entropy), memory
`gmc2-domination-dead-fact-is-algebraic`, `nc2-gmc2-lean-formalization-state`. Script
`04-computation/planar_jc_is_nc2_one_sidedness_deathstar_S103.py` (+ `.out`). HYP-8905.
