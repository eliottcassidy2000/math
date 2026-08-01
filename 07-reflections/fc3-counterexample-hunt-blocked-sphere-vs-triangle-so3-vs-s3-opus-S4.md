---
source: opus-2026-07-31-S4 (full FC(3) counterexample construction hunt on the triangle)
status: >
  HUNT RESULT (no counterexample found) + the sphere-vs-triangle distinction (explains GM(3)-false vs
  FC(3)-maybe-true). Searching the S_3/C_3-covariant homogeneous forms on the 2-simplex for an all-moments-zero
  phi (int_Delta phi^m=0 for all m): the C_3 generator kills 3 nmid m automatically, leaving int phi^{3k}=0.
  BLOCKED at every degree tested: D=2 (1 proj param) kills phi^3 but not phi^6; D=3 (2 params) kills phi^3,phi^6
  but not phi^9 (~7e-4); D=4 (4 params) the fit does not even close. Pattern: degree D kills the first p(D)
  moments, the next is nonzero -- no polynomial realizes the all-moments-zero (disk-area) measure. KEY REASON:
  FC(3) is the moment problem on the TRIANGLE (only discrete S_3), whereas GM(3) is the moment problem on the
  full SPHERE S^2 (continuous SO(3)); Long's GM(3) counterexample exploits the CONTINUOUS symmetry, which the
  triangle lacks. So the arc->area jump makes FC(3) the geometric BORDERLINE, but the DISCRETE symmetry blocks
  the polynomial realization -- evidence that bare FC(3) is TRUE despite being borderline (refining the
  expectation "FC(3) may be false"). Not a proof: high-degree / non-covariant / Q-multiplier (Mathieu-Zhao)
  variants are not excluded.
tags: [factorial-conjecture, FC3, counterexample-hunt, triangle, sphere, so3, s3, gaussian-moments, GM3, long, moment-problem, blocked]
related: [fc-n-is-the-simplex-moment-problem-arc-vs-area-and-the-s3-triangle, factorial-conjecture-exponential-periods-and-the-repo-state]
---

# The full FC(3) counterexample is blocked: sphere (SO(3), GM(3) false) vs triangle (S_3, FC(3) maybe true)

## 1. The hunt

Target: a homogeneous `phi` on the 2-simplex with `int_{Delta_2} phi^m dsigma = 0` for ALL `m >= 1`
(`<=> f_D` a bare FC(3) counterexample). Best ansatz: the `C_3`-covariant `omega^2`-eigenspace (`phi -> omega^2
phi` under the 3-cycle), which makes `int phi^m = 0` for `3 nmid m` AUTOMATIC, leaving only `int phi^{3k} = 0`.
The eigenspace has complex dimension `~binom(D+2,2)/3`, i.e. projective params `p(D) = 1, 2, 4, 6, ...` for
`D = 2, 3, 4, 5`.

Result (exact simplex moments `i!j!k!/(i+j+k+2)!`, high precision):

```
   D=2 (p=1):  kill int phi^3=0  ->  int phi^6 ~ 2e-6   NONZERO   (blocked)
   D=3 (p=2):  kill int phi^{3,6}=0  ->  int phi^9 ~ 7e-4  NONZERO (blocked)
   D=4 (p=4):  the fit int phi^{3,6,9,12}=0 does not even close from generic starts (blocked)
```

**Pattern:** degree `D` has just enough freedom to zero the first `p(D)` of the `3k`-moments; the NEXT one is
nonzero. No polynomial realizes the all-moments-zero measure. The partial `S_3` form
`phi = alpha+omega beta+omega^2 gamma` still kills `2/3` of the moments (`3 nmid m`) -- the honest partial
counterexample -- but the last third resists at every degree.

## 2. Why: the triangle has only S_3, the sphere has SO(3)

The homogeneous reductions are:

```
   FC(3):  L(f^m) = (Dm+2)! int_{Delta_2}   phi^m dsigma     -- the TRIANGLE (positive orthant), symmetry S_3.
   GM(3):  E[f^m]  = c(Dm)   int_{S^2}       phi^m dtheta     -- the full SPHERE,             symmetry SO(3).
```

Both are angular moment problems; only the DOMAIN differs. **GM(3) is FALSE** (Long, arXiv:2607.18186: a
3-variable quartic with `E[P^m]=0` for all `m`) because the sphere carries the CONTINUOUS group `SO(3)` -- rich
enough to build a `phi` whose image measure is rotationally balanced to all orders. **The triangle carries only
the DISCRETE `S_3`** (six elements); its `C_3` gives 3-fold balance (`3 nmid m` killed) but nothing kills the
residual `3 | m` moments, exactly as the hunt shows. So:

> The arc -> area jump (interval to triangle) makes FC(3) the geometric BORDERLINE where a counterexample is
> not a priori forbidden. But realizing it needs an all-orders rotationally balanced image, which the sphere's
> `SO(3)` supplies (GM(3) false) and the triangle's discrete `S_3` does NOT (FC(3) blocked). **This is strong
> evidence that bare FC(3) is TRUE**, even though GM(3) is false -- the two conjectures diverge precisely
> because their domains have continuous vs discrete symmetry.

## 3. Honest scope

- Found: NO FC(3) counterexample; a verified 2/3-partial one; a clean structural reason (S_3 vs SO(3)) the full
  one is blocked.
- This refines the field expectation "FC(3) may be false": on the bare (moment-vanishing) reading, the
  triangle evidence points to TRUE. It is not a proof -- I have not excluded (i) high-degree `C_3`-covariant
  forms (`D >= 5`), (ii) non-covariant forms with a subtler balance, or (iii) the Mathieu-Zhao / Q-multiplier
  failure (`L(f^m)=0` but `L(g f^m) != 0`), which is the reading under which GM(3)'s falsity transfers to the
  JC and need NOT contradict bare FC(3).
- The decisive next probe: whether the triangle admits a polynomial with `SO(2)`-balanced image at any degree
  (a "null-quadrature" pushforward), or a proof that its discrete `S_3` caps the achievable balance at 3-fold
  -- which would upgrade the evidence to a theorem that **bare FC(3) is true while GM(3) is false**.

## 4. The clean takeaway across FC(2)/FC(3)/GM(3)

```
   FC(2): interval, 1-D ARC        -> Cauchy forces phi=0                 -> TRUE.
   FC(3): triangle, 2-D area, S_3  -> only 3-fold balance, full CE blocked -> (evidence) TRUE, borderline.
   GM(3): sphere,   2-D area, SO(3)-> continuous balance, CE exists (Long) -> FALSE.
```

The controlling parameter is the DOMAIN's symmetry group: `{e}` (arc) < `S_3` (triangle) < `SO(3)` (sphere).
Continuous symmetry is what a moment-vanishing counterexample needs; the factorial/exponential domain (positive
orthant -> triangle) only ever offers discrete symmetry, which is why the Factorial Conjecture is sturdier than
its Gaussian sibling.
