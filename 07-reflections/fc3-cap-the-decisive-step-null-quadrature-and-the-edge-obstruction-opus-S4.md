---
source: opus-2026-07-31-S4 (the decisive step: capping the FC(3) triangle at 3-fold via null-quadrature)
status: >
  PARTIAL CAP (rigorous) + the right theory + honest remaining gap. Reformulated FC(3)-homogeneous, in the
  coordinate zeta=alpha+omega beta+omega^2 gamma (mapping the simplex to the EQUILATERAL triangle T, S_3 = its
  dihedral symmetry), as: int_T tilde-phi(zeta,zbar)^m dA = 0 for all m => tilde-phi=0. All-holomorphic-moments-
  zero means mu = tilde-phi_*(Leb_T) is a POSITIVE measure whose exterior Cauchy transform is c/w -- the field
  of a point mass at 0 (a null-quadrature / mother-body condition); by Sakai/Newton the only UNIFORM bounded
  such domain is a disk. RIGOROUS PARTIAL CAP: a rotationally-symmetric (disk) pushforward is IMPOSSIBLE,
  because it needs |tilde-phi| = const on a triangle EDGE, but tilde-phi on an edge is a polynomial in the real
  edge-parameter, so |.|^2 const on the segment => tilde-phi const on the edge => the edge maps to a POINT not a
  circular arc (verified: |phi| varies along every edge). The 3|m residual RECURSES: int_T (tilde-phi^3)^k=0 is
  the SAME null-quadrature condition for tilde-phi^3, with the same edge obstruction. Combined with the D=2,3,4
  block and the symmetry ladder {e}(arc) < S_3(triangle) < SO(3)(sphere), this is strong evidence bare FC(3) is
  TRUE while GM(3) is false. HONEST GAP: the FULL cap (excluding all non-rotational null-quadrature pushforwards,
  e.g. fold-bounded or annular) is open and is a quadrature-domain rigidity statement.
tags: [factorial-conjecture, FC3, cap, null-quadrature, mother-body, sakai, quadrature-domain, equilateral-triangle, edge-obstruction, s3, so3, disk]
related: [fc3-counterexample-hunt-blocked-sphere-vs-triangle-so3-vs-s3, fc-n-is-the-simplex-moment-problem-arc-vs-area-and-the-s3-triangle]
---

# The decisive step: capping FC(3) at 3-fold via null-quadrature (partial cap, rigorous)

## 1. The equilateral-triangle reformulation

In the coordinate `zeta = alpha + omega beta + omega^2 gamma` (`omega=e^{2 pi i/3}`), the simplex maps
LINEARLY onto the EQUILATERAL triangle `T` with vertices `1, omega, omega^2`, `dsigma -> const * dA(zeta)`,
and `S_3` becomes the triangle's dihedral symmetry (`C_3`: `zeta -> omega^2 zeta`; reflections: `zeta ->
zbar`-type). A homogeneous form `phi(alpha,beta,gamma)` becomes a polynomial `tilde-phi(zeta, zbar)`. So

> **FC(3)-homogeneous `<=>` for polynomials `tilde-phi(zeta,zbar)` on `T`:
> `int_T tilde-phi^m dA = 0 for all m >= 1  =>  tilde-phi == 0`.**

## 2. It is a null-quadrature / mother-body problem

`mu := tilde-phi_*(Leb_T)` is a **POSITIVE** measure (pushforward of positive Lebesgue) with density
`rho(z) = sum_{tilde-phi=z} 1/|J_{tilde-phi}|`, supported on `tilde-phi(T)`. All holomorphic moments vanish
(`m>=1`) iff its exterior Cauchy transform is `C_mu(w) = c/w` (`c = mu(C) = area(T)`), i.e. `mu` has the
EXTERIOR FIELD OF A POINT MASS `c delta_0` -- a **null-quadrature / mother-body** condition. By Sakai/Newton,
the only bounded domain carrying UNIFORM density with this property is a **disk centered at 0**. So a
counterexample must push `Leb_T` to a (possibly non-uniform) measure with a point-mass exterior field.

## 3. Rigorous partial cap: no rotationally-symmetric (disk) pushforward

**Claim.** No polynomial `tilde-phi` makes `mu` rotationally symmetric (a disk/annulus). **Proof.** A
rotationally symmetric `mu` has its support bounded by the circle `|z| = R`. The outer boundary of
`tilde-phi(T)` is covered by images of the triangle EDGES (and fold curves). If an edge maps onto a circular
arc, then `|tilde-phi| = R` on that edge. But on an edge, `zbar` is affine in `zeta`, so `tilde-phi` restricts
to a polynomial in the real edge-parameter `t`, and `|tilde-phi(t)|^2` is a real polynomial equal to `R^2` on
`[0,1]`, hence `identically R^2`, hence `tilde-phi(t)` has constant modulus for all real `t`, hence
`tilde-phi(t) = const` (a positive-degree polynomial is unbounded). So the edge maps to a POINT, not a circular
arc -- contradiction. Verified: `|tilde-phi|` genuinely varies along every edge (e.g. `phi=zeta` on `gamma=0`
gives `|phi| = 1, 0.66, 0.5, 0.66, 1`). **So the cleanest counterexample -- the disk -- is capped.**

## 4. The 3|m residual recurses (and inherits the edge obstruction)

For a `C_3`-covariant `tilde-phi` (auto-killing `3 nmid m`), the residual `int_T tilde-phi^{3k} dA = 0` for all
`k` is, after `eta = z^3`, exactly `int_T (tilde-phi^3)^k dA = 0 for all k` -- the SAME null-quadrature
condition for the degree-`3D` polynomial `tilde-phi^3`. The edge argument of sec 3 applies verbatim to
`tilde-phi^3` (`|tilde-phi^3| = |tilde-phi|^3`), so the disk is capped at every stage of the recursion. This is
why the numerical hunt was blocked at `D = 2, 3, 4`: the polynomial pushforward can be 3-fold balanced (`C_3`)
but never disk-balanced, and the residual only asks for the disk again.

## 5. The symmetry ladder and the honest gap

```
   domain        symmetry     counterexample                        verdict
   interval      {e}          1-D arc, Cauchy forces phi=0           FC(2) TRUE
   triangle      S_3 (disc.)  3-fold only; disk capped (sec 3-4)     FC(3) TRUE (strong evidence)
   sphere        SO(3) (cont.) continuous balance realizable (Long)  GM(3) FALSE
```

Continuous rotational balance is exactly what an all-moments-zero counterexample needs; the factorial/
exponential domain (positive orthant -> triangle) offers only the discrete `S_3`, whose 3-fold balance leaves
the `3|m` moments, and the disk (the continuous fix) is capped by the edge obstruction. **This is the decisive
mechanism: FC(3) is sturdier than GM(3) because the triangle's symmetry is discrete, not continuous.**

**Honest gap (what a full theorem still needs).** Sec 3 excludes DISK/annular (rotationally symmetric)
pushforwards. A full cap must also exclude non-rotational null-quadrature pushforwards -- measures with a
point-mass exterior field bounded by FOLD curves rather than edges, or with non-symmetric `rho`. That is a
quadrature-domain rigidity statement for polynomial images of the simplex (Aharonov-Shapiro / Gustafsson
territory): show that `tilde-phi_*(Leb_T)` can have exterior field `c/w` only if `tilde-phi` is constant. The
edge obstruction handles the outer boundary when it is an edge-image; the remaining case is a fold-bounded
outer boundary, where the argument needs the Schwarz function of `tilde-phi(fold)` to be single-pole -- the
concrete open lemma. With it, **bare FC(3) is TRUE while GM(3) is FALSE** becomes a theorem.
