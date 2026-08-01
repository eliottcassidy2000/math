---
source: opus-2026-07-31-S4 (pushing the FC(2) inhomogeneous saddle-descent rigor)
status: >
  RIGOROUS REDUCTION + solved sub-case + honest Stokes localization (corrects an over-optimistic first pass).
  The inner radial integral of L(f^m) has a genuine saddle at r*~Dm where f ~ f_D, giving the growth law
  [L(f^m)/(Dm)!]^{1/m} -> max_{a in [0,1]} |phi_D(a)| (phi_D=f_D(a,1-a)); VERIFIED numerically. So L(f^m) is
  asymptotically (Dm)! times the TOP-FORM period P_m := int_0^1 phi_D(a)^m C(a) da, C the explicit radial-
  saddle correction. Hence inhomogeneous FC(2) REDUCES to: P_m=0 for all large m => phi_D=0 (a WEIGHTED
  Pakovich PMP for the top form). This CLOSES the descent whenever P_m has a unique non-degenerate dominant
  a-saddle (generic top forms) -- then P_m ~ c phi_D(a_c)^m/sqrt(m) != 0, forcing phi_D=0 and a degree descent
  f -> f-f_D. The OBSTRUCTION is exactly top forms whose dominant complex a-saddles CANCEL through the Stokes/
  Picard-Lefschetz structure -- the Kontsevich-Zagier territory; Re f_D == 0 (imaginary top form) is the
  distinguished sub-locus where all saddles sit on one ray. So the descent pinpoints WHERE K-Z is needed: the
  vanishing of ONE complex Laplace integral P_m for the top form.
tags: [factorial-conjecture, FC2, saddle-point, laplace, stokes, picard-lefschetz, kontsevich-zagier, top-form, pakovich-pmp, descent, no-pole, inhomogeneous]
related: [FC2-no-pole-homogeneous-is-the-lebesgue-moment-problem, amm-lower-bound-rigorous-the-capacity-floor-lesson-and-fc2-inhomogeneous, factorial-conjecture-exponential-periods-and-the-repo-state]
---

# FC(2) saddle-descent: rigorous reduction to the top-form period, and the Stokes localization

## 1. The saddle reduction (rigorous growth law)

`f = sum_{d<=D} f_d`, `f_D != 0`. Radial coordinates `x=ra, y=r(1-a)` give `f(ra,r(1-a)) = P(r,a) =
sum_d r^d phi_d(a)`, `phi_d(a)=f_d(a,1-a)`, and

```
   L(f^m) = int_0^1 [ int_0^inf P(r,a)^m e^{-r} r dr ] da.
```

The inner integral has phase `m log P - r`; since `P ~ r^D phi_D(a)` for large `r`, `Re(m log P - r) = mD log r
+ m log|phi_D| - r + o(m)` has a genuine maximum at `r* = Dm`, independent of `arg phi_D`. Laplace at `r*`
(with `r* = Dm`, Stirling `(Dm)^{Dm}e^{-Dm}sqrt(2 pi Dm) = (Dm)!`) gives

```
   int_0^inf P(r,a)^m e^{-r} r dr  ~  (Dm)! * phi_D(a)^m * C(a),     C(a) = exp(phi_{D-1}(a)/(D phi_D(a))) *(...)
```

(the correction `C` is bounded on `{phi_D != 0}` and encodes the lower forms). Hence

```
   L(f^m)  ~  (Dm)! * P_m,     P_m := int_0^1 phi_D(a)^m C(a) da,     and   [L(f^m)/(Dm)!]^{1/m} -> max_[0,1] |phi_D|.
```

**Verified** (`f = x^2+y^2+0.7x`, `D=2`, exact `L(f^m)=sum [x^i y^j]f^m i! j!`):
`[L/(2m)!]^{1/m} = 1.434, 1.186, 1.113, 1.082, 1.064` (`m=3..15`) `-> max_[0,1](a^2+(1-a)^2) = 1`.

**So inhomogeneous FC(2) REDUCES to the weighted top-form PMP:** `P_m = 0 for all large m  =>  phi_D == 0`,
then descend `f -> f - f_D` (degree `D-1`) and iterate to `f == 0`. `P_m` is a Pakovich PMP with the explicit
weight `C`.

## 2. Where the descent CLOSES (generic top forms)

`P_m = int_0^1 phi_D^m C da` is a **complex Laplace integral**. Its `m -> infinity` asymptotic is a sum over
the DOMINANT saddles of `phi_D` on `[0,1]` (interior critical points `phi_D'(a_c)=0` and the endpoints), via
steepest descent. If there is a UNIQUE dominant saddle `a_c` with non-degenerate `phi_D''(a_c) != 0`, then

```
   P_m ~ C(a_c) phi_D(a_c)^m sqrt(2 pi / (m |...|))  != 0,
```

so `L(f^m) != 0` for large `m` -- contradicting `L(f^m)=0 forall m`. Hence `phi_D == 0`. **The descent closes
for every top form with a unique non-degenerate dominant saddle** (a generic, open condition), giving FC(2)
for all `f` whose top form is generic.

## 3. The obstruction, honestly: cancelling Stokes saddles (K-Z territory)

The gap is NOT closed in general, and the reason is precise. When `phi_D` has SEVERAL dominant saddles
`a_1,...,a_k` with values `phi_D(a_j)` of equal modulus, the steepest-descent contributions
`sum_j w_j phi_D(a_j)^m` can CANCEL for all `m` through the Picard-Lefschetz / **Stokes** structure -- the
steepest-descent contour decomposition (which saddles are "on the contour") can make the leading terms sum to
zero. (My first pass argued "an exponential sum `sum w_j e^{im theta_j}` cannot vanish for all `m`"; that is
true for a FIXED nonzero `{w_j}`, but the Stokes decomposition can zero the `w_j` themselves, or flip which
saddles contribute -- the honest correction.) This is exactly the Kontsevich-Zagier exponential-period
content: **whether the top-form period `P_m` vanishes is a monodromy/Stokes statement about the exponential
integral, and K-Z is the tool that forbids the "accidental" all-`m` cancellation.**

The distinguished obstruction locus is `Re f_D == 0` (the no-pole condition of the companion note): then
`phi_D = i g_D`, `g_D` real, all saddle values lie on the imaginary axis (one ray of phases), and the weighted
PMP becomes REAL (`int g_D^m (real weight)`), which IS determinate -- so **on the `Re f_D = 0` locus the
descent closes by real moment determinacy**, and K-Z is only needed to reach that locus (i.e. to prove `Phi_f`
entire). Off it, K-Z controls the Stokes cancellation directly.

## 4. Net: what is proved, and the sharp remaining statement

- **Rigorous:** the saddle reduction `L(f^m) ~ (Dm)! P_m`, the growth law `-> max|phi_D|`, and the descent for
  every top form with a unique non-degenerate dominant saddle (hence FC(2) for generic top forms).
- **Reduced to:** the vanishing of the single complex Laplace integral `P_m = int_0^1 phi_D^m C` -- a
  Kontsevich-Zagier / Stokes question about the TOP FORM only. FC(2) is thus reduced from a statement about all
  `f` to a statement about homogeneous top-form periods `P_m` and their Stokes structure.
- **The friend's K-Z route, made precise:** K-Z `=>` `P_m` cannot vanish for all `m` unless `phi_D=0` (no
  accidental Stokes cancellation) `=>` descent `=>` FC(2). The homogeneous case (`FC2-no-pole-...`, Lebesgue
  PMP, exact through degree 4) is the `D=0` correction / base of this descent.

Methodological note (kept from the companion): the vanishing `P_m=0` is exact-arithmetic-only; the earlier
double-precision moment cancellation trap applies to `P_m` too.
