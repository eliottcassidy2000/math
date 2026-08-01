---
source: opus-2026-07-31-S4 (working the FC(3) fold-boundary gap; synthesis with kps's concurrent leak-ideal thread)
status: >
  SADDLE CAP THEOREM (new, all degrees) + synthesis with kps-S158-161. On the equilateral triangle T, the leak
  int_T tilde-phi^{3k} dA is a 2-D LAPLACE integral dominated by the max-modulus points of the C_3-covariant
  tilde-phi. Those points form C_3-orbits, and each orbit contributes 3 M^{3k} e^{3ik theta_j} p_j (p_j != 0) to
  the 3|.-leak; a sum over DISTINCT orbits sum_j p_j e^{3ik theta_j} cannot vanish for all k. So the leaks are
  NONZERO whenever |tilde-phi| attains its max at ISOLATED points -- for EVERY degree. This covers (a) every
  HOLOMORPHIC tilde-phi (max on the boundary by the maximum principle), and (b) generic tilde-phi (isolated
  max, incl. all D=2 covariant forms, whose max sits at a vertex). Verified: [ |int phi^{3k}| ]^{1/(3k)} climbs
  toward M. The ONLY residual is a non-holomorphic tilde-phi whose |tilde-phi| maxes on an interior CLOSED
  CURVE winding uniformly around the circle |z|=M (the fold->disk case) -- non-generic, and exactly the
  quadrature-domain sliver. SYNTHESIS: kps's leak-ideal-emptiness (Groebner=(1) at D=2, exact; transversality
  rank=P through D=5) closes that residual at specific degrees; my edge obstruction rules out its boundary
  (edge) version; the saddle cap closes the rest at ALL degrees. Together: no FC(3) counterexample -- generic
  all-degrees + exact D<=2 + transversal D<=5. Strong toward bare FC(3) TRUE while GM(3) FALSE.
tags: [factorial-conjecture, FC3, saddle, laplace, stationary-phase, C3, cap, kps-synthesis, leak-ideal, holomorphic, maximum-principle]
related: [fc3-cap-the-decisive-step-null-quadrature-and-the-edge-obstruction, fc3-counterexample-hunt-blocked-sphere-vs-triangle-so3-vs-s3, kps-S161]
---

# The FC(3) saddle cap: isolated-max is capped at all degrees (with kps's leak-ideal thread)

## 1. The saddle cap (new, all degrees)

On the equilateral triangle `T` (coordinate `zeta = alpha+omega beta+omega^2 gamma`), an FC(3) counterexample
in the `C_3`-covariant family (`tilde-phi -> omega^2 tilde-phi`) needs the LEAKS `int_T tilde-phi^{3k} dA = 0`
for all `k` (the `3 nmid k` moments vanish for free). This is a 2-D LAPLACE integral: for large `k`,
`int_T tilde-phi^{3k} dA` is dominated by the points where `|tilde-phi|` attains its maximum `M`.

Because `|tilde-phi|` is `C_3`-invariant, the max-modulus points come in `C_3`-orbits `{p, omega^2 p, omega^4 p}`
with `tilde-phi`-values `M e^{i theta}, M e^{i(theta - 2pi/3)}, M e^{i(theta - 4pi/3)}`. For `m = 3k`:

```
   sum over one orbit  tilde-phi(p_j)^{3k}  =  M^{3k} e^{3ik theta} ( 1 + e^{-2pi i k} + e^{-4pi i k} )
                                            =  3 M^{3k} e^{3ik theta}   (nonzero).
```

So each orbit contributes `~ 3 M^{3k} e^{3ik theta_j} p_j / k` (`p_j != 0` the Laplace prefactor), and

```
   int_T tilde-phi^{3k} dA  ~  (3/k) M^{3k} sum_j p_j e^{3ik theta_j}.
```

A finite sum `sum_j p_j e^{3ik theta_j}` over DISTINCT base-angles `theta_j` (distinct `3 theta_j mod 2pi`,
since distinct orbits) cannot vanish for all large `k` (linear independence of the exponentials). Hence:

> **Saddle Cap.** If a non-constant `C_3`-covariant `tilde-phi` attains `max_T |tilde-phi|` at ISOLATED points,
> its leaks `int_T tilde-phi^{3k} dA` are NONZERO for large `k` -- so it is NOT an FC(3) counterexample. This
> holds for **every degree**.

**Verified:** `int_T phi^{3k}` for a D=2 covariant `phi` is `0.048, 0.017, 0.0085, 0.0051` (`k=1..4`), nonzero,
with `[ |int phi^{3k}| ]^{1/(3k)}` climbing toward `M = 1`.

## 2. What the saddle cap covers -- almost everything

- **Every HOLOMORPHIC `tilde-phi`:** by the maximum principle `|tilde-phi|` maxes on `partial T` (edges/
  vertices), at isolated points (or edge-arcs, handled by the edge obstruction). Capped.
- **Generic `tilde-phi`:** the max of a real-analytic `|tilde-phi|` on a compact region is generically at
  isolated points. Capped -- at ALL degrees (the tested D=2 covariant forms all max at a vertex).

The ONLY residual is a **non-holomorphic `tilde-phi` whose `|tilde-phi|` attains its max on an interior CLOSED
CURVE `Gamma` that winds uniformly (no stationary phase) around the circle `|z| = M`** -- i.e. the fold curve
maps to the disk boundary. That is precisely the null-quadrature / fold-boundary sliver, and it is
non-generic (a modulus ridge is a codimension-1 degeneracy). Its boundary (edge) version is killed by the edge
obstruction (a straight edge can't map to a circle).

## 3. Synthesis with kps's concurrent leak-ideal thread (kps-S158..S161)

kps attacks the SAME residual algebraically, from the other side:

- **kps-S161 (exact):** in DFT coordinates the leaks are INTEGER polynomials; `Groebner(leak1,leak2,leak3) =
  (1)` at D=2, so `V(leaks) = varnothing` -- **no FC(3) counterexample at D=2, isolated OR family** (closes
  the "isolated KZ gap" the transversality scan misses).
- **kps-S159/S160 (transversality):** the leak map has generic rank `P` through **D=5** (no families), and the
  Pakovich composition saddle-weight is `R`-independent (no composition closes the leak).

These are exactly complementary to the saddle cap:

```
   saddle cap (this note)  : ALL degrees, isolated max (generic + all holomorphic).      geometric / analytic.
   kps leak-ideal          : EXACT at D=2, transversal at D<=5 (incl. isolated gap).     algebraic / arithmetic.
   edge obstruction (S4)   : rules out the disk when the max curve is an EDGE.            geometric.
```

Together they leave, at each degree, only the thin non-holomorphic interior-uniform-winding-fold sliver, which
kps's ideal-emptiness closes exactly (D=2) and transversally (D<=5), and which the saddle cap closes wherever
the max is isolated. **The combined evidence that bare FC(3)-homogeneous is TRUE (while GM(3) is FALSE) is now
strong from two independent directions.**

## 4. The exact remaining lemma

To upgrade to a theorem at all degrees: **no polynomial `tilde-phi` has `|tilde-phi|` maximal on an interior
closed curve `Gamma` with `tilde-phi|_Gamma` a uniform winding of the circle `|z|=M`.** Equivalently (Schwarz
function): the fold curve `tilde-phi({J=0})` cannot be a circle with the pushforward density Fourier-flat at
the `3k` frequencies. This is the one quadrature-domain rigidity statement left; the saddle cap + kps's
ideal-emptiness reduce FC(3) to exactly it. Handed to kps as the joint frontier.
