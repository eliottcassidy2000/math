# The fixed-`a` failure is a triangle period, not a rank accident

**Date:** 2026-08-24
**Status:** research reflection; THM-4058 is the truth source for promoted claims

## Portfolio and inheritance

- **Anchor:** explain whether THM-4054's cutoff-five failure for `H=t+t^2`
  is isolated or the first member of a simple-zero obstruction family.
- **Niche:** replace retained-matrix elimination by the intrinsic geometry of
  the three affine branch surfaces.
- **Wildcard:** conjugate every simple-zero displacement back to the affine
  source and ask what nonlinear information survives only in the right-hand
  side.

The closest proved mechanism was THM-4054's exact response `rho`; the
canonical hostile was its surviving full mixed tangent; the corrected near
miss remained the illegal import of THM-3629's even-fold mixing conclusion.
The least-used sidecar was the ordinary-triple normalization viewpoint of
THM-3635/3696: a three-branch restriction algebra carries one compatibility
law beyond branchwise primitives.

## Concept board after the pull

| object | representation | invariant | operation | lost coordinate |
|---|---|---|---|---|
| fixed first output `a` | derivation on three graph surfaces | normalized triangle period | integrate around the branch cycle | freedom to move `F` |
| exceptional response `rho` | fifth-order vertex splitting | signed triangle area | pass from vertices to period | higher vertex jets |
| stable power `w^r` | common value plus a tiny edge deviation | `-(r/2)rho` | binomial expansion | nonlinear mixed-pair terms |
| `H=t+gamma*t^m` | affine source plus inverse derivative `h'` | `gamma binom(m,2)rho` | filtered source reparametrization | global polynomiality |

The matrix pattern first suggested a recurrence: `w^r` survived through
cutoff `r+3`, failed at `r+4`, and its response was `r` times the response on
`w`. The graph representation explains all three facts at once.

## The normalization coordinate that exposed the triangle

Writing `s=xD` turns the local chart into

```text
D=1+a s^2,       x=s/D,       q=aD^2,       c=3s+a s^3.
```

At the common target, the branches are the three roots `s=2,0,-2`.  On the
affine source, their images are graphs

```text
w=R_i(A,c)=A+lambda_i c+...,
lambda=(3/4,-1/3,-3/4).
```

The fixed-`a` Jacobian operator is exactly `-3` times differentiation in
`c` along each graph.  A common target primitive therefore exists precisely
when the three branch primitives close around their intersection cycle.
This turns a dense retained relation into the intrinsic period

```text
Pi(f)=integral_20^01 f_0 dc+integral_01^12 f_1 dc
                              +integral_12^20 f_2 dc.
```

The three pairwise vertices agree through `A^4`.  Their first edge vector is

```text
delta(5/13,-18/13,1),          delta=(26/15)rho.
```

Its `w` increments are obtained by multiplying by the three tangent slopes.
The scaled polygon has signed area `-15/52`, so division by the universal
opening factor `delta A^5` gives the response `-rho/2` on `w`.

## Why the period is complete

An annihilator alone would not prove the successful cutoffs.  The tangent
cone closes the other direction.  In target degree `d=n+1`, restriction to
the union of the three distinct planes has dimension `3d`. Tangential
differentiation has only the common `A^d` line as kernel. Its rank is
therefore `3n+2` in a `3n+3` dimensional branch-triple layer: exactly one
cokernel coordinate per degree.

The leading period row is

```text
(5/13,-18/13,1),
```

which kills both the constant and slope rows. Thus the normalized period
supplies that unique coordinate in every graded layer. This converts
vanishing of the period coefficients from a necessary test into an iff.

## The all-power and all-m laws

Along the tiny triangle, `w=w_*+O(A^5)` with `w_*=A+O(A^2)`. Expanding
`w^r` about `w_*`, the constant term integrates to zero, the linear term is
`r w_*^(r-1)Pi(w)`, and quadratic edge deviations arrive later. Hence

```text
Lambda(w^r)=-(r/2)rho A^(r+4)+higher.
```

For `H=t+gamma*t^m`, put `u=H(t)` and let `t=h(u)`. The target substitution
`w=h(u)` is filtered and invertible, while a unit source density becomes the
affine right-hand side

```text
h'(u)=1-m gamma u^(m-1)+O(u^(2m-2)).
```

The constant has zero period, the first nonlinear term hits the `r=m-1`
ladder, and every higher inverse term arrives after the decisive cutoff.
This gives the exact first response

```text
gamma binom(m,2)rho
```

at cutoff `m+3`, after success through `m+2`.

## Redundant checks and boundaries

The exact companion reconstructs the triangle over `Q(alpha)` without a
retained matrix. The independent companion instead rebuilds the complete
fixed-`a` matrices over all four split roots modulo `137`; for every
`m=2,...,12`, it sees the predicted pass/fail boundary and response
polynomial. Agreement matters because the two paths discard different
information.

The new theorem is stronger than THM-4054 only inside the fixed-output
gauge. It says nothing against a pair with both outputs moving. In fact the
quadratic mixed tangent survives at the very order where the fixed gauge
fails. The next honest question is therefore not another fixed-`a` cutoff:
it is whether the mixed first-order escape integrates to second order, and
whether the answer changes when corrections are restricted from the local
ring `K[a,c,w]` to the global target ring.

No new meta-pattern card is promoted. The triangle-period mechanism has one
strong instance; a second mathematically distinct thread is still needed
before treating it as a general research move.
