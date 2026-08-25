# Exceptional affine/simple-zero packet: what survives after the parity repair

**Date:** 2026-08-24  
**Status:** research reflection; THM-4054 is the truth source for promoted claims

## Portfolio and inheritance

- **Anchor:** understand the first boundary not covered by THM-4046's
  critical-displacement obstruction, namely `H'(0)!=0` on the same four
  exceptional folds.
- **Niche:** compare arbitrary two-forms, polynomial exact two-forms, and the
  nonlinear Darboux cone after passing to one complete retained multigerm.
- **Wildcard:** deform the explicit affine local pair through
  `H_epsilon=t+epsilon t^2` and ask whether its failure is transverse or only
  a bad fixed-coordinate gauge.

The closest proved mechanism was THM-4046's retained cokernel, the canonical
hostile was its nonzero order-eight response for `H in t^2 C[t]`, and the
least-used sidecar was the global/local decomposition information discarded
by a linear two-form packet. The corrected near miss was an attempted import
of THM-3629's even-fold mixing theorem.

## Concept board

| object | predicate | operation | lost coordinate | cheapest test |
|---|---|---|---|---|
| exceptional fold `Q_alpha` | simple-zero boundary | evaluate retained derivatives | parity | test `Q_alpha'(0)` |
| local pair `(a,-4c)` | constant Jacobian | perturb `H` | globality of `a` | compute `Jac_(x,q)(a,c)` |
| arbitrary two-form packet | constant density | pull back and truncate | closedness/decomposition | exact rank and augmentation |
| polynomial exact-form packet | same retained image? | exterior derivative | global Darboux factorization | compare contained ranks |
| fixed-`a` and mixed pair jets | nonlinear continuation | linearize in `epsilon` | higher deformation obstructions | cutoff-five response, then `epsilon^2` |

Every pull changed the board. The derivative evaluation removed the proposed
global mixing inference. Exact-form saturation showed that closedness is not
the missing finite coordinate. The fixed-`a` response isolated a genuine
obstruction, while the full tangent bank showed that the obstruction belongs
to the gauge `F=a`, not to the whole mixed-pair deformation problem.

## The repaired anatomy

The exceptional polynomial satisfies

```text
Q_alpha(-1,0,1)=(-3,-3/4,-3),
Q_alpha'(-1,0,1)=(-9/2,1,9/2).
```

Thus it retains the triple collision but is not even. The shifted sign
collision used by THM-3629 is unavailable. What survives without parity is
the local algebra

```text
a=e/(b+4)=q/D^2,        Jac_(x,q)(a,c)=-3.
```

Hence `(a,-4c)` has source Jacobian `12` for `H=t`. The two-form
`da wedge dc` extends globally and exactly; the first coordinate `a` does
not. This is the first separation to keep visible:

```text
global exact form  !=  global polynomial Darboux pair.
```

At total retained order five, the three branches provide `63` coordinates.
For each of `H=t` and `H=t+t^2`, the complete arbitrary bank has `168`
columns and the complete polynomial-exact bank has `231`; both ranks are
`59`, and the constant packet lies in both images. Containment plus equal
rank proves exact-form saturation separately for each displacement.

The two rank-`59` images are nevertheless different: their relation spaces
have combined rank `6`. Equal dimensions did not identify the subspaces.
This was the session's second important repair.

For `H=t+t^2`, the fixed first coordinate `a` survives cutoff four but fails
at cutoff five. The fixed bank then has rank `57` in `63` coordinates, and
its canonical nonzero response has nonzero quartic-field norm. However, over
dual numbers the full pair

```text
F=a+epsilon f,        G=-4c+epsilon g
```

can absorb the required `-24t` correction. Its `166`-column tangent bank has
rank `59` and equals the affine arbitrary-form image. This proves a mixed
first-order escape, not nonlinear integration.

## Connection contract

The source is the arbitrary/exact/pair-tangent jet universe at the common
target germ. The target is the `63`-coordinate retained source packet. The
map is pullback followed by Taylor truncation. It preserves equality of
densities through cutoff five and exactness of the selected form. It destroys
higher jets, coherent choices across cutoffs, global polynomiality,
decomposability by two global functions, and behavior at infinity.

The missing sidecar is therefore not another scalar terminal jet. It is the
characteristic foliation/conormal class of a decomposable two-form together
with global denominator control. The cheapest decisive continuation is:

1. compute the second-order `epsilon^2` deformation obstruction;
2. repeat the tangent calculation with corrections restricted to the global
   ring `K[b,c,e,w]`, rather than local `K[a,c,w]`;
3. if that survives, test compatibility at cutoff six before increasing the
   target-degree bank; and
4. keep the separate surface-only `2 x 5` and `3 x 4` weight cells from
   THM-3592 as the first honest global Darboux targets.

## Method cards used

- **Type every analogy and implication:** the triple values matched, but the
  parity hypothesis did not.
- **Compute the repair quotient before testing the residual defect:** exact
  and arbitrary retained images had to be compared before interpreting the
  fixed-coordinate failure.
- **The same representation is not the same carrier:** equal rank `59` did
  not identify the affine and nonlinear images.
- **Use redundant paths as detectors:** exact quartic-field elimination and
  an independent all-four-root modular implementation agree in normal and
  optimized modes.

No new method card is proposed: these observations reinforce existing cards
rather than supplying evidence from a second distinct thread.

## Honest frontier

THM-4054 is a finite three-multigerm theorem plus one exact local pair and one
dual-number tangent result. No global pair is constructed. There is no
all-order lift, algebraization, Keller map, or consequence for `JC(2)`. The
strongest live signal is that the fixed-coordinate obstruction is real but
not deformation-invariant; the next obstruction, if any, must see nonlinear
decomposability or globality rather than only the retained arbitrary-form
image.
