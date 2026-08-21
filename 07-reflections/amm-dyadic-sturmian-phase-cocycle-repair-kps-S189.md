# AMM dyadic Sturmian phase cocycle: repairing the scalar golden target

Status: **EXACT REFRAME + REFUTED FINITE-SCALE LAW + OPEN PHASED LIMIT**
(2026-08-21).

## 1. The hostile datum changes the state space

THM-3601 pre-registered the last failing Rule-A trace at `R=4096` as the
hostile test of the tentative scalar law

```text
q/R -> phi^(-2)/4=(3-sqrt(5))/8.                         (1)
```

THM-3616 computes

```text
(R,D0,j,q)=(4096,88,1014,382),
q/R=191/2048=0.09326171875,                              (2)
```

so

```text
q/R-(3-sqrt(5))/8=(-577+256sqrt(5))/2048
                  =-0.002229784062526... .               (3)
```

At `R=2048` the error is `(-573+256sqrt(5))/2048`, so the absolute
error increase is exactly

```text
(577-256sqrt(5)-573+256sqrt(5))/2048=1/512.             (4)
```

This refutes monotone finite-scale sharpening and the simplest scale-only
interpretation.  It does not refute an eventual limit.

## 2. The coordinate that the scalar limit discarded

Put

```text
alpha=log_5(phi^2)=2 log(phi)/log(5).                    (5)
```

The degree word at epoch `R` is

```text
d_i=floor(alpha(R+i))+D0,
delta_i=d_i-d_(i-1) in {0,1}.                            (6)
```

Consequently the entire binary step word is the mechanical word with
intercept

```text
rho_R={alpha R}.                                         (7)
```

At dyadic epochs this sidecar evolves by the exact doubling cocycle

```text
rho_(2R)={2 rho_R}.                                      (8)
```

This is not optional representation data: the compressed adjoint recurrence
chooses a different transfer operator at every `delta_i`.  The five audited
intercepts are

```text
R       rho_R
256     0.084783530353...
512     0.169567060705...
1024    0.339134121411...
2048    0.678268242821...
4096    0.356536485643... .                              (9)
```

Thus consecutive dyadic epochs are not repeated samples of one stationary
boundary condition; they are an orbit of the doubling map.

The number `alpha` is irrational.  If `alpha=m/n` were rational, then
`phi^(2n)=5^m` would be rational.  But

```text
phi^k=F_k phi+F_(k-1),
```

and `F_(2n)>0`, so `phi^(2n)` has a nonzero irrational part.  This proves the
claim.  No normality or density assertion about the particular doubling orbit
is needed here.

This also sharpens, rather than weakens, the golden connection.  Audited
THM-3027 proves that the archimedean capacity saddle has

```text
gamma*=alpha=log(phi)/log(sqrt(5)),
tau*=phi^(-2),                 (1-tau*)^2=tau*.         (10)
```

So the mechanical-word driver `alpha` and the proposed horizon fraction
`tau*/4` are two exact coordinates of the same saddle.  The new point is that
flooring the `alpha`-profile retains the intercept orbit `(8)`; the saddle
value alone does not determine a finite dyadic word.

## 3. Correct continuum object

The top-distance adjoint uses the two exact operators

```text
T_0 M=(1+z)M,
T_1 M=Pi_(>=0)[z^(-1)(1+z)^2M],                        (11)
```

driven by `delta_i(rho_R)`.  Its charges additionally depend on the Rule-A
front and clamp.  A continuum object that forgets `(7)` is therefore
undertyped.  The repaired target is a skew-product profile

```text
B(theta; epsilon, rho, kappa),                          (12)
```

where at minimum

- `theta=s/R` is normalized backward time;
- `epsilon=D0/R` (or an equivalent front-slack coordinate) records the offset;
- `rho={alpha R}` records the Sturmian intercept; and
- `kappa` records clamp saturation or the high-cell front state.

The horizon should be sought as a phase-dependent zero

```text
Theta(epsilon,rho,kappa)=first theta with B(theta;...)=0, (13)
```

not assumed in advance to be one universal constant.  A constant such as
`phi^(-2)/4` could still emerge after proving phase independence, averaging
over phase, or restricting to a returning subsequence.  None is presently
proved.

## 4. Cheap decisive experiments

1. **Phase-matched epochs.**  Compare horizons at scales whose `rho_R` are
   close instead of merely doubling `R`.  In `256<=R<=8192`, the closest
   phase matches found without any AMM trace computation are

   ```text
   rho_2048=0.678268242821...,
   rho_6023=0.678325012946...,       distance 5.6770e-5;

   rho_4096=0.356536485643...,
   rho_8071=0.356593255767...,       distance 5.6770e-5. (14)
   ```

   Terminal-offset scans at these pairs would separate scale drift from
   intercept drift.

2. **Same-scale phase sweep.**  At a fixed moderate `R`, cyclically shift the
   mechanical word or vary an abstract intercept `rho` while retaining the
   same clamp law.  Compute the adjoint zero as a function of `rho`.  This is
   cheaper than another factor-two epoch and directly tests `(13)`.

3. **Return-time subsequences.**  Search the doubling orbit `(8)` for close
   returns and compare the corresponding normalized horizons.  This tests a
   phase-conditioned limit without assuming the orbit is equidistributed.

4. **Front sidecar ablation.**  Hold the degree word fixed and vary only the
   clamp/front initialization.  If the horizon moves, `rho` is necessary but
   not sufficient and `kappa` must remain in `(12)`.

## 5. Cross-frontier mechanism

The same repair now appears independently in three active lanes.

- **LRC:** THM-3615 restores root difference and discovers a rank-two
  augmentation; marginalization had collapsed one line and created two false
  dimension anomalies.
- **AMM:** the scale-only horizon discarded the Sturmian intercept; the next
  hostile dyadic epoch then looked anomalous.
- **Jacobian cylinder:** audited THM-3611 retains the localized quotient and
  both arms before algebraic descent, while audited THM-3612 retains branch
  jets at a common collision.  Dividing by the central coordinate or keeping
  only one branch would erase the contradiction.

The reusable rule is not “never quotient.”  It is:

```text
identify the quotient kernel, carry the smallest sidecar on which the target
predicate is not constant, and only then ask for a scalar or marginal law.  (15)
```

For AMM the next honest object is the phase/front cocycle `(12)`, not another
unqualified fit of `q/R`.
