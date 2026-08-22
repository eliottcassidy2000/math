# Cross-frontier synthesis: support passports, laminar cuts, and branch labels

**Status: SYNTHESIS OF PROVED / VERIFIED-EXACT / FINITE-EXACT RESULTS.**  This
note connects THM-3308--3313 and the concurrent `z1=216` audit.  It introduces
no reduction between LRC, HFC/FC, and JC and closes none of the conjectures.

**Subsequent continuation.**  THM-3321 completes the proposed cyclic-quartic
support-four exclusion, so every surviving candidate now uses all five
coefficients; THM-3323 gives the exact degree-21 support-five quotient
dimension `1670` without closing that chart.  THM-3319 releases both JC clutch
slopes `(d,k)` as an algebraic etale germ, without transporting THM-3312's
particular cofactor ratio.  THM-3320 completes the next four LRC prefix rows,
moving `353/31` wall rows/families to `349/29`.  The statements below retain
their original scoped proposal boundaries.

## The common pattern

The source snippet's gauge warning is the right abstraction: two coordinate
systems may carry isomorphic invariant data without naming the same physical
objects.  Across the three active lanes, the useful move is

```text
labelled object -> invariant passport -> native operation response,
```

followed by an explicit audit of what the passport destroys.

| lane | labelled object | invariant passport | native operation | destroyed coordinate |
|---|---|---|---|---|
| cyclic HFC quartic | five complex Fourier coefficients | support face, shortest closed coefficient-phase covering width, modulus ratios | next surviving moments `M3,M6,M9,...` | coefficient phases if support alone is kept |
| alternating HFC sextic | projective point `[A:B:C]` | even-moment ideal | add the next even moment | individual common zeros after elimination |
| JC exceptional fibre | one of two roots in `B/A` | trace, norm, difference square | quadratic conjugation / strict transform | branch label |
| LRC projected wall | labelled body and ray states | four-bit marginals and nested capacity tails | exact status screen | physical entry, endpoint, owner, phase, current |

In every row the passport is strong enough for a one-way obstruction, but
reversing it would identify gauge data with a labelled physical system.

## Three versions of the transverse-information tax

### 1. Moment support needs phase and one more scale

THM-3310 says a cyclic quartic candidate must use at least four of its five
Fourier monomials.  That support statement alone discards phases.  The cubic
coamoeba sidecar restores part of them: shortest closed coefficient-phase
covering width must exceed `pi/3`, and
the second-largest modulus must exceed `1/23` of the largest.  These two
coordinates cut different directions around every projective coordinate axis.

THM-3311 supplies a scale analogue.  The alternating sextic equations
`I2=I4=0` retain a genuine degree-eight survivor; adding `I6` makes the ideal
unit.  The lesson is not “three moments always suffice,” but that a survivor
field should be carried to the next native moment before its elimination is
scalarized.

### 2. Trace/norm needs the sheet generator

On the JC exceptional quadratic, trace and norm completely describe an
unordered conjugate pair.  They cannot choose a branch.  THM-3309 proves the
fixed-slice pointed deck and true-gradient obstruction; THM-3312 universalizes
its trace/norm algebra and independently replays the anti-descent identity.
The projective elimination-cofactor ratio `rho` satisfies

```text
(rho-sigma(rho))^2=256 P^4 delta/Nm(W)^2 != 0.             (1)
```

Thus `rho` itself is exactly the missing quadratic sidecar.  Symmetrizing it
restores descent but destroys the label needed to continue one root.  This is
the same type error as treating gauge-identified `C12` coordinates as one
physical vertex system.

There is a second type separation: the elimination pair is unimodular, but the
gradient pair vanishes.  “Cofactor” without its source module is too coarse;
only the latter could feed a Keller mate.

### 3. Marginals need two or three nested thresholds

Most LRC status states die by one coordinate-union inequality.  THM-3308 proves
that the eleven
previously weighted cores do not.  They become elementary only after retaining
two nested capacity thresholds, with one genuine three-threshold case.  The
law is

```text
sum_(t in A) 1_{c(P)>=t} <= c0+sum_i c_i 1_{i in P}
  => sum_(t in A) T_t <= c0 q+sum_i c_i m_i.               (2)
```

Equation `(2)`, applied canonically in THM-3308 and to the new rows in THM-3313,
is an operation-response passport: it remembers how one status
pattern behaves at several capacity scales.  A single marginal or a
solver-selected scalar dual hides that structure.  Even the enriched passport
still loses physical entry and current, so its conclusion remains projected.

## A reusable research rule

When a quotient obstruction nearly closes a family, do not immediately enlarge
the quotient or brute-force the whole ambient space.  First ask for the
smallest response passport under the next native operation:

1. identify the labelled source and the quotient map;
2. list the invariant coordinates already retained;
3. apply one new moment, threshold, normal direction, or conjugation;
4. record the smallest sidecar separating surviving labels; and
5. test a hostile object known to survive the coarser passport.

The rule has three counterindications.

- Do not use it when the next operation is not intrinsic to the target.
- Do not call an unordered passport a root/runner/vertex selection.
- Do not infer physical sufficiency from emptiness of a necessary projected
  relaxation.

## Underreported next tests

The synthesis changes the ranking of follow-up work.

### HFC/FC

- On each cyclic-quartic support-four hyperplane, combine the strict coamoeba
  and lopsided-axis exclusions with exact saturation before attempting a full
  three-parameter Groebner basis.  The failed direct exact and mod-67 runs show
  that raw elimination is already the wrong first representation.
- Intersect tropical initial forms of `M3,M6,M9,M12` with the five coefficient
  valuation chambers.  A missing chamber would be a proof-grade support-four
  exclusion; a surviving chamber supplies a sharply normalized hostile.
- In the next `S3` representation cell, compute the two-moment survivor field
  first and ask which later moment is load-bearing, rather than assuming the
  first pair is decisive.

### JC(2)

- Release one of `d,k` while carrying `rho`, not only `Tr(rho),Nm(rho)`, and
  test whether the quadratic cover ramifies, splits, or meets an owner wall.
- Pull the inverse-graph cofactor module to the cover and keep elimination and
  gradient modules separate.  The cheapest decisive gate is gradient
  unimodularity before any divergence-class computation.
- Use trace/norm only for descent statements; use the cross-square `(1)` as the
  anti-descent control.

### LRC(14)

- Compile all four-bit monotone capacity events into minimal one/two/three-
  threshold modular covers, with explicit counterexamples outside each
  template.
- Run the next complete ruler prefix while retaining index `64` as a nonempty
  control.
- For physical progress, return to the horn with a lawful two-axis square that
  retains endpoint origin and current.  Another scalar status refinement cannot
  manufacture those coordinates.

These tests are deliberately asymmetric: the algebraic lanes can profit from
exact elimination after the passport shrinks the chart, while LRC must restore
physical ancestry before a projected certificate can affect the conjecture.
