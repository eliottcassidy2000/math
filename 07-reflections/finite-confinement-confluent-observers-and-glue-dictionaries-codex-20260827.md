# Finite confinement, confluent observers, and glue dictionaries

**Session:** codex-padic-zeta-density-cartier-20260826 continuation  
**Portfolio:** Anchor -- p-adic-zeta density source audit; Niche -- planar-JC
`W=0` attachments; Wildcard -- local/nonlocal pivot architecture.  
**Status:** synthesis of THM-4255, THM-4258, THM-4259, THM-4260, and
THM-4264, plus the finite-exact THM-4265 scout.  The gate-interior `W=0`
class sets are now empty.  The external p-adic density claims, boundary
walls, off-fibre exact `M=12`, `JC(2)`, and `DC(2)` remain open.

## The shared issue is not specialization; it is what was proved before it

This session exposed three genuinely different injectivity statements that
should not share one informal phrase.

1. **Ambient coefficient injectivity is false.**  For
   `ev_g:A[u][[f]]->A[[f]]`, THM-4255 proves
   `ker(ev_g)=(u-g)`.  One chosen section cannot recover unrestricted
   coefficients or coefficientwise order.
2. **Finite-candidate reduction can be injective.**  THM-4259 first confines
   a hidden Hom coefficient pair to 24 exact Eisenstein vectors.  Reduction
   modulo `397` is then injective on that set and may identify the vector.
   The characteristic-zero confinement precedes the modular computation.
3. **Differentials are injective on Hom classes in characteristic zero.**
   They do not recover translations of curve-map representatives.  The
   divisor observer `[tau Q-Q]` cancels those translations exactly.

These are different source objects, kernels, and sidecars.  The invalid move

```text
specialized value is zero  =>  every ambient coefficient is zero
```

does not become valid merely because a finite-field calculation is exact.
The valid THM-4259 move is instead

```text
characteristic-zero theorem gives a finite shell
  -> reduction is injective on that shell
  -> a nonsingular reduced system selects one shell element.
```

That distinction is the strongest conceptual bridge between the p-adic audit
and the new Jacobian calculation.

## Confluence replaces a twelve-step orbit by an initial-value problem

On the `a_u=0` residual Hom lane, THM-4258 imports the exact annihilator

```text
P_0(T)=(T-omega^2)(T^2+omega).
```

For consecutive attachment differences `delta_j`, this gives

```text
delta_(j+3)-omega^2 delta_(j+2)
  +omega delta_(j+1)-delta_j=O.
```

Thus three consecutive zero differences suffice for all twelve to vanish.
The gain is not merely a factor-four speedup from `18,144` to `4,536` group
values.  It changes the proof object from a cycle-wide equality test to a
confluent initial-value certificate.  The equality boundary is sharp at the
ambient recurrence level: two samples have four-to-one fibres in the
conditional `E_0[2]~=F_4` sidecar, while three samples determine one of the
64 abstract patterns.  At most 64 of those patterns need be geometrically
realized; surjectivity was never claimed.

The p-adic analogue is the complete Hasse jet, not a repeated scalar
specialization.  In both cases a recurrence or degree cap proves that a
finite observer window is complete.  Without that cap, more samples may
still leave a nontrivial kernel.

THM-4264 subsequently uses the actual visible incidence ideal rather than
the whole ambient recurrence module.  Its visible projector lowers the
nilpotence from `N^3=0` to `N^2=0`, so two zero differences suffice on every
one of the `1,512` inherited incidences.  Exactly `16` projected words remain,
with period profile `1:1,3:3,6:12`; a one-edge hostile proves sharpness.
This is a second illustration of the same rule: an observer shortens only
after a lawful sidecar removes the lost eigenspace.

## The half-sum obstruction was a basis problem, not a missing map

THM-4241 had already constructed the degree-four mixed map `H_lambda` and
proved that the visible-hidden direct sum has index four in the full Hom
lattice.  What remained unknown was the exact change of basis from this map
to the census generator

```text
2h=v+omega^2 f+g.
```

THM-4259 now fixes one algebraic branch and proves

```text
H_lambda+H_lambda o iota=-v,
H_lambda-H_lambda o iota=omega^2 f+(omega^2-omega)g,
h=H_lambda+v-omega^2g,
Th=omega^2h-omega f.
```

The first identity comes from exact elliptic addition.  The second comes from
the finite-shell mechanism above.  Consequently every residual row
`m=bf+cg+dh` has the executable Hom formula

```text
m=dH_lambda+dv+bf+(c-domega^2)g.
```

There is still no canonical literal translation for every curve-map
representative.  None is needed for
`F_m(Q)=m([tau Q-Q])`, because it evaluates representatives by point
differences.  This is a useful design rule: normalize only the quotient data
that the final observer preserves.

## Local/nonlocal pivot language: useful architecture, no theorem transfer

Cerf's arXiv `2608.23661v1` is a percolation paper.  Its local/nonlocal
higher-pivot split suggests a proof architecture: a local observer closes
only after a nonlocal residual has been separately controlled.  In the
present Jacobian lane, the three-sample recurrence is one local observer;
THM-4260 controls the full gate-interior residual instead through an
orthogonal reciprocal-denominator obstruction.  In the p-adic lane, Hasse
jets are local normal data; an unrestricted section kernel is the residual.

This is only an architectural analogy.  Percolation pivots are Boolean events
in random spatial configurations.  Hom-lattice recurrences and completed
coefficient modules are deterministic algebraic objects.  There is no
probability-preserving map, no shared measure, and no imported exponent or
density conclusion.

## Typed connection table

| source | target | map | preserved predicate | lost data | required sidecar |
|---|---|---|---|---|---|
| completed bivariate series | one formal arc | `u->g(f)` | specialized value | ideal `(u-g)`, coefficientwise order | Hasse jets or degree cap |
| exact hidden degree-12 shell | two residues in `F_397` | coefficient reduction | identity of one of 24 candidates | unrestricted algebraic coefficients | coercive norm bound and injectivity on shell |
| curve morphism | Jacobian Hom class | Albanese functor | degree-zero divisor values | target translation | evaluate `[tau Q-Q]` |
| twelve cyclic differences | three initial differences | cubic recurrence | all-zero orbit predicate | no orbit data once recurrence is fixed | annihilator `P_0(T)` |
| abstract `F_4` recurrence | geometric attachment orbit | realization map, not needed by THM-4260 | necessary two-torsion pattern | geometric incidence constraints | explicit map-ratio table for an orthogonal replay |

## Incoming closure and the new sharp front

THM-4260 landed concurrently after the observer and dictionary were proved.
It does not evaluate the `4,536` observer values.  Instead, for every residual
map class it applies the necessary anti-invariant projection, proves the
sharp characteristic-zero denominator bound `deg(d_ell)<=4K-1`, attains that
bound at a good prime, and uses a monic local-DVR argument to lift

```text
gcd(d_ell(t),t^deg(d_ell)d_ell(-1/t))=t^2-1.
```

Those two roots are exactly the excluded `Z=0` wall.  Therefore all `280`
classes and `1,512` incidences are impossible on the already-specialized
gate-interior `W=0` fibre.  The explicit THM-4259 dictionary is an orthogonal
coordinate certificate, not a dependency of this closure.

The finite-exact THM-4265 scout sharpens the wall geometry without yet
promoting a theorem: at both `397` and `577`, every reduced denominator has
nonzero derivative at both `t=+1` and `t=-1`.  Thus the unavoidable wall
factor is reduced in every one of the 280 audited rows.  What remains unknown
is the transverse coefficient: the natural local test is the determinant of
the `t` and `W` derivatives of a denominator and its reciprocal at
`(W,t)=(0,+/-1)`, but its `W` column is not defined until the off-fibre
deformation dictionary exists.

The first transverse obstacle is type-theoretic.  The basis `u,f,g,h`, its
`280`-class lattice, and `d_ell(t)` currently exist only after setting
`W=0`.  There is no proved flat off-fibre Hom basis or
specialization-compatible representative system.  A formal expression such
as “the `W` derivative of `d_ell`” is therefore not yet defined.

The next proof-producing object is a deformation or boundary certificate:

1. construct the bounded-degree Hom/incidence scheme over the full
   exact-`M=12` base and determine whether it is proper near the gate-interior
   `W=0` fibre; properness plus the empty special fibre would already give a
   neighbourhood exclusion;
2. if properness fails, construct a flat local deformation dictionary before
   forming any denominator/resultant packet, and retain its exact `W`-order;
3. only then compute the first nonzero transverse Hasse coefficient,
   explicitly quotienting the specialization kernel identified by THM-4255,
   and test whether a common factor can enter from infinity;
4. re-enter the homogeneous pre-division equations on `U=0`, `Z=0`, and
   `U+Z=0`; in the `t` chart the corresponding limits are `0/infinity`,
   `+/-1`, and `t^2=-1`, so division by `t`, `t^2-1`, or `t^2+1` is forbidden
   on the matching wall;
5. retain seam-entry hypotheses: exclusion on a specialized fibre is not
   entry into the full exact-`M=12` gate.

This is where the p-adic specialization lesson becomes load-bearing for the
Jacobian frontier.  `W=0` emptiness cannot be read coefficientwise in `W`.
A valid neighbourhood theorem needs a transverse jet, a Weierstrass/monic
factor argument, or an exact resultant valuation.  Neither the new fibre
closure nor the explicit dictionary by itself proves seam entry or `JC(2)`.
