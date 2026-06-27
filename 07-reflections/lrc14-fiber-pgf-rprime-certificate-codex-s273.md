# LRC14 Fiber-PGF Rprime Certificate

The useful new move is to stop treating `Rprime` as either a raw scalar or a
bare Fourier covariance.  For a residual row `S=R union 14Q`, pass to
`u=14t`.  The Q-side only sees `u`; the R-side sees all 14 lifts
`(u+a)/14`.  Therefore the right generating function is the sheet-count PGF

```text
F_R(y) = sum_n meas{u: N_R(u)=n} y^n,
N_R(u) = #{a: (u+a)/14 is R-safe}.
```

The exact identity is:

```text
Rprime = E[N_R | Q-lonely] / E[N_R].
```

This is the coefficient-level child of HYP-3136: that note factors the
multi-far floor as `L=Rprime*meas(R-safe)*meas(Q-lonely)` and leaves the
`Rprime` factor to the signed-SPEC constant chase.  The fiber PGF packet gives
that factor a finite generating-function object before any scalarization.

This is a better proof object because it is finite, exact, and keeps the
information that the Fourier/SPEC transform uses.  It also respects the
controlled-forgetting pattern from HYP-3134: individual sheets may be forgotten
only after their coefficient packet is named.

The first hard row is now sharply visible:

```text
R={1,...,12}, Q={1,2}
F_R(y)   = 7243/13860*y^0 + 6617/13860*y^1
F_R,Q(y) = 7243/13860*y^0 + 521/1980*y^1
Rprime   = 51058/72787 = 0.701471...
```

So a pointwise theorem like "Q-safe always has a positive R-safe lift" is
false.  The right statement is a conditional first-moment inequality: Q may
include zero-sheet mass, but it cannot depress the first derivative of the
fiber PGF too much.

This makes several older niche threads operational rather than decorative:

- Lee-Yang PGFs supply the coefficient/root language.
- Delsarte/MacWilliams weight enumerators supply the transform viewpoint.
- q-Pochhammer/modular-tail work says tails need a named principal packet.
- Moser/fibbinary partial cubes say finite coordinate words are legal only as
  typed sidecars, not scalar shortcuts.
- A000568 edge-envelope work supplies the quotient discipline.

The next proof-facing task is to compute the finite legal residual packet
family and prove

```text
F_R,Q'(1) / F_R,Q(1) >= c * F_R'(1) / F_R(1)
```

before transforming back to HYP-3129's signed SPEC language.
