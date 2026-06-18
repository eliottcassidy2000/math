# LRC14 Colored Resonance Discrepancy

Codex 2026-06-18.

The component bound was honest but dumb: every micro-interval endpoint paid a
unit of discrepancy.  The new question was whether those endpoints are really
independent.  They are not.

The useful object is the color-summed quadrature error.  With

`F_b(x)=prod_p h(px) prod_e h(ex-b/14)`,

the colored count samples `F_b` at `x=t/V+b/(14V)`.  After Fourier expansion,
only terms satisfying both

`sum_p p*n_p + sum_e e*n_e = ell*V`

and

`sum_e n_e == ell mod 14`

survive.  That second congruence is the color gate.  It is the algebraic
reason the thousands of phase endpoints are mostly ghosts.

Boundary values matter.  The clean Fourier identity uses the half-boundary
indicator.  The actual CRT count is closed, so binding-boundary residues add a
bonus.  In one row the open deficit was `75.457`, but the closed actual
deficit was `-44.543` because `120` boundary hits came back.  That is not a
rounding nuisance; it is part of the arithmetic.

## What Survived

A universal constant `100` is false.  The row

`P=(1,2,11)`,
`E=(0,84,293,301,355,416,485,665,843,886)`,
`V=1203`

has actual deficit about `135.435`.

But the additive bound

`Delta <= 8*(k+cGP)+1`

survived the deterministic rows and the reproducible stress bank.  This is
much more plausible than the old product-style or raw-`K` bounds.  It says the
cost is cluster size plus the number of coarse small-speed safe components,
not the number of micro-components created by the large offsets.

If true, it is strong.  The maximum possible `cGP` for `P⊆{1,...,13}` is `32`;
with `k<=13`, the deficit is at most `361`.  Together with the observed
`Sigma` floor `14249/28028`, that gives a large-`V` cutoff of `711`.

So the remaining proof would become:

1. Prove the uniform colored mass floor.
2. Prove the resonance discrepancy lemma.
3. Check `V<711`.

That is still real work, but it is a finite-looking proof architecture now.

## The Next Lemma

The lemma I would try to prove is:

For every coarse component `I` of `G_P`, the color-summed quadrature error of
`I cap C_b(E)` is bounded by `O(k+1)`, after closed-boundary correction.

Summing over the coarse components naively gives `O(k*cGP)`, which is already
useful.  The data hints at the stronger additive `O(k+cGP)` form, which would
require cancellation between coarse components too.  The Fourier condition is
the likely way to see that cancellation cleanly.

The tournament piece was useful mostly as a discipline check: candidate bounds
formed a transitive proof-obligation tournament.  Raw `K` wins robustness by
being huge; among usable bounds, `8(k+cGP)+1` is the current best target.
