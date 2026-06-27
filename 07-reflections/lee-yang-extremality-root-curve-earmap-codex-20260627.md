# Lee-Yang Extremality, Root Curves, And Ear Maps

The prior PGF-zero signal was right, but still too scalar.  Counting real roots
is the first useful shadow; the whole root curve is the object.

The anchored scout makes that concrete.  In the `C(13,7)=1716` bounded bank,
`consec_8` is the unique `p0` and `L_yK8` leader, and all high rows live in the
`#real=0` stratum.  The correlations are strong enough to matter:
`corr(p0,min|z|)=+0.899`, `corr(p0,dist(roots,[-1,0]))=+0.681`, and
`corr(p0,-#real)=+0.483`.  So the Lee-Yang reading is not cosmetic: roots
near the negative segment are genuinely low-coverage warnings.

The caveat is equally important.  Nearest-root radius alone is not the theorem:
the row `(0,2,4,6,7,8,10,12)` beats `consec_8` on nearest root radius but not
on `p0`.  The proof signal has to retain the full root locus, including root
product, angles, and zero-free segment clearance.

The phi4 cue gives a new diagnostic.  Symmetrizing around `N=3` and fitting
`V(s)=c+b*s^2+lambda*s^4` makes `consec_8` a clean double-well:
`b=-0.10232`, `lambda=+0.00566`, well radius about `3.005`.  That says the
sector-empty law has boundary/extreme preference, exactly what gK8 rewards.
But by `k=9..13` the naive quartic fit destabilizes, so the honest next model
is a moving-center or sextic potential, not a forced phi4 story.

The most useful surprise is graph-theoretic.  The anchored `#real=0` stratum is
a single one-swap component with `290` vertices, `1048` edges, cycle rank `759`,
and no articulation points.  Bidirect the edges and it is strongly connected,
so it has an ear decomposition.  The proof question is sharper than that
theorem: can the whole zero-real component be generated from a small AP/consec
seed by labelled root-deformation ears before crossing the root-collision wall
into `#real=2`?

That is where the user prompt's ear hints should land.  Strong connectivity is
the zero-real plateau.  Odd-ear/factor-critical structure should be tested on
minimal collision cycles at the `#real=0` to `#real=2` boundary.  Nested ears
should be tested after quotienting AP symmetries; if present, they would
explain how exchange optimality can be non-transitive but still finite and
controllable.

The two useful maps are the HYP-3043 lens map and the HYP-3104/HYP-3106
signal/functor maps.  Together they say when the PGF root sidecar is legally
forgotten.  The answer after this scout is: almost never before the root
collision coordinate has been discharged.
