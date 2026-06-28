# LRC14 Equioscillation Contact Graph Case Split

The AP safety profile touching `1/14` at exactly six points should be treated
as a contact graph, not only as a plotted curve.

For AP and Goddyn-Wong, the six contacts are the units
`{1,3,5,9,11,13}/14`.  Each unit `a` is bound by the complement pair
`{+a^{-1},-a^{-1}} mod 14`, so the six unit contacts collapse under
`t -> 1-t` to three binder slots:

```text
{1,13} -> {1,13}
{3,11} -> {5,9}
{5,9}  -> {3,11}
```

That is the graph hiding in the equioscillation picture: three `K_{1,2}`
contact stars, or a perfect matching after antipodal quotient.  The topology
view says these are zero-measure holes of the open danger cover.  The
geometry view says they are active constraints of a max-min saddle.  The
graph view says the active constraints have only three independent
complement-pair slots.

The strongest immediate payoff is the case split.  At a unit `a`, a row has
`f_S(a/14) < 1/14` exactly when it contains a speed `0 mod 14`.  Thus every
14-free row is already safe at all six unit times, while every row with a
multiple of 14 kills the whole unit-contact graph and must find its witness
off the units.  Tightness means the contact graph survives globally and no
off-unit peak rises above it.

This sharpens HYP-3243.  The finite chamber theorem should not merely ask
whether there is a safe interval.  It should classify how the contact graph
behaves:

```text
survives and is global -> AP/GW equality core
survives but is not global -> strict open row
is killed by 0 mod 14 -> promoted Phi_{14d} / covering-floor route
does not fit -> finite chamber discharge or named residual debt
```

The non-geometric proof routes are projections of this contact graph.
Kolmogorov/Chebyshev sees active-gradient convex hulls.  Borsuk-Ulam sees the
odd index of the three antipodal pairs.  The Cech route sees boundary holes
of the danger cover.  Green/Toeplitz and root-motion routes become cell labels
that discharge rows after the contact graph is classified.

The post-rebase HYP-3250/HYP-3251/HYP-3252 evidence makes the burden cleaner.
HYP-3250 supports the finite tight-locus plus uniform-margin shape that this
chamber classifier wants.  HYP-3251 and HYP-3252 warn that the
index/Gauss-sum data is ambient: it explains the AP saddle and the three
antipodal contact pairs, but the `S`-dependent work still has to come from the
covering floor, decorrelation margin, or finite discharge when the contact
graph is killed or non-global.

The later S81/S82 and Q(sqrt(-7)) Thread A/B updates fit this reading.  They
strengthen the unit/dilation construction and bounded-margin side, reframe
resonant multiples of `14` as off-grid bulk after an on-grid core hit, and
separate the residue layer from the magnitude-level census.  They do not turn
the ambient index into the proof.  The contact graph remains the
equality-branch carrier plus a dispatcher to the floor/margin machinery.
S83 adds a useful symmetry refinement: the three binding-pair slots form one
`C_3`/real-cyclotomic Galois orbit, so a proof can aim for one local
binding-pair lemma plus equivariant transport around the contact graph.
HYP-3260 gives the complementary nullspace audit: the contact graph is the
right equality-branch carrier, but its unit active gradients have rank `3`
and forget the nonunit residue/height and covering sidecars.
HYP-3300 turns that into the next proof-program shape: make the contact graph
and nullspace sidecars columns in an observability matrix, then use them as
Morse boundary data in the finite chamber descent.

The challenged assumption is that runners are the natural graph vertices.
For this route the vertices are unit contacts, antipodal pairs, binder pairs,
boundary holes, covering kill-switches, and proof obligations.  A scalar
safety value is only a shadow; the proof-grade object is the contact graph
plus the sidecars that say whether it survives, opens, or is killed.

-> HYP-3265, HYP-3300, HYP-3260, HYP-3259, HYP-3258, HYP-3256, HYP-3255, HYP-3253, HYP-3252, HYP-3251, HYP-3250, HYP-3248, HYP-3246, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3238,
HYP-3237, HYP-3236, HYP-3218, HYP-3214, HYP-3132, HYP-2928, HYP-2909,
THM-523, THM-530, T1355, LTI-355, LTT-255, OPEN-Q-108.
