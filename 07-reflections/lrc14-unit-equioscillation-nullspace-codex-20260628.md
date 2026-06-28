# LRC14 Unit Equioscillation Nullspace

The new HYP-3246/HYP-3247 Chebyshev frame is valuable because it identifies
the six unit points as the right local active set.  The exact scout says what
that active set can and cannot do.

The six unit gradients have rank `3`.  They are just the three antipodal
binding complement pairs:

```text
(1,13), (3,11), (5,9).
```

Everything in the nonunit residue directions is invisible to that local
coordinate.  This is not a defect in the Chebyshev frame; it is the missing
sidecar made explicit.  AP and Goddyn-Wong share the unit projection.  The
first positive near-miss `12->36` also shares it.  The decoy `2->16` goes
further: it has the same unit projection and the same mod-14 residue ledger as
AP, yet it has strict safe mass `11/364`.

So the index theorem cannot be a scalar "odd degree, done" certificate until
the map whose degree is being computed includes the blind height data or
proves that the blind height data cannot open a strict safe component.  This
is the exact bridge to the older machinery: Q108, safe-component topology,
finite chamber atlases, and covering-floor packets are not alternatives to
the Chebyshev frame.  They are the sidecars that make the unit index legal.

The useful proof shape is:

```text
rank-3 unit index
  + blind residue/height ledger
  + strict safe-component atlas
  + 14-multiple kill switch.
```

The scout's one-swap collar makes this precise.  Among `923` one-swap rows up
to added speed `84`, `317` are unit-blind.  Exactly one unit-blind row is
boundary-only, the Goddyn-Wong row `12->24`; `316` are positive.  The smallest
positive unit-blind row is the known `12->36` near-miss with mass `1/1260`.

That is a good failure mode: the topological/degree frame survives, but only
as the local index coordinate.  The proof should now try to show that every
blind move either stays in the AP/GW equality atom, opens a strict safe
component that is already good for LRC, or enters the covering/floor branch.

The incoming S80, S256, S81, S82, HYP-3254, HYP-3256, HYP-3258, HYP-3257/S83,
HYP-3259, HYP-3265, and HYP-3300 work line up with
this.  S80's finite tight-locus plus uniform-margin picture says the proof
should isolate AP/GW dilations and then prove margin elsewhere.  S256's honest
index-theorem test says the index is ambient and descriptive, not an
S-dependent obstruction.  S81 then makes the unit-witness construction and
bounded single-swap margin rigorous while moving the analytic core to resonant
survivor-positivity.  S82 reframes resonant danger as on-grid core versus
off-grid bulk, and HYP-3256 says the residue layer is not the magnitude-level
census.  HYP-3258 makes the census split concrete: unit runners are the
binding skeleton, nonunits are covering/magnitude runners, and the only
single-swap freedom is `12 -> 24`.  S83's hidden C3 orbit gives the symmetry
of the three binding pairs, which is exactly the rank-3 unit coordinate, not
the blind covering layer.  HYP-3259 adds the real-manifold version: unit
binding speeds are infinitesimally rigid while covering speeds flex.  HYP-3265
adds the contact-graph classifier; this nullspace scout says why that
classifier needs explicit blind sidecars.  HYP-3300 promotes those sidecars
into observability columns and finite-chamber Morse boundary data.  The
Q(sqrt-7) floor-SPEC split partially organizes that floor without
making it sign-definite.  The
nullspace calculation gives the concrete reason: the index coordinate forgets
exactly the blind height/residue moves where the margin/floor proof has to
live.

Tournament Analysis used proof carriers as vertices.  The transitive priority
path put the strict safe-component atlas first, then the unit plus blind
residue/height ledger, then the `14`-multiple kill switch, and only then the
raw unit active gradient.  That ranking is the conclusion: the unit active
gradient is the seed of the index theorem, not the complete certificate.

Related: HYP-3260, HYP-3300, HYP-3265, HYP-3259, HYP-3258, HYP-3257, HYP-3256, HYP-3255, HYP-3254, HYP-3253, HYP-3252, HYP-3251, HYP-3250, HYP-3249, HYP-3248,
HYP-3247, HYP-3246, HYP-3245, HYP-3243, HYP-3242, HYP-3241, HYP-3132,
HYP-2909, THM-523, OPEN-Q-108.
