# LRC14 Green-Current Conductance Graphs

HYP-3223 was right to ask for effective resistance, but the scout shows the
right unit is a graph of proof obligations, not just the electrical graph on
sectors.

Push-time rebase clarified the numbering.  Incoming HYP-3225 is the local
Green-current / Lorentzian trap-fingerprint scout on the HYP-3202 trap
neighborhood; incoming HYP-3226 reserves the small-pattern motif atlas.  This
reflection, now HYP-3227, is the full-bank conductance graph layer.  The three
fit together: HYP-3225 names trap sidecar types, HYP-3227 measures how those
trap debts connect to certificate coordinates, and HYP-3226 should treat both
as payload atoms rather than loose pattern numerology.

On the sector side, the news is clean but not sufficient.  Consecutive and
the doubled AP dilation have no beaters for positive-covariance total
conductance, algebraic connectivity, minimum degree, Kirchhoff index, and the
corresponding precision-graph lambda2/Kirchhoff/killing metrics.  That is a
real extension of HYP-3202 and HYP-3224: the covariance-layer dominance is
visible as electrical capacity, not just as three sums.

The failure mode is equally useful.  The precision M-matrix defect has `181`
primitive beaters, and precision minimum degree has `3`.  So the proof cannot
say "invert the covariance matrix and AP is best in every Green-coordinate."
The inverse kernel is a grounded conductance sidecar with a named positive
off-diagonal defect.  That matches the old repo warning from the raw
conductance route: conductance is a signal, but scalarized conductance is the
wrong theorem.

The better object is the trap-discharge graph.  HYP-3224 made every exchange
trap discharge through Toeplitz lambda-min under its six-coordinate normal-fan
weighting.  HYP-3227 changes the question: if the traps are vertices and the
coordinates are terminals, does the debt graph stay connected when we remove
or restrict coordinates?  It does.

```text
all-coordinate trap graph: lambda2 = 2.719948208
no-Toeplitz trap graph:   lambda2 = 2.537866286
green-only trap graph:    lambda2 = 1.208613477
```

That is the important proof-frontier improvement.  Toeplitz remains the
strongest moment-cone chart, but it is not carrying the finite trap manifold
alone.  The Green-only graph connects every non-AP trap by
`cov_positive_lambda2`, `cov_positive_kirchhoff`, `negative_covariance_debt`,
`precision_lambda2`, `precision_kirchhoff`, or `precision_mmatrix_defect`.
This makes the electrical language more than metaphor: it names which
networks leak, which networks bottleneck, and which finite rows are in the
same obstruction island.

The Fiedler split is especially suggestive.  Four traps sit with
`negative_covariance_debt` and `precision_mmatrix_defect`:

```text
(0,2,4,6,7,8,10,12)
(0,1,2,3,7,8,9,10)
(0,2,5,7,9,10,12,14)
(0,1,4,5,7,8,11,12)
```

This looks like the first finite Schur-complement subcase.  The remaining
traps stay on the covariance/precision-capacity side.  A reasonable next
lemma is not "all traps have the same bottleneck"; it is "the trap manifold
has a small conductance cut, and each side has a different discharge rule."

The proof route I would push next is:

```text
1. Model C_E^{-1} as a grounded response matrix Q_E = L_E + K_E + D_E.
2. Show legal exchange moves are Schur-complement or star-mesh edits of Q_E.
3. Prove AP is the only row where all capacity coordinates are tight and
   all named defect coordinates vanish or are sidecar-equivalent.
4. For the 11 traps, use the trap graph as the finite certificate ledger:
   each trap has positive conductance to at least one boundary coordinate.
```

The algebraic-connectivity idea is strongest after this translation.  Lambda2
on the six-sector graph is a good AP-tight capacity coordinate.  Lambda2 on
the trap graph is a finite Helly/connectivity statement: the trap debts are
not isolated anomalies; they are connected to the same certificate boundary
that HYP-3205 and HYP-3224 already identified.

This also clarifies how to synthesize the recent incoming work.  HYP-3222
gives the Joukowski/Hermite-Biehler/Perron exact spectral sibling.  HYP-3224
gives the normal-fan payload cube.  HYP-3227 supplies an electrical
interpretation of two pieces of that cube:

```text
covariance layers     = conductance channels
Toeplitz/Perron side  = grounded response and spectral capacity
exchange traps        = finite bottleneck networks
AP support direction  = terminal normal for the discharge graph
M-matrix defect       = named precision-side leakage, not a scalar proof
```

The later HYP-3214 Fejer-kernel magic-function packet gives a sharper
candidate for the dual object above this graph: if the Fejer/Delsarte slack
dominates the Green-only trap weights, the conductance graph becomes a finite
electrical shadow of the same magic-function certificate rather than a
separate sidecar.

The creative extension is to use algebraic connectivity twice: once in the
usual graph-theoretic sense on sectors, and once as a proof-obligation
connectivity invariant on traps and certificate coordinates.  The second use
is more likely to become formal, because it can be reduced to a finite
weighted inequality ledger over the `11` non-AP traps.

Next tests:

1. Emit exact characteristic-polynomial certificates for the AP-tight
   `cov_positive_lambda2` and `precision_lambda2` comparisons, at least on
   the trap rows.
2. Compute Schur-complement reduction words for the four Fiedler-positive
   defect traps and compare them to the remaining seven traps.
3. Add the Green coordinates to the HYP-3205 dictionary and test whether the
   AP/dilation two-row skyline survives all coordinates that have no primitive
   beaters.
4. Search for the HYP-3214 Fejer/Toeplitz/Verblunsky slack that dominates the
   Green-only trap graph weights, not just the Toeplitz gap.
5. Formalize the finite trap graph as a Helly-style lemma: every non-AP trap
   has positive normalized deficit in a connected certificate cover after
   deleting any one nonessential coordinate.

-> HYP-3227, HYP-3226, HYP-3225, HYP-3224, HYP-3223, HYP-3222, HYP-3221,
HYP-3214, HYP-3213, HYP-3212,
HYP-3211, HYP-3210, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3201,
HYP-3200, HYP-3163, HYP-3162, HYP-3161, HYP-3160, THM-577, T1325, LTI-325,
LTT-225, OPEN-Q-108.
