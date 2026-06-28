# LRC14 Covariance Laplacian And Associator-Ear Angles

This pass tries to turn the current dense Lee-Yang/Worpitzky/quartic web into
two proof tasks that are less scalar and more checkable.  Incoming HYP-3161
changed the tone: consecutive already maximizes total covariance over all
`3432` bounded k=8 clusters, and the sign change at `k=5 -> 6` makes the
ferromagnetic mechanism concrete.  So the right next question is certificate
shape, not more confirmation.

The later HYP-3153 scout gives this note a concrete input table rather than a
wish list.  It supplies the exact `L_y <= cap` margins, the k=8
`-9S3+6S4` split, the Worpitzky `-1/3` local mode, the biquadratic fold, and
the odd-ear witness grammar.  The next scout should join those columns to the
covariance-kernel and associator-cocycle columns proposed here.

Final rebase also made HYP-3162 available: the Joukowski/de Moivre sidecar is
not just root-curve vocabulary, but a cubic 7th-cyclotomic calibration for why
`n=14` is the first open apex after the degree `1,2` cases.  That changes the
scout slightly: preserve `cyclotomic_cubic_defect`, `de_moivre_angle_slack`,
and the ferromagnetic/antiferromagnetic bridge status while decomposing rows
into covariance slack and associator debt.

The later HYP-3200 bounded-bank audit makes the proof target cleaner again.
Work in primitive normal form: consecutive is the exact `Sigma kappa_2` leader
there, while the all-bank exception is the dilation twin.  It also kills the
temptation to prove a `1/7` law and downgrades `kappa_4` to a stabilizer
sidecar.  Those are bookkeeping fields, not the theorem.

HYP-3201 adds a useful correction to the information-theory language.  Row
entropy is not the extremal quantity, but conditional entropy of a target
after a quotient is a good illegality detector.  Add `H(target|quotient)` to
the scout as a sidecar: zero means the compression is legal; positive residual
names the ordered, bracket, action, root, or cumulant payload still owed.

KPS-S31al adds the strongest concrete spectral formulation so far: treat
`q_t` as a trigonometric moment sequence and maximize the minimum eigenvalue
of its Toeplitz matrix.  Consecutive reportedly wins this `lambda_min(T)` test
over all `3432` bounded k=8 rows, which makes Szego/Caratheodory-Fejer tools a
serious proof route.  Its Griffiths lead also warns that naive speed-path
monotonicity fails; the ferromagnetic route needs a random-current or
coupling-manifold order.

The first angle is to stop treating HYP-3160/HYP-3161's total covariance as
only a number.  The covariance matrix of empty-sector indicators should be
measured as a kernel.  If the consecutive row is the ferromagnetic extremal
state, then the difference between its kernel and any competitor should show
up as a PSD, Laplacian, Monge, or conditionally positive type certificate
after the known sidecars are retained.  This would convert the even face into
a degree-two rearrangement problem.

The second angle is to let the odd Worpitzky residue stay odd.  Forcing it
into the covariance proof is probably the wrong compression.  The better
object is a third-cumulant/associator cocycle whose boundary is visible in
the n=3 minority-edge kernel and whose recursive witnesses are odd ears.
Then the target becomes an exchange inequality: odd associator debt is bounded
by even covariance slack plus named finite sidecar debt.

The most important practical detail is the vertex choice.  Runners, arcs, and
score classes are not the right first vertices here.  The tournament should
rank proof obligations and sector-pair/triple certificates by what LRC
predicate they preserve.  Raw `L_y` comes last, because HYP-3153, HYP-3154,
HYP-3161, and HYP-3199 have already made clear that the scalar only becomes
proof-facing after root curves, real-rootedness-defect data, canary/filler
status, minority-edge data, and sidecar exits are named.

Next computation: emit `covariance_kernel_distance_profile`,
`monge_four_point_defect`, `PSD_dual_slack_vector`,
`associator_triple_cocycle`, `odd_ear_surplus`, and
`even_covariance_odd_associator_exchange` over the same bounded k=8 banks used
by HYP-3138/HYP-3139/HYP-3142/HYP-3160.  The desired outcome is not another
confirmation that consecutive wins.  It is a decomposition of every
competitor into even kernel slack and odd associator/ear debt.
