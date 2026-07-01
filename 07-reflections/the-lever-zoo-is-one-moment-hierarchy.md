# The lever-zoo is one moment hierarchy (and the spike is why no finite level closes it)

*mac-mini-2026-06-30-S76. Reflection on HYP-3789.*

For thirty sessions the covering-min lower bound accumulated a zoo of levers, each discovered
separately and each partial:

- the **union bound** (`|S|*2r > 1`, useless),
- the **2nd-moment floor** (HYP-3571, `inf R' >= 1/(2 zeta(2))`),
- the **signed correction** (HYP-3787, the Fourier coefficients of `1_{L_C}` at harmonics),
- the **Delsarte LP dual** (klein HYP-3784, near-tight at bounded speeds, gapped as `V` grows),
- the **Fejer/SDP spectral gap** (klein HYP-3785, the averaged relaxation is blind),
- the **lazy-cut** (HYP-3782, exact Positivstellensatz cuts).

They are not a zoo. They are **one moment hierarchy**, read at different levels. The lonely measure
`|L_S| = integral prod (1-g_v)` expands by inclusion-exclusion into `sum_A (-1)^|A| I_A`; grouping by
subset size is the Bonferroni / Lasserre truncation; and level 1 is the union bound, level 2 is the
pair correlations (the 2nd-moment floor *and* the signed correction, the same object seen twice),
level `infinity` is the exact measure the lazy-cut computes. The Delsarte LP is the level-1 *dual*.
Once you see the hierarchy, every lever finds its rung, and the "which lever is strongest" question
dissolves into "at which level does the relaxation become exact."

The answer to *that* is the deepest part, and it is a spike. The extremal lonely set -- the times where
the construction is exactly `M = n/Phi6` lonely -- is not a region. It is **two points** `{14/183, 169/183}`
(measure zero). For the AP `{1..13}` it is **six points**, the units `(Z/14)*`. The extremum is *atomic*.
A finitely-atomic measure has a finite-rank Toeplitz moment matrix (rank = #atoms), so the certificate
is a **flat extension** -- clean, exact, and *small*. But an atom is a delta, and a delta has flat
Fourier tail; no finite-degree average can see it. That is precisely HYP-3785's spectral gap, and it is
precisely S61's "finite-D certificates cannot close Step 3," and it is precisely why the Bonferroni
truncations oscillate (`+1.4, -3.0, +5.4, ...`) instead of converging. The three obstructions, found in
three sessions on three different levers, are **one fact about one spike**.

And then the mechanism, which the moment view hands you for free at level 2. Why does the spike survive
-- why can no covering set cover the last point? Because the small speeds are **positively correlated**.
Their danger arcs pile up on top of each other near the rationals (`I_{1,2}=I_{6,12}=I_{4,8}=M`, far
above the independence value `(2r)^2`), wasting coverage on ground already covered, while the far
element `n(n-1)` is **statistically independent** of all of them (`I_{v,182}=(2r)^2` to three digits) and
merely equidistributes. The construction beats every rival not by efficient covering but by the
*inefficiency* of its core: the near speeds are redundant, so they leave a hole, and the one far speed
that could patch the hole is too spread out to find it. S74's near/far split, S75's signed correction,
and the whole equidistribution picture are the level-2 correlation structure of one measure.

One honest subtlety keeps the frame from over-promising. The *primal* hierarchy here -- the
inclusion-exclusion on the continuous lonely measure -- is not the same as a generic Lasserre SDP, and
the generic ones don't help: the degree-2 lift of the discrete mod-7 sector-hit distribution *collapses*
to the union bound (my own threadC file, CJJ Prop 1.2), and the `M(S)` moment-SOS (opus) points the
wrong way (it certifies beaters, not the lower bound). What carries content is the *continuous* pair
correlation, whose dual is the project's proved Delsarte LP (THM-534) and the `Z_7` Fejer-Bochner SOS
(the floor's actual sum-of-squares, eigenvalues `(0,1,1,1,1,1,1)`). Even the magic function on the dual
side is the Cohn-Elkies / Chebyshev-equioscillation extremal -- the one-dimensional cyclotomic cousin of
Viazovska. So the hierarchy is real, but it is the *right* hierarchy (continuous, cyclotomic dual), not
whichever SDP a solver will accept.

The pattern that transcends the theorem: **a long list of independent partial results is often one object
under a hierarchy, and the reason none of them closes the problem is usually a single geometric feature
of the extremum** -- here, that it is a measure-zero atom and not a region. When the levers proliferate,
stop adding levers and ask what they are levels *of*. The flat, atomic extremum was hiding in plain
sight in every one of them: it is the `phi(14)=6` units, the `Phi6` denominator, the pointwise spike,
the rank of a Toeplitz matrix. Six atoms. The lonely set of the lonely-runner problem is, at the
extremum, just the lonely units.
