# The seven-sector cover is a moment problem, and its dual is an integer-root Bonferroni polynomial

*mac-mini-2026-06-18-S7 — Angle D (LP / SOS dual certificate) for LRC(14)*

## The reframe

The seven-sector reduction (HYP-2603/THM-532) asks to bound `meas(S7(E)) ≤ cap_k`, where
`S7(E)` is the set of `x` at which all seven length-1/7 sectors are hit by the cluster
phases `{frac(e_i x)}`. THM-532 tried to write `meas(S7) = M7(k) + corr(E)` and bound
the relation-lattice tail `corr` by `C·W(E)`. That product bound is the wrong shape — the
tail is not a small perturbation (at consec it is 0.30 vs a main term M7(8)=0.024), and
the crude absolute bound `C*·W(consec_8)=0.384` overshoots the margin `0.357`.

The right object was hiding in plain sight. Let `N(x)` = the number of UNHIT sectors
among `{1,…,6}` (sector 0 is always hit by `e=0`). Then `S7` is exactly `{N=0}`, and

    meas(S7) = p_0,    p_t = meas{N = t},    t = 0,…,6.

So `meas(S7)` is the value-at-0 of a **probability distribution on {0,…,6}**, and the
only data we can compute exactly and uniformly are its **factorial moments**
`S_r = E[C(N,r)] = Σ_{|A|=r} J(A,E)` (J = multi-sector avoidance measure). This is a
**truncated moment problem**: bound `p_0` from the moments `S_1,…,S_R`. That is a linear
program, and LP duality hands you the certificate for free.

## The certificate is a polynomial that hugs an indicator

The dual of "max `p_0` given moments `S_r`" is: find the cheapest polynomial
`g(t) = Σ_r y_r C(t,r)` that dominates the indicator `1[t=0]` on the lattice `{0,…,6}`.
Then for every E,

    meas(S7(E)) = Σ_t 1[t=0] p_t ≤ Σ_t g(t) p_t = Σ_r y_r S_r(E) =: L_y(E).

The beautiful part: the optimal `g` for each dangerous row factors with **integer roots
sitting inside {1,…,6}**:

    k=11,12,13:  g(t) = (t−3)(t−4)/12
    k=9,10:      g(t) = −(t−2)(t−3)(t−6)/36
    k=8:         g(t) = (t−1)(t−2)(t−4)(t−5)/40

Each is `≥ 0` at every integer `0≤t≤6` (zero only at the roots) and `g(0)=1`. Dual
feasibility — the entire rigor of the per-E bound — collapses to "this polynomial is
nonnegative at seven integers," a fact a human checks in one line. The `g` is the
**minimal-degree Bonferroni majorant** of the indicator. The classical even-truncation
Bonferroni bound `1 − S_1 + S_2 − … ` is the *crude* member of this family (its `g` is the
truncated alternating sum, which over-hangs the indicator); the LP picks the cheapest one,
and the cheapest one is what closes the dangerous rows where the crude one fails.

## Why this is the natural coordinate

THM-532 already learned that the seven sectors are "the natural coordinate" of the
problem. The moment view says more: the right *function* of those coordinates is the
**occupancy count** `N`, and the right *data* is its low factorial moments. Everything
about LRC's residual then lives in a 7-dimensional moment body, and `cap_k` is a single
hyperplane the body must stay below. The dual `g` is the supporting functional. This is
the same baby-Hodge / moment-SOS machinery (THM-509) the project built for the
cycle-count realizable region — the LRC sector cover turns out to be another truncated
moment problem, and it yields to the same dual-certificate hammer.

## Where it stops, honestly

The LP gives a **proved** per-E inequality `meas(S7(E)) ≤ L_y(E)` at L=7 sectors — exactly
the step the C·W bound could not take. What remains is the SCALAR extremal statement
`L_y(E) ≤ cap_k for all E`. Computation shows the clean form — "consec maximizes L_y" —
holds for the dangerous k=8,9,10 (exhaustive bounded-spread, zero violations) but FAILS
at k=11,12,13 with harmless AP-beaters far below cap. The mechanism is a sign-coupling:
consec MINIMIZES the single-miss moment S_1 (most uniform phases) and MAXIMIZES the
higher even moments S_2,S_4; the dangerous-row duals put a negative weight on S_1 and
positive on S_2,S_4, so both pulls point at consec. That coupling is the thing a final
proof must capture — a rearrangement inequality for the alternating moment combination,
not a separate bound per moment. The cap slack (huge at k=12,13, tiny 0.0014 at k=9)
absorbs the rest. LRC(14) is not proved; the analytic surface has been cut from a
two-sided product bound down to one scalar moment-rearrangement inequality.

## The transfer

Any "cover by L equal arcs, want full coverage" measure is `p_0` of an occupancy count,
hence a truncated moment problem with an integer-root Bonferroni dual. The finer-cover
route (THM-533, L=14 arcs) and this moment-LP route (L=7, cheapest dual) are the two
faces of the same coin: refine the partition, or refine the functional. The functional
refinement is cheaper — it closes the dangerous rows without enlarging the cover.
