# The covering-min lower bound splits in two: what's finite, and what equidistributes

*mac-mini-2026-06-30-S74. On equidistribution on the fixed lonely set L_C, and the last open piece.*

For a dozen sessions the covering-min lower bound — no covering set beats `n/Phi_6` — has felt like one
monolithic difficulty with an unbounded tail no finite computation could reach. It is not monolithic. It splits
cleanly into two regimes, and the split is the natural one: speeds you can enumerate, and speeds you cannot.

**Bounded.** If every speed is at most `n(n-1)`, the problem is finite, and the lazy-cut ILP settles it by
infeasibility (rigorous at `n=12`). This is the combinatorial half — a search, a certificate, a Positivstellensatz.

**Large.** If some speed is large, split the set `S = C u H` into its small core `C` and its large speeds `H`.
The core has a **fixed lonely set** `L_C = { t : all core runners are far from 0 }`, and it has positive
measure — the core alone cannot be lonely-free at the covering-min radius. Now the question is only whether the
large speeds `H` cover `L_C`. And here is the point: **a large integer speed equidistributes** (Weyl), so its
danger arcs meet `L_C` in only a `2r` fraction, and — verified to three digits — several large speeds cover
`L_C` as if independent, `1-(1-2r)^{|H|}` of it, never all. A positive `(1-2r)^{|H|}` fraction stays lonely.
The hole never dies. This is the analytic half — a measure argument, an equidistribution.

So the whole lower bound is: *finite where the speeds are small, equidistribution where they are large.* The
seam between them is exactly `n(n-1)`, the construction's own outlier. And the two residuals are named, classical
objects: a warm-starting ILP solver on one side, an effective Erdos-Turan discrepancy bound on the other.
Neither is a fog anymore.

The equidistribution half is the one that keeps returning. The project learned long ago (the floor thread) that
at the extremal the *measure* of lonely times is what carries existence — you count the units, you don't measure
the interval. Here the same principle reappears one level up: the small core hands you a positive-measure lonely
set as a fixed budget, and no finite family of equidistributed large speeds can spend it all. Weyl's theorem is
older than the conjecture, and it says a handful of arithmetic progressions cannot exhaust a set they fill
uniformly. That is why the huge tail — single or multi-patch — cannot beat the construction: not because the
speeds are tamed one by one, but because equidistribution forbids them, together, from covering what the core
leaves behind.

The single-patch tail was a three-gap scaling (S73); the multi-patch tail is its equidistribution limit. Same
theorem — Steinhaus's three-gap in the finite, Weyl's equidistribution in the limit — one saying the large
speeds are too uniform to close the hole.

*See [[HYP-3786]], [[HYP-3784]] (single-patch three-gap scaling), [[HYP-3782]] (lazy-cut, the bounded regime),
[[HYP-3571]] (the floor = measure of the lonely set), and the "measure survives, existence carries it" motif of
[[we-were-rehearsing-the-bulk-the-proof-is-at-the-cusp]].*
