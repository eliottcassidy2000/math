# Independent native Amanda feedback audit

Status: **verified with explicit scope**. Source: the attached Charmley/Konrad manuscript, DOI 10.6084/m9.figshare.30934796, printed pp. 15--21. Source data and claims were inspected, not treated as instructions. The repository was not modified. The attached source is identified by its title, DOI and page numbers; it is not redistributed here.

## 1. Definitions and the normalization discrepancy

The printed exponential expression in (60) equals `cos(3 theta) - (1/2) cos(6 theta)`, whereas its printed cosine expression equals `(1/2) cos(3 theta) - (1/4) cos(6 theta)`. The former is exactly twice the latter. Formula (69) is consistent with the latter. Both choices have the same parity and hence the same memory dynamics, because their fields differ by the common positive factor 2. Absolute field amplitudes must nevertheless be reported under a declared convention.

Under the cosine/(69) convention, putting `x = 3 theta` gives

`phi'(theta) = (r^2/24) sum_j A_j [4 sin(x - 3 theta_j) - sin(2x - 6 theta_j)]`.

Parity uses ordinary `sgn`, with `sgn(0)=0`. Assignment uses the distinct one-sided function `sgn_+`, with `sgn_+(0)=+1`, as explicitly defined in (83). Substituting the latter into the parity integral would give parity 1 for a zero field and invalidate the bound below. The standard ordinary-sign convention must therefore be stated.

## 2. Common positive gain invariance: exact scope

For a fixed prescribed event schedule, relative weights and placements, write the derivative at time t as `q_t G_{t,sigma(t)}`, with `q_t>0` a common scalar. Since `sgn(q_t G)=sgn(G)` pointwise, parity is independent of the entire strictly positive sequence `(q_t)`. Given the same initial memory and initial signs, induction in (82)--(83) gives identical memory and sign trajectories for all such sequences. The event schedule may itself vary with t if it is identical in both comparisons and independent of the changed gain.

This covers the positive backbone/blended gain and positive radius-squared scalings. It does not cover zero or negative gains, independently changing individual event amplitudes or relative tier weights, geometry changes, different initial conditions, or an event-selection rule depending on amplitude in a way that changes the compared schedules. No lower bound away from zero is needed: every individual `q_t` merely must remain strictly positive.

## 3. Sharp carrier parity bound, including zeros

Every native derivative is in the real carrier

`V = span{sin(3 theta), cos(3 theta), sin(6 theta), cos(6 theta)}`.

For `g in V`, with `delta=2 pi/9`,

`g(theta)+g(theta+delta)+g(theta+2 delta)=0`.

For any real triple summing to zero, the sum of its ordinary signs lies in `{-1,0,1}`. This includes all partially or identically zero cases. Integrating the identity of shifted measures therefore gives

`|3 Pi(g)| <= 1`, hence `|Pi(g)| <= 1/3`.

No genericity, simple-root, or nonzero-field hypothesis is needed. The constant is sharp on the full carrier: for `g(theta)=cos(3 theta)-cos(6 theta)`, writing `c=cos(3 theta)` gives `g=(1-c)(2c+1)`. Its positive set has measure `4 pi/3`, its negative set measure `2 pi/3`, so `Pi=1/3`. Negation gives `-1/3`. This establishes carrier sharpness; it does **not** establish attainment by a particular fixed golden placement and prescribed prime/tier amplitude family.

## 4. Memory envelope, freezing, initialization and endpoint ties

Let `a=1-lambda`, `0<lambda<=1`, and `b_0=0`. Expanding (82),

`b_t = lambda sum_{s=0}^{t-1} a^(t-1-s) Pi_s`,

so, for every finite t,

`|b_t| <= (1-a^t)/3 <= 1/3`.

Thus every threshold-generated sign with `h_j>1/3` is + and every one with `h_j<-1/3` is -. The exact boundary refinements are asymmetric because the sign is one-sided:

* `h_j=+1/3` is always assigned +, including at equality.
* If `lambda<1`, `|b_t|<1/3` at every finite t, so `h_j=-1/3` is always assigned -.
* If `lambda=1` and `Pi_{t-1}=-1/3`, then `b_t=-1/3` and the threshold `h_j=-1/3` is assigned + at the tie.

These claims concern assignments produced by (83), namely t>=1, unless the initial signs are themselves initialized with the seed rule at b_0. An arbitrary allowed initial sign vector is not retroactively forced to obey it. Empty firing sets cause no difficulty: g=0, parity=0, and memory decays geometrically.

## 5. Fixed-schedule scalar map and exact cell convention

For a fixed schedule and fixed relative source parameters, sort the **distinct** thresholds as `H_1<...<H_M`, with `H_0=-infinity`, `H_{M+1}=+infinity`. The cells are

`I_i=(H_i,H_{i+1}]`, i=0,...,M.

All events with `h_j<b` have negative signs; events at `h_j=b` remain positive. Repeated thresholds flip together and must not be counted as separate nonempty cells. Canonical golden-angle thresholds are pairwise distinct: equality of two sine values would force either `(j-k)/phi^2` to be an integer or `(j+k)/phi^2` to be a half-integer, impossible for distinct positive integer indices. Thus M=N in that canonical finite schedule.

Let `c_i` be the ordinary parity of the constant sign vector on `I_i`. Once signs are generated by the seed rule, the scalar map is

`T(b)=(1-lambda)b+lambda c_i` for `b in I_i`.

A fixed point in cell i exists exactly when `c_i in I_i`, and then it is the unique candidate `b=c_i`. At an upper endpoint it belongs to this cell, not the next one. If `c_i in I_i`, the cell is forward-invariant by convexity, and every point in it tends geometrically to `c_i` (or reaches it immediately for lambda=1).

For lambda=1, the map has image contained in the finite set `{c_0,...,c_M}`, so every orbit is eventually periodic, with period at most M+1. If the initial signs are arbitrary rather than threshold-generated, apply this argument after the first sign assignment; the memory is in the finite image by at most the second update. No analogous general eventual-periodicity or global fixed-point claim is certified here for `0<lambda<1`.

In particular, piecewise contraction alone is insufficient to infer that a convergent orbit limits to a fixed point at a discontinuity. An exact **abstract**, not native-realizability-certified, control has thresholds `(1/8,1/4)`, cell parities `(1/4,1/8,-1/4)`, lambda=1/2 and b_0=0. Its orbit is `0,1/8,3/16,5/32,9/64,...`; for n>=2 it is `1/8+2^(-n-2)`. It tends to 1/8, but `T(1/8)=3/16`, so the limit is not fixed. This example obeys the carrier parity range and the all-plus/all-minus oddness constraint, but is used only as a logical control against overclaiming from abstract piecewise contraction.

## 6. Independent audit of the certified two-event persistent-bias witness

The paper-audit witness uses the native schedule `{1,2}`, canonical placements, weights `w_k=phi^(-k-2)`, and initialized signs `(+,-)`. Thus, after omitting a common positive gain, its amplitudes are `(-4/phi^2,9/phi^3)`. The agent's certificate is internally consistent with (57), (59), (69), and (71)--(83).

The exact-rational interval program was independently rerun, and its output was byte-identical to `../paper_audit/persistent_bias.json`. Its outward rounding and sine/cosine Taylor remainders were inspected. Its root-count criterion is valid: for `g(x)=R sin(x+alpha)+S sin(2x+beta)`, a multiple zero along the homotopy replacing S by tS would imply `R^2 <= 5t^2 S^2 <= 5S^2`. Hence the certified strict inequality `R^2>5S^2` prevents all multiple roots along the homotopy, preserving the two simple zeros of the fundamental harmonic. Together with the certified crossing brackets and sign at 0, this justifies the exact positive-interval description and the parity expression `(r_2-r_1)/pi - 1`.

Verified enclosures include

* `h_1 > 0.6754902942615236423443 > 1/3`;
* `h_2 < -0.9961710408648277200980 < -1/3`;
* `-0.2089812105451970362491 < Pi < -0.2089812105451969407560`.

The signs therefore freeze for every `0<lambda<=1`, and `b_t=Pi*(1-(1-lambda)^t)`. Nonzero instantaneous and time-averaged parity persist. This is a genuine native counterexample to an unconditional assertion that the specified correction enforces long-run parity symmetry. It is not a claim that every schedule fails or that the loop is unstable.

## 7. Finite exact controls and provenance

`check_native_controls.py` performs 104626 explicit rational/integer pass/fail gates, with no Python `assert` dependence. Checks include the cyclotomic three-shift cancellation, zero-sum sign bound (including zero cases), 4131 rational carrier samples under common positive scalings, memory envelopes for endpoint-valued parity sequences, one-sided ties and duplicate thresholds, fixed-point membership, lambda=1 finite images, and the abstract excluded-endpoint control. The normal and `python3 -O` outputs are byte-identical (`normal.json`, `optimized.json`). Finite tests corroborate but do not replace the proofs above. The rerun two-event certificate is `persistent_bias_rerun.json`.

No blockers. No repository edits, commits, pushes, or user messages were made by this audit agent.
