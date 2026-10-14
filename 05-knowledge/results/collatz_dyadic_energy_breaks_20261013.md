# Collatz discrepancy through parity breaks and the unresolved integer clock

**Session:** opus S28, 2026-10-13 (worktree `collatz-posets-zeta5-20260927`).
**Owner's directive:** "keep working for a long session on proving the dyadic
energy bound, using the Q–Pell recursion and golden-tree controls to build
source-preserving refinement rules. pull in recent repo work and reroute
creatively around obstructions. think THM-4630 (a Q covariance update with an
exact two-generation Pell matrix) and on the Collatz amplitude problem (a
conditional amplitude certificate) ... open items: the white-sampling
hypothesis for Q, the limit of `A_k√k`, and the asymptotic decay of `ε_M`.
consider how for Conway, `2C(n)-n` even counts the recursion tree's terminal
leaves labelled 1. Hofstadter Q introduces an additional signed defect. For
BSD: the additional Tunnell test certifies 724 more unchanged curves,
increasing coverage from 2886 to 3610 of 4056, with 446 remaining. and prior
work surrounding `φ^8 + 1 = 7φ^4`, `φ^10 = 1 + 11φ^5`."

**Status after source-clock audit, 2026-10-13:** PROVED scoped break/jump identities and fixed-event-count covariance algebra; FINITE-EXACT census, with full gap counts corrected by the independent native census; NUMERICAL short-lag patterns; OPEN deterministic dyadic energy and amplitude. The original Theorem E's constants and nonzero-pole implication are withdrawn. See [the source-clock repair](collatz_break_clock_repair_20261013d.md) for the exact positive weighted criterion, factor-four control, and unchanged-current/different-energy hostile. The binary break identities below require a binary-labelled prefix; the general ROOT-supported sequence may also have zero labels. Independent audit here is focused on these stated claims, not a global audit of THM-4581.

Script `04-computation/experiments/collatz_dyadic_energy_breaks_20261013.py`; its original `.out` is retained as historical pre-audit output. The corrected function bodies have finite independent controls, and a separate native census supplies the replacement full counts and source-clock moments through `2^27`; see the repair note. No full rerun of the corrected Numba producer is claimed.

**What the repo had, what is new.** The Codex amplitude note defines
`t(n) = (-1)^{τ(n)}`, `w_n = 3·1_{n≡2(3)} - 1`, `H(N) = Σ_{n≤N} w_n t(n)`, proves
the boundary-current identity (6) `H(N) = W(N)t(N) - Σ_{n<N} W(n)(t(n+1) - t(n))`
with `W = 0, -1, +1` by `n mod 3`, verifies `|H(N)| < 1.5√N` and the dyadic
energy `∫_X^{2X} H(⌊x⌋)² dx < 2X²` to `2^34`, and proves (Theorem B3) that the
dyadic energy hypothesis (9a) for all `X ≥ 2^34` makes the martingale series
holomorphic in `Re s > 1/2` and the Mercator pole nonzero with a certified
residue interval. THM-4630 is the model of a *conditional closure*: an exact
identity plus measured stationarity hypotheses give a fixed point. New here:
(i) (6) read as a current over the **breaks** `β(n) = [t(n+1) ≠ t(n)]`, with
the break density census; (ii) the exact first-step identities for breaks and
the exact octave transfer, which exhibit a two-sided dependence and do not by themselves give a closed two-coordinate covariance map; (iii) the twisted-residue formula: the whole autocovariance of the
current is `(-1)^l` times the residue correlation of break positions, so the
Collatz "white-sampling hypothesis" is a statement about gaps between breaks
modulo 3; (iv) the measured six-lag memory, the failure of the renewal and
Markov closures, and the candidate energy closure with measured twelve-lag factor `0.84` (its original theorem is repaired below); (v) `ε_M = -H(3M)/(6M)`, which makes the open "decay of
`ε_M`" a sampled discrepancy question, distinct from the mean-square energy bound and explains the range-dependent
slopes; (vi) the Haar-coalescence route to the gap law.

## 1. `H` is a signed count of breaks (PROVED + FINITE-EXACT)

**Proposition 1 (binary-labelled prefix).** Let `β(n) = 1` iff `t(n+1) ≠ t(n)` (a *break*: the stopping
times of `n` and `n+1` differ by an odd number). Then
`H(N) = W(N) t(N) - 2 Σ_{n<N, β(n)=1} W(n) t(n+1)`,
so `|H(N)| ≤ 1 + 2 B_{1,2}(N)` where `B_{1,2}` counts breaks at `n ≢ 0 (mod 3)`.
*Proof.* In (6), `t(n+1) - t(n) = 2 t(n+1) β(n)`. ∎ (Identity (6) itself is
re-verified at every `N < 2^27`: no violation.)

**Census (all `n` through `2^27`; `p(n) = τ(n) mod 2` by descent and lookup).**

| octave `k` | breaks `B_k` | `B_k/2^k` | even-`n` share | `J_k = Σ W t(n+1)` over breaks | `J_k^2/(B1+B2)` | `max|H|/√(2^{k+1})` | `∫H²/X²` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 4537 | 0.277 | 0.376 | -85 | 2.41 | 1.12 | 0.64 |
| 18 | 69291 | 0.264 | 0.375 | 103 | 0.23 | 0.70 | 0.27 |
| 22 | 1053961 | 0.251 | 0.375 | 822 | 0.96 | 0.69 | 0.37 |
| 24 | 4122503 | 0.246 | 0.375 | -1403 | 0.72 | 1.33 | 1.77 |
| 26 | 16152641 | 0.241 | 0.374 | -2519 | 0.59 | 1.19 | 1.08 |

In these finite octaves about a quarter of the integers are breaks. The measured density decreases, with a fitted log2 slope -0.017 per octave over k=14..26; this is not an asymptotic law. Every binary parity break precludes synchronous merging of the pair (n,n+1), but the converse need not hold. Transferring a Haar no-merge tail at a growing stopping-time horizon needs a separate uniform estimate. Proposition 1 alone gives only an O(N) bound, so the measured square-root scale of H depends on signed cancellation. The octave energy ratios fluctuate between roughly 0.1 and 1.8 in the reported large octaves; their mean over octaves 8..26 is 0.45. These are finite summaries.

**Proposition 2 (first-step identities for breaks).** For every `m ≥ 2`:
`β(2m) = [t(m) ≠ t(3m+2)]` and `β(2m-1) = [t(3m-1) ≠ t(m)]`.
*Proof.* `τ(2m) = 1 + τ(m)`, `τ(2m+1) = 1 + τ(3m+2)`, `τ(2m-1) = 1 + τ(3m-1)`. ∎
(Verified at `44,739,239` indices.) So a break at an even `n = 2m` is a parity
disagreement between `m` and `3m+2`, at an odd `n = 2m-1` between `3m-1` and `m`:
every break compares a number with one `≡ 2 (mod 3)` in a *higher* octave.
The even-`n` share `0.375` is stable.

**Proposition 3 (the octave transfer, with boundary terms).** With
`lo = 2^k`, `hi = 2^{k+1}`,
`H(2hi-1) - H(2lo-1) = W(2hi-1)t(2hi-1) - W(2lo-1)t(2lo-1) + Σ_{m=lo}^{hi-1} [W(2m)(t(3m+2) - t(m)) + W(2m-1)(t(m) - t(3m-1))]`.
*Proof.* (6) over the octave, split by the parity of `n`, Proposition 2. ∎
(Verified exactly at `k = 12, 16, 20`; an earlier version without the boundary
terms was off by one.) The increment of `H` over octave `k+1` is a signed sum
over octave `k` of the cross-parities `t(3m±·) - t(m)`. Since `3m+2` and `3m-1`
lie in octaves `k+1` and `k+2`, the recursion is **two-sided**: this formula does not express the breaks of octave `k+1` solely in terms of the break data of octave `k`. This is the structural
difference from Hofstadter's `Q` (earlier indices whenever the recurrence is defined; the previous-octave approximation is part of the covariance model), so this calculation does not provide THM-4630's two-coordinate Pell covariance closure. Exact source-preserving refinements do exist: see [the ternary-wavelet hierarchy](collatz_wavelet_energy_20261013d.md). A transverse gap/current model is one possible route, with its source-clock obligation retained.

## 2. The twisted-residue formula and the measured memory (PROVED + FINITE-EXACT)

Index the breaks `n_1 < n_2 < ...` (all classes) and put `c_i = n_i mod 3`,
`v_i = W(c_i) t(n_i + 1)` (so `v_i = 0` on class `0` and `H` is `-2 Σ v_i` up to
the boundary term).

**Proposition 4 (binary labels; uncentered second moments).** The parity flips at every break, so `t(n_{i+l}+1) =
(-1)^l t(n_i+1)` and
`γ(l) := E[v_i v_{i+l}] = (-1)^l E[W(c_i) W(c_{i+l})]`:
the uncentered second moment of the current is the `(-1)^l`-twisted correlation of
the residues mod 3 of break positions, nothing else. *Proof.* As stated. ∎
(The direct products and the formula agree to every printed digit at lags
`1..12`.)

**Measured (breaks in `[2^22, 2^27)`, corrected total `31,569,545`; the old `31,484,229` counted only short gaps).** `γ(0)` is approximately `2/3`;
`γ(1..12) = 0.1406, -0.0619, -0.0676, -0.0325, -0.0183, -0.0149, -0.0017, -0.0030, 0.0027, -0.0013, 0.0043, -0.0013`.
In the convention of lags among nonzero-class breaks the normalized values
`0.127, -0.100, -0.059, -0.040, -0.014, -0.004, 0, 0` show a similar short-lag pattern on octaves `16, 20, 24, 26`, with some differences up to about0.008. This is finite evidence, not stationarity or an infinite six-lag cutoff.
Using the measured `γ(0)=0.6666953547...`, the empirical twelve-lag factor is `F=(γ(0)+2Σ_{l≤12}γ(l))/γ(0)=0.8356757365...`; the infinite covariance sum is not supplied by this statistic.

**The historical gap histogram below conditions on gaps<64; its normalization as the full law was wrong.** The repaired full probabilities modulo3 are `(.1926924,.5357171,.2715906)`;85,315 gaps at least64 were omitted. The old short-gap values were: `P(g = 1..8) = 0.408, 0.180,
0.117, 0.064, 0.048, 0.039, 0.032, 0.018`; `P(g ≡ 0, 1, 2 mod 3) = (0.192, 0.536,
0.271)`. In this finite census break positions are approximately equidistributed mod 3, and the transition
matrix is approximately circulant with rows
`(0.193, 0.536, 0.272)`: the next break prefers the class one step up.

**Proposition (renewal model; PROVED for the model, REFUTED as a closure).**
If the gaps were i.i.d., the initial residue uniform and independent, and `|ŵ|<1`, with `(a, b, c) = P(g ≡ 0, 1, 2 mod 3)`, then with
`ŵ = a + bω + cω²` (`ω = e^{2πi/3}`), `γ(l) = (-1)^l (2/3) Re(ŵ^l)` and
`F = Re[(1 - ŵ)/(1 + ŵ)]`. *Proof.* `E_c[W(c)W(c+s)] = 2/3, -1/3, -1/3` for
`s = 0, 1, 2`, and `P(S_l ≡ 0) = (1 + 2Re ŵ^l)/3`; sum the geometric series. ∎
With the measured `(a,b,c)`: `ŵ ≈ -0.2109614 + 0.2287403i`, `γ(1) ≈ 0.141` (close finite agreement; exact equality from the gap marginal also needs uniform phase independent of the gap), but `γ(2..6) ≈ -0.005,
-0.016, -0.006, -0.001, 0` against the measured `-0.062, -0.068, -0.033,
-0.018, -0.015`, and `F_renewal ≈ 1.338225` against `0.836`. A first-order Markov
chain on the classes (measured `3×3` matrix) fails the same way. **The observed low-order gap statistics do not fit the iid renewal model. No infinite short-memory assertion follows from twelve measured lags.** The analogy with THM-4630 is a proposed separation between an exact deterministic relation and an additional sampling model. Here alternation `(-1)^l` is exact on binary labels, while an asymptotic law for the residue sequence and its memory remains unproved; no six-lag cutoff is asserted.

## 3. Repaired fixed-index model theorem and the remaining clock obligation

The original Theorem E is superseded by [the exact source-clock audit](collatz_break_clock_repair_20261013d.md#3-repair-of-the-incoming-stochastic-theorem). It omitted the factor4 in H=-2J+boundary, mixed all-break and nonzero-break normalization, and inferred a deterministic amplitude certificate from a model expectation.

For a specified mean-zero stationary process with summable covariance gamma, and **deterministic event count m**, the valid identity is

    E[J_m^2]=m gamma(0)+2 sum_(l=1)^(m-1) (m-l)gamma(l).

Hence `E[(-2J_m)^2]=4m sigma^2+o(m)`, where `sigma^2=gamma(0)+2 sum_(l>=1)gamma(l)`. With gamma(0)=2/3 and sigma^2=F gamma(0), the coefficient is `(8F/3)m`. Weighted covariance summability is needed for an O(1) remainder; a bounded boundary can contribute an O(sqrt(m)) cross-term.

The measured twelve-lag factor F≈0.836 does not prove summability, stationarity, mean zero, or that this formula remains valid after substituting the dependent clock m=B(N). Even with an additional deterministic clock B(N)≈dN and uniform moment control, the heuristic expected dyadic coefficient is `4Fd≈0.806`, not0.20. Expectations do not bound a designated integer source.

A uniform expected bound on a suitable probability space would imply an eventual almost-sure log-squared bound by Markov and Borel–Cantelli. It supplies neither the explicit cutoff and constant required by the existing amplitude certificate nor an actual law for the fixed Collatz sequence. **The original nonzero-pole conclusion is withdrawn.** The prior conditional certificate remains unchanged.

The deterministic replacement is precise: retain plateau lengths L_i and the current J_i inherited from source1. If `S_X=sum L_i J_i^2<=alpha X^2`, alpha<1/2, then `E_X<=2X^2` for `X>=(sqrt(2)-2sqrt(alpha))^(-2)`. The premise is OPEN; it retains the source clock rather than replacing it with event-index stationarity.

## 4. `ε_M` (PROVED identity; the open decay item re-read)

**Proposition 5.** With `A(x) = (1/x) Σ_{n≤x} t(n)`, Codex's quantity
`ε_M = A(2M) + A(M)/2 + A(3M)/2` satisfies `ε_M = -H(3M)/(6M)` exactly.
*Proof.* `t(2m) = -t(m)`, `t(2m-1) = -t(3m-1)`, and the images `3m-1 ≤ 3M-1`
are exactly the `j ≡ 2 (mod 3)` below `3M`; `Σ_{j≤3M, j≡2} t(j) = (3M A(3M) +
H(3M))/3`. ∎ (Verified at every dyadic `M ≤ 2^25`.)

The identity directly relates epsilon to sampled discrepancy; a pointwise rate and a dyadic mean-square bound remain different quantified statements. With the additional deterministic-clock model of section3, the corrected heuristic is `rms(epsilon_M)≈sqrt((8F/3)d(3M))/(6M)≈0.212 M^(-1/2)` at F=.84,d=.24; this is not a consequence of the measured short-lag data alone. Historically, measured `rms(ε_M)·√M` over 64 prefixes per octave runs `0.06–0.24`
(mean `0.13`, high in the excursion `2^22–2^26`), and the `log_2` slopes are
`-0.44` (`2^12..2^22`), `-0.41` (`2^18..2^25`), `-0.42` (`2^22..2^25`). The
range-dependent slopes of S27 (`-0.45`, `-0.7`, `-1.2`; MISTAKE-599's caution)
are the walk's own excursions seen on short ranges; the model's exponent is
`1/2`, conditional, and the data bracket it no better than `[0.4, 0.5]` on
`2^12..2^25`.

## 5. The obligation, and the route through Haar coalescence

A binary parity break at `n` precludes synchronous merging of `(n,n+1)`. If both sources reach 1, their first-hit time difference is odd; conversely that odd difference is precisely a parity break. Failure of synchronous merging alone is not an equivalent test for a break. THM-4581 (mac-mini) proves
synchronous merging almost surely for `u = 3^k v + e`, hence for `y, y+1`
2-adically, with no-merge probability of order `T^{-1/2}` by time `T`; THM-4627
transfers this to integers with density one. Truncated gap events and finite-lag products can be probed through joint statistics of a growing window of pairs `(n+g,n+g+1)`. A fixed window such as `g<=20` cannot determine the full gap law: the present census has gaps as large as 1412. The coupled-reader Haar model supplies finite-horizon approximations (THM-4565's kernel is the two-reader case), and the computation is the Collatz version of "source-
preserving refinement": alongside the exact two-sided refinement hierarchy,
refine the *transverse* state — the joint merge pattern of a window of
consecutive pairs — until the omitted correlations are charged to a summable
remainder. Reproducing a finite gap law and six lags in a Haar model would be a control, not a proof of H_C or of the deterministic energy bound. A transfer needs uniform growing-horizon tail control, the integer residence lengths, and an authenticated source measure. That is the repaired D57.

## 6. Against the one-sided controls, and the cross-citation board

* **Conway's recursion** (`a(n) = a(a(n-1)) + a(n-a(n-1))`): the Codex leaf
  charge `e(n) = 2a(n) - n = e(j) + e(n-j)` counts the leaves labelled 1 of
  the computation tree (golden/nested note (4)); it is *nonnegative*, so the
  one-sided control has a current with no cancellation at all, the opposite
  of `H`. **Hofstadter's Q** adds the signed defect `E_Q(n) = E_Q(p) + E_Q(q) +
  (p + q - n)` ((5) there): a local mass defect, the first place where a sign
  enters. Collatz's current is signed from the start (the twist `(-1)^l`) and
  two-sided.
* **The pinned Conway generation map to `k = 400`** (S27's open item on
  `A_k√k`): `A_k√k = 0.1451, 0.1446, 0.1439, 0.1437, 0.1440, 0.1445, 0.1443,
  0.1431, 0.1416, 0.1387, 0.1372` at `k = 20, 25, 29, 40, 60, 80, 100, 150, 200,
  300, 400`; the local exponent of `A_k ~ k^{-e}` is `0.50–0.52` to `k = 100`
  and `0.54–0.55` from `150` to `400`; neither `C/√(k + k_0)` (rel. rms `5%`)
  nor `C/√(k log(k/k_1))` (`2.7%`) fits the eleven points well. The pinned map does not decide the actual sequence's limit. [The concurrent classical-Pascal continuation](conway_golden_memory_20261013d.md) now closes the integer limit at0.138014902399..., with an O(1/k) scaled-maximum error and an explicit nonlinear eigenprofile. Only the lagged model's approximation error and its own long-time behaviour remain open.
* **White sampling for Q (H1):** untouched here; Proposition 4 shows what its
  proposed Collatz analogue could be (a residue statistic of binary breaks), which is the kind of
  object a finite computation can reach.
* **BSD Tunnell terminals** (`bsd_tunnell_terminals_20261013c.md`, table
  line "Total 1,170 / 724 / 446"): the additional Tunnell gate certifies 724
  more unchanged curves, coverage `2886 → 3610` of `4056`, `446` remaining —
  recorded here for cross-citation only.
* **The Lucas identities** `φ^8 + 1 = 7φ^4`, `φ^10 = 1 + 11φ^5` (eleventh S15
  note; the circulant note's Pillai–Lucas census): the Pell unit of Codex's
  (16), `3 + 2√2`, is the silver companion of the golden `φ`; the eleventh
  note's "unit minus one" clocks are the common frame (`|N(u^j - 1)|`
  torsion orders), but no Collatz quantity in this note is a Pell or Lucas
  number — the Pell structure belongs to Q's one-sided covariance map.

## 7. Verdicts

| claim | status |
|---|---|
| Propositions 1–5 with binary-label scope for the break shortcut and twist; generalized jump formula retains zero labels | PROVED scoped; finite identities checked |
| finite break counts and raw lag products; rounded correlations and truncated twelve-lag F≈0.836 | FINITE-EXACT counts, NUMERICAL summaries; full gap histogram corrected |
| renewal formula `γ(l) = (-1)^l (2/3) Re ŵ^l`, `F = Re[(1-ŵ)/(1+ŵ)]` | PROVED for the model; REFUTED as a closure (lags `≥ 2`) |
| Repaired model second moment; original Theorem E constants and amplitude conclusion withdrawn | PROVED at deterministic event count under stated model assumptions; source-clock transfer OPEN |
| `ε_M` slopes `-0.44/-0.41/-0.42`; model exponent `1/2` | NUMERICAL / CONDITIONAL |
| Conway pinned map: `A_k√k` drifting `0.145 → 0.137`, exponent `0.50 → 0.55` | NUMERICAL (open) |
| LRC(14), Collatz | OPEN |

## 8. Directions

* **D57 (the Haar gap law).** Compute in the Haar coalescence chain
  (THM-4581/4565) the joint law of synchronous-merge failures for the pairs
  `(n+g, n+g+1)`, `g ≤ 20`, hence the truncated gap statistics and lag products; retain the omitted long-gap and horizon tails, and compare finite approximations with
  `(0.192, 0.536, 0.271)` and `0.141, -0.062, -0.068, -0.033, -0.018, -0.015`.
* **D58 (the census to `2^34`).** Re-run the break census on the Codex packed
  census (`collatz_discrepancy_20261013c`) to test the finite scale stability of
  `γ(l)` and the density slope by seven octaves.
* **D59 (odd-offset merges).** Separate the breaks into "never merge" and
  "merge with odd offset" using THM-4565's offset classes; the density of
  the second kind may have its own rate.
* **D60 (the pinned map's error).** Compare against the now-proved classical integer limit and exact Pascal maxima through generation1024. A uniform error/stability estimate is needed for the lagged model, not another numerical attempt to decide the integer limit.
