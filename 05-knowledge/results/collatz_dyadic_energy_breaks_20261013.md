# The dyadic energy of the signed parity discrepancy, reduced to the break process: `H` is a signed count of the parity breaks between consecutive integers, the breaks of an octave are the cross-parities of the previous one (exact, two-sided, no octave-to-octave closure), the autocovariance of the break current is the `(-1)^l`-twisted residue correlation of break positions (exact), it is scale-invariant with six lags of memory (measured), the energy theorem is conditional on exactly that (the Collatz analogue of THM-4630's H1–H3), `ε_M = -H(3M)/(6M)` exactly, and the gap law of breaks is the obligation THM-4581's coalescence chain can address

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

**Status: PROVED (elementary) for Propositions 1–5 (the break current, the
first-step break identities, the octave transfer with its boundary terms, the
twisted-residue formula, the `ε_M` identity) and for the renewal-model
formula; FINITE-EXACT for the census to `2^27` (every identity verified at
every index; the gap law; the autocovariance to twelve lags on five octaves);
CONDITIONAL for Theorem E (the energy bound under the stationary
short-memory hypothesis H_C, measured and not proved); REFUTED: the
i.i.d.-gap renewal closure and the first-order Markov closure beyond lag 1;
NUMERICAL: the `ε_M` slopes and the pinned Conway map to `k = 400`; CITED:
the Codex amplitude note (identities (6), (7), the census to `2^34`,
hypothesis (9a), Theorem B3), THM-4630, THM-4581, THM-4627, THM-4565, the
Codex golden/nested note (leaf charge (4), Q defect (5), the Pell conjugacy
(14)–(17)), the Tunnell note. LRC and Collatz OPEN. Independent audit OWED.**
Script `04-computation/experiments/collatz_dyadic_energy_breaks_20261013.py`
(output `.out`, 6 s with numba to `2^27`).

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
the exact octave transfer, which show the process is two-sided (no
octave-to-octave covariance map, hence no Pell structure) and locate the
"source"; (iii) the twisted-residue formula: the whole autocovariance of the
current is `(-1)^l` times the residue correlation of break positions, so the
Collatz "white-sampling hypothesis" is a statement about gaps between breaks
modulo 3; (iv) the measured six-lag memory, the failure of the renewal and
Markov closures, and the conditional energy theorem with the measured
factor `0.84`; (v) `ε_M = -H(3M)/(6M)`, which makes the open "decay of
`ε_M`" the same question as the energy and explains the range-dependent
slopes; (vi) the Haar-coalescence route to the gap law.

## 1. `H` is a signed count of breaks (PROVED + FINITE-EXACT)

**Proposition 1.** Let `β(n) = 1` iff `t(n+1) ≠ t(n)` (a *break*: the stopping
times of `n` and `n+1` differ by an odd number). Then
`H(N) = W(N) t(N) - 2 Σ_{n<N, β(n)=1} W(n) t(n+1)`,
so `|H(N)| ≤ 1 + 2 B_{1,2}(N)` where `B_{1,2}` counts breaks at `n ≢ 0 (mod 3)`.
*Proof.* In (6), `t(n+1) - t(n) = 2 t(n+1) β(n)`. ∎ (Identity (6) itself is
re-verified at every `N < 2^27`: no violation.)

**Census (odd `n` through `2^27`; `p(n) = τ(n) mod 2` by descent and lookup).**

| octave `k` | breaks `B_k` | `B_k/2^k` | even-`n` share | `J_k = Σ W t(n+1)` over breaks | `J_k^2/(B1+B2)` | `max|H|/√(2^{k+1})` | `∫H²/X²` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 4537 | 0.277 | 0.376 | -85 | 2.41 | 1.12 | 0.64 |
| 18 | 69291 | 0.264 | 0.375 | 103 | 0.23 | 0.70 | 0.27 |
| 22 | 1053961 | 0.251 | 0.375 | 822 | 0.96 | 0.69 | 0.37 |
| 24 | 4122503 | 0.246 | 0.375 | -1403 | 0.72 | 1.33 | 1.77 |
| 26 | 16152641 | 0.241 | 0.374 | -2519 | 0.59 | 1.19 | 1.08 |

A quarter of all integers are breaks, and the density falls only slowly
(`log_2` slope `-0.017` per octave over `k = 14..26`, i.e. about `k^{-0.22}`):
breaks are the pairs `(n, n+1)` whose orbits do not merge synchronously in
the sense of THM-4581 (whose no-merge-by-time-`T` probability is of order
`T^{-1/2}` up to logarithms, `T` here being the stopping time, of order
`log n`), so the trivial bound of Proposition 1 is `O(N)` and useless; the
`√N` size of `H` is entirely signed cancellation among `0.24 N` terms. The
octave energies `∫_X^{2X} H²/X²` (`0.1`–`1.8`, Codex's `1.77` at `2^24` reproduced)
fluctuate as a random walk's do; their mean over octaves `8..26` is `0.45`.

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
lie in octaves `k+1` and `k+2`, the recursion is **two-sided**: the breaks of
octave `k+1` are not a function of octave `k`'s breaks. This is the structural
difference from Hofstadter's `Q` (one-sided: both backward indices are in the
previous generation), and it is why THM-4630's octave-to-octave covariance
map `r ↦ (r-1)/(3-r)` with its two-generation Pell matrix (Codex golden/nested
§6, (14)–(16)) has **no Collatz counterpart**: there is no generation
operator to iterate. The closure, if any, must be *transverse* (stationary
along `n` at a fixed scale), which is what §2–3 test.

## 2. The twisted-residue formula and the measured memory (PROVED + FINITE-EXACT)

Index the breaks `n_1 < n_2 < ...` (all classes) and put `c_i = n_i mod 3`,
`v_i = W(c_i) t(n_i + 1)` (so `v_i = 0` on class `0` and `H` is `-2 Σ v_i` up to
the boundary term).

**Proposition 4.** The parity flips at every break, so `t(n_{i+l}+1) =
(-1)^l t(n_i+1)` and
`γ(l) := E[v_i v_{i+l}] = (-1)^l E[W(c_i) W(c_{i+l})]`:
the autocovariance of the current is the `(-1)^l`-twisted correlation of
the residues mod 3 of break positions, nothing else. *Proof.* As stated. ∎
(The direct products and the formula agree to every printed digit at lags
`1..12`.)

**Measured (breaks in `[2^22, 2^27)`, `31,484,229` of them).** `γ(0) = 2/3`;
`γ(1..12) = 0.1406, -0.0619, -0.0676, -0.0325, -0.0183, -0.0149, -0.0017, -0.0030, 0.0027, -0.0013, 0.0043, -0.0013`.
In the convention of lags among nonzero-class breaks the normalized values
`0.127, -0.100, -0.059, -0.040, -0.014, -0.004, 0, 0` are identical to three
decimals on octaves `16, 20, 24, 26`: a **scale-invariant six-lag memory**.
The empirical long-run variance factor `F = (γ(0) + 2Σ_{l≤12} γ(l))/γ(0) = 0.836`.

**The gap law.** Gaps between consecutive breaks: `P(g = 1..8) = 0.408, 0.180,
0.117, 0.064, 0.048, 0.039, 0.032, 0.018`; `P(g ≡ 0, 1, 2 mod 3) = (0.192, 0.536,
0.271)`. Break positions are equidistributed mod 3, and the transition
matrix of consecutive break classes is circulant with rows
`(0.193, 0.536, 0.272)`: the next break prefers the class one step up.

**Proposition (renewal model; PROVED for the model, REFUTED as a closure).**
If the gaps were i.i.d. with `(a, b, c) = P(g ≡ 0, 1, 2 mod 3)`, then with
`ŵ = a + bω + cω²` (`ω = e^{2πi/3}`), `γ(l) = (-1)^l (2/3) Re(ŵ^l)` and
`F = Re[(1 - ŵ)/(1 + ŵ)]`. *Proof.* `E_c[W(c)W(c+s)] = 2/3, -1/3, -1/3` for
`s = 0, 1, 2`, and `P(S_l ≡ 0) = (1 + 2Re ŵ^l)/3`; sum the geometric series. ∎
With the measured `(a,b,c)`: `ŵ = -0.2115 + 0.2293i`, `γ(1) = 0.141` (exact
agreement, as lag 1 depends only on the gap law), but `γ(2..6) ≈ -0.005,
-0.016, -0.006, -0.001, 0` against the measured `-0.062, -0.068, -0.033,
-0.018, -0.015`, and `F_renewal = 1.339` against `0.836`. A first-order Markov
chain on the classes (measured `3×3` matrix) fails the same way. **The gaps
carry memory over about six breaks (clustering of non-merging pairs); the
current is short-memory but not renewal.** This is the Collatz form of the
"source is white except for the shared index" finding of THM-4630: here the
exact structure is the alternation `(-1)^l` (deterministic), the residual
randomness is the residue of the next break, and the memory is six lags
instead of one.

## 3. The conditional energy theorem (CONDITIONAL; the Collatz analogue of THM-4630)

**Hypothesis H_C (measured, not proved).** Within every octave the current
`(v_i)` is stationary with absolutely summable autocovariance, with
`γ(0) = 2/3` and `F = (γ(0) + 2Σ_{l≥1} γ(l))/γ(0) = 0.84 ± 0.02` (the values
above, stable from `2^16` to `2^27`), and the break density `d(X) = B(X)/X` is
non-increasing.

**Theorem E.** Under H_C, `E[H(N)²] = (4F/3)·B_{1,2}(N) + O(1) = (2F/3) B(N) + O(1)
≤ (2F/3) d(N) N`, so for every dyadic `X`
`E[∫_X^{2X} H(⌊x⌋)² dx] ≤ F d(X) X² + O(X) ≈ 0.20 X²`,
ten times inside Codex's hypothesis (9a) (`≤ 2X²`); and, under H_C with
Gaussian-type fourth moments, `∫_X^{2X} H² ≤ X² (log_2 X)^2` for all but
finitely many dyadic `X` almost surely (Chebyshev and `Σ_k k^{-2} < ∞`), which
is reroute 1 of the Codex note with `A = 2`, hence (Theorem B3 there) the
Mercator pole is nonzero with the certified residue interval. *Proof.* `H`
is `-2` times a sum of `B_{1,2}(N)` stationary short-memory terms plus a
bounded boundary term; the variance of such a sum is `F γ(0)` per nonzero
term; the octave integral is bounded by its right endpoint. ∎

**What is and is not proved.** Everything deterministic in §1–2 is proved;
H_C is exactly as unproved as THM-4630's H1–H3, and it is *measured* with the
same precision (three decimals over ten octaves). It is a cleaner hypothesis
than Q's: by Proposition 4 it concerns only the residues mod 3 of break
positions, i.e. the joint law of the gaps `n_{i+1} - n_i` — a local
arithmetic statistic of coalescence, not a global one.

## 4. `ε_M` (PROVED identity; the open decay item re-read)

**Proposition 5.** With `A(x) = (1/x) Σ_{n≤x} t(n)`, Codex's quantity
`ε_M = A(2M) + A(M)/2 + A(3M)/2` satisfies `ε_M = -H(3M)/(6M)` exactly.
*Proof.* `t(2m) = -t(m)`, `t(2m-1) = -t(3m-1)`, and the images `3m-1 ≤ 3M-1`
are exactly the `j ≡ 2 (mod 3)` below `3M`; `Σ_{j≤3M, j≡2} t(j) = (3M A(3M) +
H(3M))/3`. ∎ (Verified at every dyadic `M ≤ 2^25`.)

So the open "asymptotic decay of `ε_M`" is the energy question itself. Under
H_C, `rms(ε_M) ≈ √(F·(2/3)·d·3M)/(6M) ≈ 0.105 M^{-1/2}` with a slowly falling
prefactor; measured `rms(ε_M)·√M` over 64 prefixes per octave runs `0.06–0.24`
(mean `0.13`, high in the excursion `2^22–2^26`), and the `log_2` slopes are
`-0.44` (`2^12..2^22`), `-0.41` (`2^18..2^25`), `-0.42` (`2^22..2^25`). The
range-dependent slopes of S27 (`-0.45`, `-0.7`, `-1.2`; MISTAKE-599's caution)
are the walk's own excursions seen on short ranges; the model's exponent is
`1/2`, conditional, and the data bracket it no better than `[0.4, 0.5]` on
`2^12..2^25`.

## 5. The obligation, and the route through Haar coalescence

A break at `n` is a pair `(n, n+1)` whose orbits either never merge
synchronously or merge with an odd time offset. THM-4581 (mac-mini) proves
synchronous merging almost surely for `u = 3^k v + e`, hence for `y, y+1`
2-adically, with no-merge probability of order `T^{-1/2}` by time `T`; THM-4627
transfers this to integers with density one. The *gap law* and the six-lag
memory of §2 are joint statistics of the merge events of the pairs
`(n+g, n+g+1)` for `g ≤ 20` — a finite collection of coupled readers of
one tape (THM-4565's kernel is the two-reader case). In the Haar model these
are computable, and the computation is the Collatz version of "source-
preserving refinement": instead of refining octave by octave (impossible, §1),
refine the *transverse* state — the joint merge pattern of a window of
consecutive pairs — until the omitted correlations are charged to a summable
remainder. If the Haar model reproduces `(0.192, 0.536, 0.271)` and the six
lags, H_C reduces to THM-4581/4627's framework and the energy bound becomes a
coalescence computation. That is D57.

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
  nor `C/√(k log(k/k_1))` (`2.7%`) fits the eleven points well. The limit of
  `A_k√k` is not settled by the pinned map; a logarithmic correction is not
  excluded (NUMERICAL; the discretization of the pinned map at `k ≥ 150` is
  itself a suspect).
* **White sampling for Q (H1):** untouched here; Proposition 4 shows what its
  Collatz analogue *is* (a residue statistic of breaks), which is the kind of
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
| Propositions 1–5 (break current; first-step break identities; octave transfer with boundary; twisted-residue formula; `ε_M = -H(3M)/(6M)`) | PROVED, verified at every index to `2^27` |
| break density `0.28 → 0.24`, slope `-0.017`/octave; even share `0.375`; gap law; circulant transitions; `γ(1..12)`; `F = 0.836` | FINITE-EXACT |
| renewal formula `γ(l) = (-1)^l (2/3) Re ŵ^l`, `F = Re[(1-ŵ)/(1+ŵ)]` | PROVED for the model; REFUTED as a closure (lags `≥ 2`) |
| Theorem E: energy `≈ 0.20 X²` in expectation, `(9a)` with `A = 2` a.s. | CONDITIONAL on H_C |
| `ε_M` slopes `-0.44/-0.41/-0.42`; model exponent `1/2` | NUMERICAL / CONDITIONAL |
| Conway pinned map: `A_k√k` drifting `0.145 → 0.137`, exponent `0.50 → 0.55` | NUMERICAL (open) |
| LRC(14), Collatz | OPEN |

## 8. Directions

* **D57 (the Haar gap law).** Compute in the Haar coalescence chain
  (THM-4581/4565) the joint law of synchronous-merge failures for the pairs
  `(n+g, n+g+1)`, `g ≤ 20`, hence the gap law mod 3 and `γ(1..6)`; compare with
  `(0.192, 0.536, 0.271)` and `0.141, -0.062, -0.068, -0.033, -0.018, -0.015`.
* **D58 (the census to `2^34`).** Re-run the break census on the Codex packed
  census (`collatz_discrepancy_20261013c`) to extend the scale invariance of
  `γ(l)` and the density slope by seven octaves.
* **D59 (odd-offset merges).** Separate the breaks into "never merge" and
  "merge with odd offset" using THM-4565's offset classes; the density of
  the second kind may have its own rate.
* **D60 (the pinned map's boundary layer).** The exact `k`-term initial run
  and the `O(1/n)` selector corrections (S27 next step 1) before deciding the
  limit of `A_k√k`; the drift `0.145 → 0.137` must be shown to be the map's
  own, not the grid's.
