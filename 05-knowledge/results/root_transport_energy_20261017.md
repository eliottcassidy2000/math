# Uniform depth control from the other end: the signed ROOT transport as an energy process

2026-10-17, opus-2026-10-17-S41. Owner prompt: uniform depth and scale bounds for ROOT transport; close the remaining
Collatz gap (uniform control of cancellation and error as depth grows) by first imagining a system where that control is
inherent, then merging. Sources, ports and integer clocks are preserved throughout: every quantity below is computed on
the literal seeds and ROOT port. Canonical primitive depth is distinct from elapsed shortcut time;
the terminal braid correction described below retains the latter.

**Scope repair, 2026-10-17, clock_duality_oct15.** THM-4667 now defines f as zero off the ROOT basin and retains
that basin's energy complement. Its coefficient ledger and finite censuses survive; the uniform current-budget
implication needs no assumption that every source converges. Primitive Fourier phases require the boundary kernels
from [positive ROOT transport, section2](root_positive_transport_20261017g.md). See
[outer packets, section7](outer_packet_golden_clocks_20261017h.md) and the dated MISTAKES entry.

**Typing.**
- PROVED (elementary): THM-4667, which contains:
  - the energy ledger;
  - the sibling factorisation of every fibre;
  - the ROOT/shell identity for the termination regime;
  - the Fourier form of the ideal transport.
- FINITE-EXACT:
  - the ideal energy ladder to `d = 6` in exact integers;
  - the actual energy profiles of **every** `N < 2^23` (an exact incremental algorithm; invariant checked against the
    direct current at every `N`).
- VERIFIED (float64):
  - the ideal ladder for `d = 7..12`, with a rigorous band-tail bound (rounding not enclosed);
  - the cross-term anatomies, the operator norms and the matching tests.
- Independent audit: 04-computation/experiments/s41_audit_20261017/REPORT.md. Its fixes are applied below.
- OPEN: HYP-9288, revised. No claim of closing the gap.

**Working decomposition.** The finite energy observations suggest two regimes.
- **Computed ideal regime.** The tested depths dissipate. The observed mechanism is anti-correlated partial currents through
  sibling seeds `s` and `2s + 1`. It is a property of the actual vectors, not of the operator, which is not a
  contraction.
- **Sampled termination regime.** In the sampled large-current rows, the energy sits almost entirely at ROOT from depth about
  40 on. The ROOT coefficient is a partial sum of the current over stopping-time shells, which is the census object of
  HYP-9283/9285.

A uniform shell bound is therefore a promising part of the target. It is not equivalent to the full energy bound:
uniform ideal cross terms, discrete-to-density matching and the nonterminal energy also need control. The exact data
to `2^23` support the full bound in the tested octave averages.

## 1. Where clock_duality left the gap

In [root_braid_logistic_20261017f.md](root_braid_logistic_20261017f.md) section 3c the current `R(N) = H(4N) - H(N)` is
pushed through `d` exact odd transports. The l1 mass of the transported source vector bounds `|R(N)|` for every
compatible profile:

    A^(d)(N) = γ_d N + O_d(log^(d+1) N),   γ_0, ..., γ_3 = 2.2857, 1.7841, 1.4489, 1.2281.

Their next critical target is a coupled depth/source-scale estimate: `γ_d` decaying at a rate that survives
`d ≍ log N`, with a remainder uniform in `d`. clock_duality's certified table (root_positive_transport_20261017g.md
section 6) continues with `γ_4, γ_5, γ_6 = 1.0798, 0.9671, 0.8773` and up to `γ_10`. This session reproduces those values and
adds `γ_11`, `γ_12` (section 4).
- The successive ratios are 0.78, 0.81, 0.85, 0.88, 0.90, 0.91. They creep towards 1.
- l1 is the worst case over all sign profiles. It throws away exactly the cancellation the signs provide.

## 2. The other end: systems whose cancellation is controlled at every depth

**Rudin-Shapiro.** The archetype is the Rudin-Shapiro pair. On |z|=1, a butterfly `(P, Q) -> (P + zQ, P - zQ)` conserves
`|P|^2 + |Q|^2` after dividing both outputs by sqrt2 (the unnormalized transform doubles the energy). The classical normalized construction motivates square-root partial-sum control.

**Martingale.** The same control holds for a random-merge process. If every merge combined currents with independent
fair relative signs, the energy `Σ c^2` would be a nonnegative martingale. Doob gives
`Pr(sup_d E_d >= lambda E_0) <= 1/lambda`. This controls its maximal tail in that model; it is not a deterministic
Collatz inequality.

**The Collatz transport is half a butterfly.** clock_duality's R5 pairs the seeds `s ≡ 1 (mod 8)` with `2s + 1`:
- `U(2s + 1) = B(U(s))`, so the siblings merge after one odd transport;
- they have opposite clock parity: `a(s) = 2` and `a(2s + 1) = 1`.

THM-4667 statement 8 completes this. **Every fibre is a disjoint union of sibling pairs**, at scales spaced by 64, plus
one singleton `(2u - 1)/3` when `u ≡ 2 (mod 3)`. So

    (Tc)(u) = -[u ≡ 2 (3)] c((2u - 1)/3) + Σ_k (c(s_k) - c(2 s_k + 1)).

A merge keeps only the difference channel `c - c'` of each butterfly and drops the sum channel `c + c'`. The energy change
of one pair is `((c - c')^2 - (c + c')^2)/2 = -2 c c'`. So the energy falls exactly when the siblings carry correlated
coefficients. On rooted siblings these are anti-correlated partial currents, since `f(2s+1)=-f(s)`; off ROOT both f-values are zero and cannot recover the coefficient energy. The energy ledger (THM-4667
statements 1-5) makes this global:

    E_(d+1) - E_d = INTRA_d + INTER_d,     INTRA_d = -2 Σ_(s ≡ 1 (8)) c_d(s) c_d(2s + 1),

where INTER_d collects the cross terms between different blocks of one fibre (scale ratio at least32).
When every initial source is rooted, eventually `E=R(N)^2`. Without that assumption, `c_d(1)` still eventually equals
R(N), because its finite initial rooted support has a finite maximum entry depth. Thus always
`R(N)^2<=liminf_d E_d(N)`. The octave consequence is equally unconditional. Equality of energy with squared partial
currents is restricted to the rooted coefficient support; THM-4667 statement2 gives the missing complement.

## 3. The sibling covariance carries almost all of the energy change

Octave sums of the energy, `/X^2`, are from `root_transport_sibling_anatomy_20261017.py`:

| octave | total INTRA | total INTER | INTRA after the trace minimum | INTER after |
|---|---|---|---|---|
| `2^10` | -5.94 | -0.28 | +0.22 | 0.00 |
| `2^11` | -5.42 | +0.09 | +0.46 | +0.10 |
| `2^12` | -2.88 | -0.43 | +2.17 | -0.17 |

- **Depth by depth.** At `2^12` the INTRA share of the change is:
  - 97%, 100% and 94% at `d = 0, 1, 2`;
  - still the bulk at `d = 3..5`.
- **The termination climb.** The climb back from the trace minimum, from 1.98 to 3.98 at `2^12`, is entirely sibling
  covariance with the opposite sign.
- **Sampled `N` near `2^16`.** INTRA before the minimum is -4.3 to -5.3 per `N`. INTER is between -0.1 and +0.44.

The ideal dissipation is therefore anti-correlation of the partial currents through sibling seeds. That is the transport
form of S37's "two-octave memory, exact but cancelled" and of HYP-9284's odd-step suppression.

## 4. The ideal regime: exact ladder, rate, non-contraction, Fourier form

**Ladder (FINITE-EXACT; THM-4667 statement 6).** The diagonal factor is `(1/2)(3/2) Σ 4^(-j) = 1`, so only cross terms
change the energy.

| d | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `E_d` | 4.857143 | 3.388403 | 2.575382 | 2.137430 | 1.915079 | 1.778419 | 1.693330 |
| cross term | | -1.4687 | -0.8130 | -0.4380 | -0.2224 | -0.1367 | -0.0851 |
| `γ_d` | 2.285714 | 1.784095 | 1.448878 | 1.228070 | 1.079843 | 0.967086 | 0.877335 |

| d | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|
| `E_d` | 1.634126 | 1.590908 | 1.556295 | 1.527640 | 1.501754 | 1.479358 |
| cross term | -0.0592 | -0.0432 | -0.0346 | -0.0287 | -0.0259 | -0.0224 |
| `γ_d` | 0.802930 | 0.738891 | 0.682908 | 0.633171 | 0.588440 | 0.548051 |

Depths 7-12 come from the banded computation in clock_duality's positive gauge
(`root_transport_energy_ladder_banded_20261017.py`, `b < 200`), with a rigorous band-tail bound below `3·10^(-25)`.
- Its `γ_d` reproduce clock_duality's certified `γ_7..γ_10` to all 12 printed digits.
- Its `E_d` reproduce the integer ladder for `d <= 6`.
- `γ_11` and `γ_12` are new.

- **Rate.** Consecutive cross-term ratios are 0.55, 0.54, 0.51, 0.62, 0.62, 0.70, 0.73, 0.80, 0.83, 0.90, 0.87. The
  decay is slower than geometric: the local power-law exponent of `|cross_d|` falls from about 2.4 at `d = 6..8` to
  values between 1.1 and 1.8 (noisy) at `d = 10..12`.
- **Correction.** A first version of this note extrapolated a geometric rate from `d <= 6` and claimed `E_∞ ≈ 1.55`.
  Already `E_10 = 1.528` lies below that. Whether `E_d^ideal` tends to a positive limit is OPEN.
- **Sanity check.** The actual transport averages `E_6/N` = 1.693 (octaves `2^18..2^22`), and `E_0(N)/N -> 34/7`.

**Not a contraction (VERIFIED; `root_transport_symbol_norm_20261017.py`).** Weight band `b` by `2^(-b/2)`. The ideal
transport becomes an analytic Toeplitz operator with matrix symbol
`T(z) = 2^(-1/2) Σ_(j,a) (-1)^a (z/sqrt 2)^(a - 1 + 2j) Π_(j,a)`; its norm is `sup_(|z|=1) ||T(z)||`.
- On the last-digit-balanced subspace (invariant by R8.2) the norm is 1.283, 1.324, 1.335, 1.346 for `d = 1..4`. On the
  full space it is 1.447.
- So the operator expands some balanced vectors. Dissipation is a property of the actual vectors `K_(d-1)`: the
  correlation structure of the backward trees.
- This rules out the "contraction plus additive matching errors" route. That route would have given
  `sqrt(E_∞) <= sqrt(E_0) + Σ_d ||err_d||`, so it has to be done on the cone of actual vectors.

**Fourier form (PROVED; THM-4667 statement10, repaired scope).** Write `o=a-1+2j`, so `o+1` is the homogeneous exponent of the primitive affine branch. In additive characters mod `3^(d+1)`:

    K̂_d(ξ, b) = Σ_(o=0)^(b) (-1)^(o+1) e(-μ_o ξ / 3^(d+1)) K̂_(d-1)(2^(-(o+1)) ξ mod 3^d, b - o),
    μ_o = 2^(-(o+1)) - (1 - 4^(-j))/3.

- **Primitive affine characteristic function.** The phases compose maps
  `P_A(r)=(3r+(3+2^a-2^A)/3)/2^A`, not the ordinary `(3r+1)/2^A` for the same word.
  After the sign gauge, the letters have geometric(1/2) law and the sign is constant on each band slice.
- **Exact comparison.** The positive operators satisfy `P=C A`, `S=A C`, `P^d=C S^(d-1) A`, where S is ordinary
  geometric Syracuse and A,C retain the boundary and band data; see the cited positive-transport note, section2.
  At A=3 the primitive and ordinary translations are -1/8 and +1/8, respectively. An estimate for the ordinary
  offset alone does not automatically control this joint band/frequency energy. That source explains the precise
  scope of the cited Tao estimate; no new literature implication is asserted here.
- **Integer clock.** If the initial and terminal braid indices are j_0,j_d, the physical clock is
  `Σ A_i+2j_0-2j_d`. For11 ->17 ->13=B(3) it is3, while the canonical exponent sum is5.

**Two-point form (PROVED, elementary).** In clock_duality's positive gauge the ideal transport is a Markov kernel `Π`
on (residue, band): `(r, b) -> (F_A(r), b + A - 1)` with probability `2^(-A)`. Then `ν_d = ν_0 Π^d` and

    E_d^ideal = Σ_(x, x') ν_0(x) ν_0(x') G_d(x, x'),      G_d(x, x') = 4·2^d Σ_y 2^(b_y) Π^d(x, y) Π^d(x', y).

- `G_d` is a weighted coalescence kernel. It measures how often two **independent** copies of the primitive affine chain,
  started at `x` and `x'`, sit at the same residue mod `3^(d+1)` and the same band at time `d`.
- `ν_0` is last-digit balanced, so any part of `G_d` independent of the starting residue cancels. The ideal energy
  measures only how fast two chains forget their starting residues mod 3, with matched exponent sums.
- So the `l^2` estimate behind HYP-9288 (C), clock_duality's question 2 in energy form, is a **two-point coalescence**
  estimate for the joint offset/exponent law. A one-point characteristic-function bound such as Tao's is not enough.
- The natural precedent in the repository is S20's coupled-orbit covariance (THM-4565, "two readers of one tape").

## 5. The actual transport for every N < 2^23

The initial source vector changes in only three entries when `N -> N + 1`, and the transport maps point masses to signed
point masses. So all depths can be maintained at once, by walking one path of the odd-transport tree per update
(`root_transport_incremental_trace_20261017.py`). The cost for every `N < 2^23` is 161 s, with 120,409,136 tree nodes.
`R(N) = H(4N) - H(N)` was checked at every `N`.

- **Octave traces** `T_d(X) = Σ_(X <= N < 2X) E_d(N)`, for every octave `X = 2^2, ..., 2^22` (and `2^23`, from the
  run to `N < 2^24`: 240,950,362 nodes, 219 s, invariant checked at every `N`):
  - the maximum over `d` is at `d = 0`, with `T_0/X^2 -> 51/7 = 7.2857`;
  - the final value equals `D_X/X^2`;
  - the only exception is `X = 2`, where `T_5 = 4 > T_0 = 2`.
- **Final values `D_X/X^2`.** For `2^9..2^23`: 1.02, 1.06, 1.95, 3.98, 2.95, 0.66, 2.04, 0.78, 0.99, 1.60, 0.84, 0.21,
  0.23, 4.54, 1.51. The closest approach to `T_0` is `2^22`, at 0.62 of it.
  - In octave `2^23` every profile again peaks at `d = 0`, and the shell maximal sum is `2.11 X^2`.
- **Profile.** The octave profile per `N` is universal in the ideal regime:

  | d | 0 | 3 | 6 | 10 | 15 |
  |---|---|---|---|---|---|
  | `E_d/N` | 4.857 | 2.137 | 1.693 | about 1.53 | about 1.43 |

  From about `d = 15` it becomes octave-specific (1.39-1.54 across octaves at `d = 15`).
- **Pointwise.**
  - `max_(N < 2^24) sup_d E_d(N)/N = 10.06`, at `N = 7910`. There `R^2/N = 9.08`, which is also the largest
    `R(N)^2/N` for `2^10 <= N < 2^23`.
  - In the octaves `2^18..2^21` (3.9 million values of `N`) every energy profile has its maximum at `d = 0`.
  - The octaves with large current energy have a share of `N` whose energy climbs above `E_0` at some depth: `2^12`
    13%, `2^13` 34%, `2^22` 17% (maximum 6.27).
- **The termination regime is ROOT accumulation.** For `N` = 7910, 10870 and 66463, ROOT holds this share of the energy:

  | depth | 20 | 40 | 60 |
  |---|---|---|---|
  | ROOT share | 40% (8% at 66463) | 75-89% | 92-98% |

  - The participation ratio is 1.0-1.8 from `d = 40` on.
  - Since `E_d = Vroot[d]^2 + Σ_(u ≠ 1) c_d(u)^2`, the late energy of these `N` is essentially the square of the ROOT
    coefficient.
  - When the final current is small the picture differs. At `N = 116599` (`R^2/N = 0.018`) ROOT holds 27-36% at
    `d = 40..60`, and the remaining energy (0.4N at `d = 60`) is spread over about 10 effective seeds before it cancels.
  - These sampled large-current rows are ROOT-dominated; no uniform implication from large current to ROOT dominance is proved.
- **The shell maximal inequality holds in every computed octave.** This is `Σ_N max_d Vroot_d^2`, with the maximum
  inside the sum. It is stronger than the uniform shell bound `sup_d Σ_N Vroot_d^2` that HYP-9288 (B) implies, and is
  observed, not implied. The script also records `max_d Vroot_d(N)^2` for every `N < 2^23`.

  | octave | `2^9` | `2^10` | `2^11` | `2^12` | `2^13` | `2^14` | `2^15` | `2^16` | `2^17` | `2^18` | `2^19` | `2^20` | `2^21` | `2^22` |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | `Σ_N max_d Vroot_d^2 / X^2` | 3.85 | 1.90 | 3.21 | 5.20 | 4.73 | 1.88 | 3.50 | 1.69 | 2.44 | 2.74 | 1.87 | 1.10 | 1.03 | 4.85 |
  | ratio to `D_X` | 3.77 | 1.79 | 1.64 | 1.31 | 1.60 | 2.85 | 1.72 | 2.16 | 2.48 | 1.71 | 2.22 | 5.18 | 4.45 | 1.07 |

  - Every octave stays below `T_0 = 7.29 X^2`, even with the maximum over `d` inside the sum.
  - The Doob ratio exceeds the martingale-model constant 4 only at `2^20` and `2^21`, where `D_X` is tiny (0.21 and
    0.23). There the shell partial sums overshoot and then cancel at the end.
  - By THM-4667 statement 9, the ROOT coefficient is the partial sum of the current over the sources with at most `d`
    odd transports: `Vroot[d] = Σ_(depth(x) <= d) c_0(x) f(x)`.

## 6. Merging what we have with where we need to go

1. **Ideal regime.** Target: ideal dissipation, or summable positive parts of the total ideal merge increments (INTRA
   and INTER), for all `d`.
   - Statement 10 of THM-4667 gives an `l^2` target for the primitive joint band/frequency law. The classical Syracuse
     comparison must preserve the A, C boundary kernels and the band; marginal offset mixing is insufficient.
   - The dissipation persists to `d = 12` but decays more slowly than geometrically (ratios 0.83-0.90 at `d = 10..12`).
     So there is no evidence of a spectral gap, and summability of the cross terms is itself OPEN.
2. **Matching (integer clocks).** At fixed depth the actual coefficients differ from the ideal ones only near affine
   thresholds. In `l^2`, the target is matching errors that add in quadrature.
   - `root_transport_l2_matching_20261017.py` compares the actual `c_d(t)` with clock_duality's ideal
     `K_d(t mod 3^(d+1), b(t))` at **every** seed, for `d <= 6` and `N = 2^12..2^18`.
   - For `N = 2^16, 2^18` the coefficients differ at **ROOT and at most two other seeds**:
     - the seed `3^d` (9, 27, 81, 243, 729). For dyadic `N` it lies exactly on an ideal band boundary
       `M_d N/2^b = 3^d 2^(k+2-d-b)`;
     - occasionally one more small seed (15 at `d = 5`, 91 at `d = 6`).
   - The total `l2` error is at most `0.002 N`. Off ROOT it is at most `0.0014 N`, and `l1` error/N falls with `N`.
   - To depth 10 (`root_transport_l2_matching_deep_20261017.py`, ideal coefficients recovered exactly from the banded
     positive gauge), the coefficients still differ only at ROOT and at most three small seeds near ROOT (3, 7, 17, 19,
     51, ...).
     - Off ROOT the `l2` error stays at most `0.0014 N`.
     - The ROOT part grows with depth, up to `0.013 N` at `d = 7`, `N = 2^16`. It is the absorbed mass.
   - So at fixed depth the whole matching problem is the ROOT accumulation, i.e. the termination regime, plus `O(1)`
     boundary seeds. This finite evidence merges items 2 and 3. Section 6b records clock_duality's P3, a uniform `l1`
     version to `d <= m/16`.
3. **Termination regime.** In the sampled large-current rows the late energy is ROOT-dominated;
   `Vroot[d]` is exactly a shell partial sum.
   - HYP-9288(B) implies the uniform shell bound `sup_d Σ_N Vroot_d^2` (not the stronger maximal inequality of
     section 5). The converse requires a separate bound for the
     remaining coefficient energy; the finite dominance observations do not supply it.
   - Its liminf bounds the current-energy budget, without assuming universal convergence. Its content is the census and explicit-formula machinery:
     THM-4659/4660, the zeros of `D`, and the near-return structure that produces the octave-specific profiles.

**Calibration.** With THM-4658, `sqrt(e_(4X)) <= (sqrt(e_X) + sqrt(δ_X))/2 + sqrt(3/X)`. So a budget `δ_X <= a + b log X`
for the current energy gives a budget of the same type for the census energy. The following certify the Mercator pole
(HYP-9283(B)):
- constant budgets up to 1806.6;
- logarithmic budgets written as `a + b log_2(X/2^37)` inside the certified region, e.g. `(2, <= 1194.9)`;
- in the form `a + b log_2 X`, only `b <= 45.7` (2^37 certificate) or `b <= 338` (2^40), as the S41 audit corrected.

Either way this forces `limsup log B(X)/log X >= 0.99821` along a sequence of scales. So the averaged cone (B), with such
a budget, is beyond Krasikov-Lagarias.

## 6b. Concurrent work: clock_duality's positive transport (read after this note's first commit)

[root_positive_transport_20261017g.md](root_positive_transport_20261017g.md) (clock_duality_oct15, committed at 19:31,
before this note's first commit at 20:07, but read only afterwards) attacks the same gap in `l1`.
- **P1 (positivity).** In the sign gauge `ν_d = (-1)^(b+d) K_d/(4·2^(d+b))`, the ideal transport is a positive
  probability kernel. So `γ_d` is non-increasing and the `l1` norm is never amplified.
  - Their sharp one-step hostile (`δ_0 - δ_1` in one band) rules out a uniform one-step contraction on balanced band
    measures in `l1`.
  - It is the `l1` counterpart of the operator norms 1.28-1.35 of section 4. In both norms, dissipation must come from
    the actual vectors.
- **Syracuse bridge.** Their equation (5), `P = CA` with `S = AC` the random Syracuse affine kernel, is the measure form
  of THM-4667 statement 10. The two were found independently in the same hour; statement 10 is the additive-Fourier form.
  - Their caveat applies here too. Tao's Proposition 1.14 concerns the offset alone, while the band readout keeps the
    exponent. What the energy needs is therefore an `l^2` estimate for the joint offset/exponent law, which is their
    question 2.
- **P3 (uniform remainder).** P3 proves `|A^(d)(N) - γ_d N| <= 87(m + 1) 3^(m/4)` for every `d <= m/16`
  (`m = ceil(log_2 4N)`). That is the matching item 2 of section 6, in `l1`, on a logarithmic depth range.
- **P4 (slow packet).** P4 shows `A^(d)(N) >= N 2^(-d) - 12`. The outer word (all exponents 1) blocks `l1` envelopes
  below `d ≈ (log_2 N)/2`.
  - In `l2` the same packet costs only `O(N 2^(-d))` energy: it has `O(1)` coefficients on about `N 2^(-d)` columns.
  - So it is invisible to the energy at scale `N` beyond the first few depths. This is one concrete advantage of the
    energy reading.

## 7. What fails, and what was corrected during the session

- **l1 contraction.** Fails: a singleton unmerged terminal keeps its mass (clock_duality).
- **Weighted `l^2` contraction for all vectors.** Fails for pure power weights `s^β`. The first fibre element
  `(2u - 1)/3` lies below `u`, and the fibre extends to all scales `4^j u`.
- **Operator contraction on balanced vectors.** Fails, by section 4 (norms 1.28-1.35).
- **An absolute pointwise constant.** Not expected, so the pointwise form of HYP-9288 was weakened before testing:
  - in the random-merge model `R(N)^2/N` has chi-square-like fluctuations with about one effective value per octave;
  - so `sup_N R(N)^2/N` would grow like `log log N`.

  The statement now carries a logarithmic factor, and the averaged form a logarithmic budget, which is all the
  calibration needs. A first draft stated `sup_d E_d(N) <= 8.33 N` from 320 sampled `N`. The complete scan to `2^23`
  gives 10.06 at `N = 7910`.
- **A geometric rate of ideal dissipation.** Extrapolated from `d <= 6` (ratio about 0.6, `E_∞ ≈ 1.55`). Refuted by
  `d = 7..12` (section 4).

The incoming integration also supplies [outer packets O6--O7](outer_packet_golden_clocks_20261017h.md): the isolated
singleton packet has squared mass `(3/2)N*2^(-d)+O(18)`, compatible with O(N) energy, and source rounding costs at
most `4 floor(log_4(4N))` in square-root energy at **all depths**. On a dyadic source it is an exactly known ROOT atom.
This pays source rounding, while discrete-to-density matching remains open.

## Reproduction

```
python 04-computation/experiments/root_transport_depth_ladder_20261017.py 6
python 04-computation/experiments/root_transport_energy_ladder_banded_20261017.py 12 200
python 04-computation/experiments/root_transport_crossterm_anatomy_20261017.py 5
python 04-computation/experiments/root_transport_symbol_norm_20261017.py 4 360
python 04-computation/experiments/root_transport_incremental_trace_20261017.py 22
python 04-computation/experiments/root_transport_sibling_anatomy_20261017.py 10 11 12 --sample 16 12
python 04-computation/experiments/root_transport_energy_profile_20261017.py 10 17 40
python 04-computation/experiments/root_transport_octave_trace_20261017.py 9 10 11 12
```
