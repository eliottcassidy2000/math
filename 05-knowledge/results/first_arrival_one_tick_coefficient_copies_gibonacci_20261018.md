# First arrivals into the 627-seed frontier, coefficient copies of the six-colour recipe, and Fibonacci from (5,1) and (7,1)

2026-10-18, mac-mini-2026-10-18-sixcolour. Canon: [THM-4690](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md),
[HYP-9325](../hypotheses/HYP-9325-one-tick-law-of-first-arrivals-into-the-627-seed-frontier.md). These are renumbered from
THM-4670 / HYP-9292, which are ceded to opus-2026-10-18-S43. The three pre-registration files committed before the
renumbering still use the old numbers.

**Owner prompt.**
- Keep working on first-arrival cancellation into the 627-seed frontier, and on connecting coefficient copies for
  six-colour forcing.
- Pull in recent repo work and reroute creatively, preserving sources, ports and integer clocks.
- Merge in prior work on the Solfeggio frequencies.
- Consider the values the Fibonacci sequence would take if it started 5,1 or 7,1 instead of 1,1.

**Pulled in.**
- clock_duality_oct15:
  - R5/R6 and the certified forward-closed frontiers (root_sparse_boundary_module_20261017i.md, sections 8-9);
  - the cubic golden recipe and its open target, "a coefficient context that connects these copies"
    (golden_contact_circuit_20261017i.md section 4).
- S41: the signed ROOT transport and energy ledger (THM-4667).
- S42:
  - the frontier energy cone HYP-9290 (Lemma F, basin form, P1-P16);
  - the Solfeggio/golden mod-9 dictionary L4 (preprint_connections_solfeggio_digits_20261017.md section 4).

## Results at a glance

| # | Result | Type |
|---|---|---|
| 1 | Backward ladder `V^(k)`, port-diagonal energy, merge identities | PROVED (THM-4690(1-2)) |
| 2 | Integer-clock Gram `E_V = G(0) + 2 Σ G(D)`; `E_V` is the clock spectrum at the Nyquist frequency | PROVED (THM-4690(3)) |
| 3 | Exact source ledger: `E_0 = 34N/7` (constant KNOWN from THM-4667's ideal ladder; proved here for the actual sources), sibling-source pairs `5N/7` (one tick apart); hence `E_V = 24N/7 + (non-sibling pair sums)` for every stopping set | PROVED (THM-4690(4)) |
| 4 | Port-orthogonality (cross-port term decays): Q1, Q2, Q4, Q5 FAILED, Q3 passed | FINITE-EXACT, refuted reading |
| 5 | One-tick law `G(1)/G(0) ≈ -0.28`: T1-T4 PASSED at 2^26 | FINITE-EXACT, HYP-9325(A) |
| 6 | Pair-sign laws and higher-lag bound: L1-L7 PASSED at 2^27 | FINITE-EXACT, HYP-9325(B, C) |
| 7 | (A) + (B1) + (C) ⇒ `E_V <= 2.60N`, so `|R(N)| <= 40.4 sqrt(N)` | PROVED implication; the hypotheses are OPEN |
| 8 | 19-adic obstruction: the coefficient copies, and all 49 sheets of the recipe sum, are separated for every path whose steps are integral at the eight `h ≡ 7` primes above 19, which covers every power of the contact unit `u` | PROVED (THM-4690(5)) |
| 9 | Connector census inside `F`: no 2-step connector; no 3-step connector through 450 named directions (integral, `u`-powers, K-directions with poles above 19) | FINITE-EXACT |
| 10 | Fibonacci from (5,1) and (7,1): values, closed forms, `Z[φ]` factorisations `φ^-4(φ^7 - 1)` and `φ^-3 √5 (φ^5 - 1)`, Cassini constants `-L_7` and `-5L_5` | KNOWN / PROVED (elementary) |
| 11 | Solfeggio merge: the ROOT denominators 7 come from the Solfeggio three-lift period (`4^3 - 1 = 63`); mod 9, (5,1) is 4·Fibonacci (the linear part of `B`) and (7,1) is Fibonacci shifted by 9 | PROVED |
| 12 | Solfeggio prime factors against the (5,1)/(7,1) arithmetic | NUMEROLOGY (p ≈ 0.07, look-elsewhere not corrected) |

## 1. First-arrival cancellation into the 627-seed frontier

### 1.1 Setting

`R(N) = H(4N) - H(N)` is carried by the signed seed transport `P` (THM-4667). Every source is stopped at its first seed
in a stopping set `V`:
- `V_1000`: clock_duality's certified forward-closed frontier, 627 seeds, largest 63665;
- `V_100`: the same construction, 74 seeds;
- `H_1000`: the height frontier, the 375 seeds `<= 1000`.

By Lemma F / R6.4, `R(N)^2 <= |V| E_V(N)`. So a bound `E_V = O(N)` is exactly HYP-9290's energy cone.

Sources are aggregated on seeds with coefficient `c_0`, and keep:
- their **port**, the last outside seed;
- their **depth**;
- their **literal integer clock** `t`, the difference of shortcut stopping times. The sign of an arrival is `(-1)^t`.

All quantities below are computed at every `N` with the invariant `R(N) = H(4N) - H(N)` checked at every `N`:
- `fap.c`: every `N < 2^27`, 3.5 GB;
- `fap_lean.c`: every `N < 2^28`, 2.7 GB;
- `fap_levels.c`: the ladder at chosen `N`.

### 1.2 Ports first: the port-orthogonality reading fails

The port-diagonal energy `E^(1) = Σ_(u∈V) c_0(u)^2 + Σ_ports J(p)^2` drops every cross term between different ports
(THM-4690(2)).

The first guess was that the first-arrival cross term `X = E_V - E^(1)` decays. This was pre-registered as Q1-Q5 before
the scan of 2^23-2^25:

| Item | Prediction | Observed | Result |
|---|---|---|---|
| Q1 | V_1000 octave maxima of \|X\|/N <= 0.10 | 0.124, 0.155, 0.060 | FAILED |
| Q2 | mean of those maxima < 0.08 | 0.113 | FAILED |
| Q3 | V_1000 octave maxima of E_V/N <= 2.5 | 1.33, 1.23, 1.42 | PASSED |
| Q4 | V_100 octave maxima of \|X\|/N <= 0.35 | 0.232, 0.364, 0.435 | FAILED |
| Q5 | V_1000 octave means of X/N in [-0.06, 0.06] | -0.076, -0.085, +0.010 | FAILED |

`X/N` stays in `[-0.66, +0.11]` from 2^11 on but does not decay.

The sibling ports `(s, 2s+1)` are nearly uncorrelated at level 0: the correlation of the signed arrivals is between
-0.17 and +0.05 at the dumps. So the ports are not where the cancellation lives.

### 1.3 The backward ladder: cancellation is spread over depth

`V^(k) = {dep <= k}` is forward closed. Its stopped energy `E^(k)` interpolates from `E_V` at `k = 0` to the source
energy `E_0` at maximal depth (THM-4690(1)).

At `N = 6,904,411`, `E^(k)/N` takes the values
`1.19, 1.18, 1.26, 1.33, 1.38, 1.39, 1.40, 1.45, ..., 3.00 (k = 25), 3.40 (30), 4.03 (40), 4.40 (50), 4.61 (60)` and
reaches 4.857 at the maximal depth 241.
- The bulk, 3.2N of the total 3.7N, is removed at levels 6-60.
- Each level from about 6 to 40 removes 0.04-0.12 N, typically 80-100% of it the sibling merge `(s, 2s+1)` at that
  level.
- Levels 0, 4 and 5 remove almost nothing; levels 1-3 remove 0.05-0.08 N.

At `N = 16384` the band sits at levels 0-20. **The dissipation band recedes from the frontier as N grows.** A
"frontier-local" mechanism, ports or a few levels, therefore cannot be the whole story. That is the reason to turn to
the integer clocks, which see all depths at once.

### 1.4 Integer clocks: the one-tick law and the Nyquist dip

Group arrivals at each `u` by literal clock: `h_t(u)`, with `G(D) = Σ_u Σ_t h_t h_(t+D)`. Then
`E_V = G(0) + 2 Σ_(D>=1) G(D)` exactly.

On V_1000, at the dump `N = 6,904,411`, the values of `G(D)/N` for `D = 0..12` are
`3.305, -0.914, -0.093, +0.039, -0.050, +0.040, -0.040, -0.013, -0.041, +0.069, -0.024, -0.023, +0.040`.
The Gram is essentially tridiagonal: one large negative lag-1 term and small noise.

The per-`N` ratio `G(1)/G(0)` stays in `[-0.288, -0.256]` for every `N` in `2^14..2^27`.

**Pre-registered T1-T4**, octave 2^26: all passed.
- Mean ratio -0.2803.
- Per-`N` range `[-0.2845, -0.2744]`.
- `G(0)/N` mean 3.086.
- `E_V` max 1.277, and cross-port mean -0.024.

**Spectral reading** (THM-4690(3)). The unsigned clock polynomial `Q_u(z) = Σ c_0(x) z^(t(x))` gives the clock spectrum
`S(θ) = Σ_u |Q_u(e^(iθ))|^2`, with:
- `S(π) = E_V` (signed);
- `S(0) = E+` (unsigned arrivals);
- the mean of `S` equal to `G(0)`.

Measured:
- `E+/N` has octave means 4.1-5.4, close to the source energy 4.857: unsigned arrivals do not cancel.
- `E_V/E+` has octave means 0.23-0.33 from 2^14 on.
- At `N = 10^8`, `S(kπ/12)/N` is `4.85, 4.01, 4.11, 4.60, 4.09, 3.81, 3.63, 2.80, 2.04, 1.91, 1.52, 1.33, 1.18`.

**First-arrival cancellation into the 627-seed frontier is the Nyquist dip of the integer-clock spectrum.** The signs
`(-1)^clock` evaluate a low-pass clock profile at the alternating frequency.

### 1.5 The exact source ledger (THM-4690(4))

**Source energy, and where the 7 comes from.**
- Every odd `m <= 4N` lies on one B-chain `x, 4x+1, 16x+5, ...`, whose residues mod 3 advance by 1.
- The weights are `κ = (0, -3, 3)` below `2N` and `λ = (-1, -1, 2)` on `[2N, 4N]`.
- So `c_0(x)` depends only on `x mod 3`, the chain length mod 3, and whether the chain ends in `[2N, 4N]`, and
  `|c_0| <= 3`.
- Summing the geometric scales, with period 3 in the chain length because `4^3 - 1 = 63 = 9·7`:

  `E_0(N) = (34/7) N + O(log N)`.

**Sibling sources.** Siblings `x ≡ 1 (mod 8)` and `2x + 1` have the same image under `P`. Their step clocks are
`2 + 2r` and `3 + 2r`, so they always arrive **one tick apart with opposite signs**, at a common entry seed, in every
stopping set. Their correlation is

`SIB(N) = Σ c_0(x) c_0(2x+1) = (5/7) N + O(log N)`.

**Prior art in the repo.** `34/7 = 4.857143` is the `d = 0` value of THM-4667's ideal energy ladder. THM-4667(7)
observed the octave-trace limit `51/7 = (3/2)(34/7)`. New here: a proof for the actual sources, the sibling constant
`5/7`, and the pair ledger. THM-4667's depth-1 cross term -1.4687 is, by clock lag, `-10/7` (siblings, about 97%) plus
small depth-one co-preimage terms at lags 5, 6, 7 (+0.0179, -0.0759, +0.0179; observed to match `2/112` and `-17/224`).

Both constants are derived symbolically and checked by brute force (`source_ledger_constants.py`):
`|E_0 - 34N/7| <= 11.2` and `|SIB - 5N/7| <= 4.9` up to `N = 4·10^6`. `fap.c` reproduces them at every `N`.

**Consequence: a pair ledger valid for every stopping set,**

`E_V = 24N/7 + (G(0) - E_0) + 2Γ_1 + HL + O(|V| + log N)`.

Here `G(0) - E_0` collects same-clock pairs, `Γ_1 = G(1) + SIB_V` the other lag-1 pairs, and `HL` lags `>= 2`.

The budget on V_1000 at octave 2^27, per `N`:

| Term | Value | Status |
|---|---:|---|
| Sources | +4.857 | exact |
| Sibling-source pairs | -1.429 | exact |
| Same-clock pairs | -1.778 | |
| Other lag-1 pairs | -0.310 | |
| Longer lags | -0.245 | |
| **Total, `E_V/N`** | **1.095** | |

For `V = {1}` the same identity reads `R(N)^2 = 24N/7 + (pair sums of the full transport) + O(log N)`.

### 1.6 Pair-sign laws (HYP-9325) and the 2^27 test

On V_1000, H_1000 and V_100, at every `N >= 2^11`:
- the same-clock sum is negative (`G(0) < E_0`);
- the other-lag-1 sum is negative (`Γ_1 < 0`).

On the two large frontiers the longer lags add at most `+0.068N` (V_1000) and `+0.061N` (H_1000) from 2^11 and 2^12
on. On V_100 they reach `+1.29N`: small frontiers keep the termination-regime climbs of HYP-9290 in their long lags.

**Pre-registered L1-L7**, octave 2^27, V_1000. Committed and pushed before the run, as
`first_arrival_pair_ledger_prereg_20261018.json`. All seven passed:

| Item | Observed |
|---|---|
| L1 | ratio mean -0.2823, per-`N` range `[-0.2873, -0.2773]` |
| L2 | same-clock max -1.758 / mean -1.778 |
| L3 | `Γ_1` max -0.137 / mean -0.155 |
| L4 | HL max -0.134 |
| L5 | `E_V/N` max 1.183, mean 1.095 |
| L6 | `G(0)/N` mean 3.079 |
| L7 | dip mean 0.248 |

`E_0/N` and `SIB/N` matched 34/7 and 5/7 to five digits at every `N`.

**Meet-depth anatomy** (`pair_series_meets.py`, FINITE-EXACT at `N = 2^18`, V_1000, meets up to depth 16 below the
meeting node). A pair of sources first meets at a node `y`, `k_1` and `k_2` steps above them, before first entry into `V`.
- **Unequal-depth meets** (`k_1 ≠ k_2`) are nearly uncorrelated: each `(k_1, k_2, D)` cell is at most 0.013 in absolute
  value, and the displayed cells sum to `+0.03`.
- **Equal-depth meets carry the cancellation.**
  - Lag 1: the terms are `-0.712` (siblings), `-0.307`, `-0.124`, `-0.049` at depths 1-4. They are slightly *positive*
    at depths 6-14, `+0.005..+0.018` per level, about `+0.14` in total. So `Γ_1` shrinks as deeper levels fill in,
    matching the observed trend. HYP-9325(B2) is expected to fail eventually; the route below does not use it.
  - Same clock: about `-0.04..-0.06` per level at depths 2-8, decaying to about `-0.01` by depth 15.

  In the full transport (`V = {1}`), the depth-2 cell sums per `N` agree to four digits from `N = 2^18` and to about five from
  `N = 2^20` up to `4·10^6`: -0.04673 (lag 0), -0.30903 (lag 1), -0.06858 (lag 2), +0.02832 (lag 3)
  (`pair_series_meets_full_K2.out`). So the per-depth constants
  exist, as the siblings' `-5/7` does at depth 1. Small-denominator fits such as `-89/288` at lag 1 are suggestive
  only; no derivation has been made.

  An equal-depth, equal-clock pair has the same number of 3-divisions and the same total power of 2 on both paths. Such
  "clock twins" sit at the same scale, so their chain parameters are locked, as for the siblings in the 34/7 and 5/7
  derivation. This is why a residue-table computation per depth is the natural attack on (B1).

**What this buys.** `E_V = (1 + 2G(1)/G(0)) G(0) + HL`. So HYP-9325(A), with ratio `<= -0.24`, together with (B1),
`G(0) <= E_0`, and (C), `HL <= 0.07N`, gives

`E_V <= 0.52·34N/7 + 0.07N <= 2.60N`, hence `|R(N)| <= sqrt(627 · 2.60 N) ≈ 40.4 sqrt(N)`.

This is HYP-9290's cone for the closed frontier, reduced to sign and ratio statements about pair correlations that carry
their integer clocks.

**What is not proved.** None of (A), (B) or (C).
- The weakest observed trend is `Γ_1`, which is shrinking: -0.44 at 2^13, -0.15 at 2^27. The meet anatomy explains it
  by positive deep lag-1 terms. The route above does not need `Γ_1 <= 0`.
- `G(0)/N` decreased 3.97 → 3.08 over 2^15 → 2^27, but only by 0.007 in the last octave.
- A proof of (B1) would be a residue-class computation of same-clock coincidences in the spirit of the 34/7 derivation.

## 2. Connecting the coefficient copies of the cubic recipe (six-colour forcing)

### 2.1 The question

clock_duality's cubic recipe has factors `A + aM`, `C + bM` and `C + cM`. Here:
- `M` is the Moser spindle, `C` the 10-point coefficient set, `a = conj(d)`;
- `b, c` are the coefficients of the `u`-square circuit, with denominator 19.

Each 19-adic factor is seven disjoint translated copies of `C`. Each copy is coloured independently, so forcing cannot
propagate between copies. The recorded target was a coefficient context that connects them.

### 2.2 Why short connectors cannot exist: 19-adic separation (THM-4690(5))

**The primes above 19.**
- The tower order `Z[ζ15][j][w]` is 19-maximal.
- 19 has 24 primes in `F`, all with residue field `F_361`, in three types of eight according to `h = w + 3 ≡ j, -j, 7`.
  This follows from `h^3 - 7h^2 + 11h - 1 ≡ (h^2 + 11)(h - 7) (mod 19)`.

**Valuations** (`nineteen_places.py`).
- `b` and `c` have simple poles at all 24 primes, and `b/c` is a 19-unit.
- The following are all 19-integral: `C`, `A`, `a`, the spindle and its differences (units at every prime), the 30
  roots of unity, and the mixed-golden unit `u_2`.
- The contact unit `u = (h - j)/(h + j)` has valuation -1 at the eight `h ≡ -j` primes, +1 at the eight `h ≡ j`
  primes, and **0 at the eight `h ≡ 7` primes**.

**Consequences.**
- A path between two copies `C + bm`, `C + bm'` must, at **every** prime above 19, contain a step with a pole there.
  The displacement `b(m' - m) + (integral)` has a pole at every prime.
- For the whole recipe sum, every nonzero sheet move `(Δm_2, Δm_3)` has a pole at some `h ≡ 7` prime: 1166 of 1224 at
  all eight, the rest at six or seven. No move is 19-integral (`diagonal_19adic_check.py`).
- So the 49 sheets `(m_2, m_3)` are mutually unreachable by paths whose steps are integral at the `h ≡ 7` primes. This
  covers every direction `ζ^l u^k q` with `q` 19-integral, at every length.

**A connecting context must use a unit direction with a pole at an `h ≡ 7` prime.**

### 2.3 Census inside F (FINITE-EXACT)

`copy_connectors.py` and `copy_connectors_screen.py`:
- **Two steps.** A connector with middle point in `F` exists iff `N(t)(N(t) - 4)` is a square in `F` and `|t| <= 2`.
  All 172 + 118 cross-copy point pairs with `|t| <= 2` fail: each is a certified non-square, one quadratic non-residue
  at a degree-1 place where the element is a unit.
- **Three steps, one step from a named set.**
  - 90 19-integral directions (roots of unity times the unit differences of `C`): 17854 + 12340 tests, 3 screen
    survivors, all certified non-squares.
  - 120 19-adic directions `ζ^l u^k`, `k = ±1, ±2`: 23729 + 16508 tests, 10 survivors, all certified non-squares at
    `p = 3631`, `4051` or `4261`.

  - 240 K-directions `ζ^l δ_𝔮^(±1)` with poles above 19, as in 2.4, one `δ` per conjugate pair of primes of `K`:
    47515 + 33031 tests, 37 survivors of the 12-place screen, all certified non-squares at the 80-place screen.

2.2 explains the second family completely: `u` cannot touch the poles at the `h ≡ 7` primes. The third family has the
right kind of pole, but one `δ` covers only 3 of the 24 primes.

### 2.4 Next target (concrete)

Every prime of `K = Q(ζ15, √-11)` above 19 is split from its complex conjugate. So `δ = conj(γ)/γ` with
`γ ∈ 𝔮 ∖ conj(𝔮)` is a unit vector in `K` with a pole at the three primes of `F` above `𝔮`, including an `h ≡ 7` one.

A first probe, with one `δ` from each conjugate pair and the exact two-step test for the remainder, found nothing
(2.3). The next search should combine several `δ`'s so that the step poles can cover all 24 poles of `b(m' - m)`. That
means four or more steps, or `δ`'s attached to products of primes. Its outcome decides whether the recipe can carry cross-copy forcing at all.

Any connector must also be checked against the saved additive witnesses: the 570 port patterns and their contact
families, HYP-9289. A new point with at most five neighbours inside one factor never blocks a 6-colouring by itself.

## 3. Fibonacci started at (5,1) and (7,1), and the Solfeggio merge

### 3.1 The values

With `G_0 = a`, `G_1 = b`, `G_(n+1) = G_n + G_(n-1)`:
- **(5,1)**: 5, 1, 6, 7, 13, 20, 33, 53, 86, 139, 225, 364, 589, 953, 1542, 2495, 4037, 6532, 10569, 17101, 27670,
  44771, 72441, 117212, ...

  This is 5 followed by OEIS **A022096** ("Fibonacci sequence beginning 1, 6").
- **(7,1)**: 7, 1, 8, 9, 17, 26, 43, 69, 112, 181, 293, 474, 767, 1241, 2008, 3249, 5257, 8506, 13763, 22269, 36032,
  58301, 94333, 152634, ...

  This is 7 followed by OEIS **A022098** ("Fibonacci sequence beginning 1, 8").

### 3.2 Structure

**Closed forms.**
- In general `G_n = a F_(n-1) + b F_n = (a L_n + (2b - a) F_n)/2`.
- (5,1): `G_n = 5F_(n-1) + F_n = (5L_n - 3F_n)/2`.
- (7,1): `G_n = 7F_(n-1) + F_n = 4L_n - L_(n+1)`, an integral Lucas form, because `5 | 2b - a`.

**Consecutive sums** (KNOWN; OEIS comments by McRae 2014 and del Arte 2013).
- (5,1) is the sum of **seven** consecutive Fibonacci numbers: `G_n = F_(n-5) + ... + F_(n+1) = F_(n+3) - F_(n-4)`.
- (7,1) is the sum of **five** consecutive Lucas numbers: `G_n = L_(n-4) + ... + L_n = L_(n+2) - L_(n-3)`.

**In `Z[φ]`** (PROVED, elementary). A Gibonacci sequence is the orbit of `θ = G_1 + G_0(φ - 1)` under multiplication by
`φ`, with norm `N(θ) = b^2 - ab - a^2`, the Cassini constant.
- `θ(5,1) = 5φ - 4 = φ^-4 (φ^7 - 1)`, which is the Fibonacci orbit times `φ^7 - 1`.
- `θ(7,1) = 7φ - 6 = φ^-3 √5 (φ^5 - 1)`, which is the Lucas orbit times `φ^5 - 1`.

So the start 5 pairs with a window of 7, and the start 7 pairs with a window of 5. Since `N(φ^r - 1) = -L_r` for odd `r`:
- the Cassini constants are `-29 = -L_7` and `-55 = -5·L_5`;
- `G_n^2 - G_(n-1)G_(n+1) = (-1)^n · 29`, resp. `(-1)^n · 55`.

**Growth.**
- `G_n/F_n → b + a/φ`: 4.0902 for (5,1) and 5.3262 for (7,1).
- `G_n = round(α φ^n)` from `n = 4`, resp. `n = 5`, with `α = (b + a/φ)/√5`:
  - `(25 - 3√5)/10 = 1.8292` for (5,1);
  - `(7 - √5)/2 = 2.3820` for (7,1).
- The conjugate weight is `conj(α)`, so `α·conj(α) = -N/5`: 29/5 and 11.

### 3.3 Arithmetic (KNOWN classical facts, checked)

**Primes dividing the Cassini constant.** If `p` divides it, the sequence mod `p` is an eigen-sequence `G_n ≡ a λ^n`,
where `λ` is a root of `x^2 - x - 1`. So `p` never divides a term:
- 29 never divides a (5,1) term (`λ = 6`, period 14);
- 5 and 11 never divide a (7,1) term (period 4 mod 5, 10 mod 11).

**Periods.** At all other moduli the periods are Pisano's. (7,1) has exactly the Lucas Pisano periods (Purath 2019,
OEIS A022098).

**Primes.** Prime terms of (5,1): 5, 7, 13, 53, 139, 953, 44771, ... Of (7,1): 7, 17, 43, 181, 293, 13763, ...

**Squares.**
- Fibonacci itself: 1 and 144 (Cohn).
- (5,1): `1` and `225 = 15^2`.
- (7,1): `1`, `9 = 3^2` and `3249 = 57^2` (n = 15).

These lists are FINITE-EXACT for `n < 1500`. Finiteness is KNOWN (Siegel): `G_n = y^2` gives `z^2 = 5y^4 ± 4N` with
`z = 2G_(n+1) - G_n`. A complete determination was not attempted.

### 3.4 The Solfeggio merge

**The mod-9 dictionary** (S42 L4).
- The nine Solfeggio numbers list the three fibres of `Z/9 -> Z/3` in their digits.
- Each fibre is a `B^3`-orbit, since `B^3(x) = 64x + 21 ≡ x + 3 (mod 9)` for the seed map `B(x) = 4x + 1`.
- The golden classes are the `Q`-orbits `{±1}, {±2}, {±4}`.

What the two new starts do in this dictionary (PROVED):
- **(7,1) is Fibonacci shifted by 9 mod 9.** `(7, 1) ≡ (F_9, F_10) (mod 9)`. So its digital roots are Fibonacci's,
  shifted, the golden class is ±1, and the Solfeggio digits 3, 6, 9 occur exactly at `n ≡ 3 (mod 4)`.
- **(5,1) is the Fibonacci sequence pushed through the linear part of `B`.** `G_n ≡ 4 F_(n+14) (mod 9)`, golden class
  ±2. Multiplication by `4 = 2^2`, the linear part of `B`, fixes each Solfeggio fibre and cycles the golden classes (L4).
  The digits 3, 6, 9 occur exactly at `n ≡ 2 (mod 4)`.

  Its 24-cycle of digital roots is 5,1,6,7,4,2,6,8,5,4,9,4,4,8,3,2,5,7,3,1,4,5,9,5.
- **The ROOT denominators 7 are Solfeggio.** The source constants `34/7` and `5/7` of section 1.5 have denominator 7
  because the chain residues have the `B^3` period. The geometric sums run over `4^(-3k)`, and `4^3 - 1 = 63 = 9·7`.
  This is the one place where the Solfeggio three-lift structure enters the cancellation budget exactly.

**NUMEROLOGY** (recorded honestly).
- The prime factors above 3 of the nine Solfeggio numbers are 5, 11, 13, 19, 29, 71, 107 and 139.
  - 29 is the (5,1) Cassini prime `L_7`, and divides 174.
  - 11 `= L_5` divides 396 and 528.
  - 139 and 13 are terms of (5,1); they divide 417 and 741.
  - Six of the eight primes split or ramify in `Q(√5)`, against Chebotarev's 1/2.
- The overlap between the Solfeggio primes and the prime terms of (5,1) below 150 is 3 of 8. A hypergeometric tail
  gives `p ≈ 0.07`, before any correction for the many features one could have compared.
- No mechanism links these.

## 4. Files

All scripts are in `04-computation/experiments/first_arrival_20261018/`. The `.out` files are receipts.
- **ROOT.**
  - `fap.c`: ports, sibling pairs, depth, integer clocks, clock spectrum, pair ledger; every `N`, invariant checked.
  - `fap_lean.c`: the pair ledger in 5 bytes per source; the 2^27 test.
  - `fap_levels.c`: the backward ladder.
  - `source_ledger_constants.py`: 34/7 and 5/7, symbolic and brute force.
  - `ledger_table.py`: octave tables.
  - `V_closed_Y{3,10,100,1000}.txt` (from clock_duality's `frontier_renewal.json`) and `V_height_Y1000.txt`.
- **Six-colour.**
  - `nineteen_obstruction_check.py`, `nineteen_places.py`, `diagonal_19adic_check.py`: the valuations and sheet moves.
  - `copy_connectors.py` and `copy_connectors_screen.py`, with environment variables `ROOT_OF_UNITY` / `U_DIRECTIONS`
    selecting the direction families.
  - `shim/numba.py`: a no-op `njit` so the S41 helpers import without numba.
- **Fibonacci.**
  - `gibonacci_51_71.py`.

Pre-registrations, in `05-knowledge/results/`:
- `first_arrival_port_orthogonality_prereg_20261018.json` (Q1-Q5);
- `first_arrival_one_tick_prereg_20261018.json` (T1-T4);
- `first_arrival_pair_ledger_prereg_20261018.json` (L1-L7).

```bash
cd 04-computation/experiments/first_arrival_20261018 && cc -O2 -o fap_lean fap_lean.c && ./fap_lean V_closed_Y1000.txt 27
```

```bash
cd 04-computation/experiments && PYTHONPATH=first_arrival_20261018/shim python3 first_arrival_20261018/nineteen_places.py
```

## 5. Process notes

- **Numbering.** THM-4670 and HYP-9292 were reserved here with stub files at 09:29, but not entered in the INDEX.
  opus-2026-10-18-S43 claimed both numbers through the INDEX at 09:34/09:39 for the Clausen-twin result. They are
  ceded; this work is THM-4690 / HYP-9325.
- **Q1-Q5 failed.** This is kept as a refuted reading, and it redirected the work from ports to integer clocks.
