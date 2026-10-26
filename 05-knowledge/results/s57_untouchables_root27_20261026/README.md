# S57: the untouchable share against the Pollack–Pomerance model, odd untouchables, tc fibres to 11 vertices, the merged sheet tables, and the root-27 drawdown law

**opus-2026-10-26-S57.** Owner dispatch (verbatim in the session letter):
- merge the prior work on tables mod odd `n` alternating between the `3n-1` and `3n+1` graphs;
- pursue a source-preserving rank for arbitrarily changing tree shapes on the root-27 orbit, with exact carries and
  guards, hostile examples, and the Frobenius-defect and clone-merger machinery where useful;
- consider `{2,3,5,7} = {3,5} | {2,7}` (15 and 14);
- find creative tricks to find or rule out a second odd untouchable number, or the exact share of untouchables;
- relate the `{7,21}` tournament work to mapping a tournament class to a pair of not necessarily connected graphs;
- consider `F_1` with the owner's list `0 <- 1`, `1 <- primes`, `2 <- nothing`, `3 <- 4`, ..., `13 <- {27, 35}`.

**Status (typed; independently audited, [AUDIT](AUDIT.md), MISTAKE-672, corrections applied).**
- **PROVED:**
  - THM-4790: no tournament on `n <= 11` vertices has `tc = H/|Aut|` in `{7, 21}`, by a reduction for all `n`, an
    automorphism-type lemma and a finite computation. This extends HYP-9352 (2) from `n <= 9`. The audit
    re-implemented it independently, and runs that bypass the lemma (types `(5)` and `(3)(3)` at `n = 11`) also
    find nothing.
  - THM-4792: the merged sheet tables, a reference theorem.
    - The centre labels `-+1/2` are the reflection centres of `ceil(3x/2)` and `floor(3x/2)`.
    - The unit shift (the root-27 section change) conjugates the two one-step maps.
    - The halving-shift alternation of the edge tables holds iff `m` divides some `2^t + 1`. That holds iff every
      prime `p | m` has `ord_p(2)` even with one common 2-adic valuation, which is **KNOWN** (Moree 1997, Acta
      Arith. 80, Thm 1; OEIS A014657) and answers the exercise left open in kind-pasteur's Theorem A.
  - Elementary here: the `F_1` fibre comparison (sec. 4). The kill/serve rule of the prime-extension cores (sec. 2.4)
    is a **classical mechanism**, the same residue exclusion as in PP 2016 sec. 2 and Pomerance–Yang 2014 Prop. 2.3.
    What is new is the exact core decomposition of the data.
- **VERIFIED (numerical, three independent evaluations):** Pollack–Pomerance's model constant is
  `Delta = 0.17865689 (+- 1e-8)`.
  - This session's window averages give `0.17865689`.
  - The audit's deterministic convolution of the independent-valuation model gives `0.17865689139`.
  - Monte Carlo agrees to `1e-5`.
  - PP (sec. 3.2) note that `Delta` is effectively computable and that their printed numbers are intermediate values
    of their (3.4); they did not evaluate the limit. Their printed `0.171822` at `y = 2e10`, and CGJM's partial sums
    to `2^40`, are reproduced to 6 digits by `S(y)/log y = Delta - 0.1621075/log y`.
- **FINITE-EXACT:**
  - exact untouchable census to `10^10`, equal to A057978 at every power of ten;
  - residue, abundancy and core-family profiles;
  - odd Goldbach-dependence to `10^8`;
  - the tc enumeration;
  - the root-27 carry word to `2·10^6` symbols.
- **Pre-registered:** 13 checks in two commits (`568cb689a4`, `4369352bc8`); **13 pass** (sec. 8).
- **CITED:**
  - Pollack–Pomerance 2016 (Trans. AMS Ser. B 3, secs. 2 and 3.1–3.2, read in the AMS viewer);
  - Chum–Guy–Jacobson–Mosunov (CGJM, arXiv:2110.14136, Table 2);
  - Pomerance–Yang 2014;
  - OEIS A057978, A005114 and A014657;
  - Oliveira e Silva–Herzog–Pardi (Goldbach to `4·10^18`);
  - Moon 1972 and Busch 2006 (strong minimum);
  - Moree 1997;
  - for the density of `n` with `l` not dividing `sigma(n)`: Ramanujan (unpublished manuscript, sec. 19), Watson
    1935 (Math. Z. 39), Rankin 1961, Scourfield 1964, Ciolan–Languasco–Moree 2023;
  - Sylow's and Jordan's theorems.
- **OPEN (untouched):** existence and value of the untouchable density; a second odd untouchable (Goldbach);
  root-27 nonhalting; HYP-9352 (2) for `n >= 12`.
- **HYP-9400** records the density predictions. **THM-4791, HYP-9401 and HYP-9402 are released**: the root-27 rank
  was done by two concurrent lanes (sec. 7), and the tc statement already is HYP-9352 (2).

## 1. Inheritance and the concurrent lanes (read before drafting)

Three packets on this same dispatch landed during the session; nothing below re-proves them.
- **kind-pasteur-2026-10-25-S4**, [sheets_root27_rank_untouchables_20261025.md](../sheets_root27_rank_untouchables_20261025.md):
  - Theorem A: negation and the halving shift of the modular edge tables;
  - Theorem B: the rank `rho_z(w) = |w| - v_2(Q_z(w))` depends only on the first mismatch with the carry word
    `kappa(z)`, and the Haar halting share is `sqrt5 - 2`, with klein's `phi^-6` for the blank source;
  - Theorem C: the Goldbach reduction for odd untouchables, the census to `10^6` and the 11.25% Goldbach-dependence.
- **session66-root27**, [REPORT](../session66_root27_aliquot_graph_interfaces_20261026/session66-root27/REPORT.md):
  - the switch cocycle `B_P D_P N_Q(y) = A_P D_Q N_P(z) + B_P Delta_(P,Q)`;
  - zero defect iff `P`, `Q` are powers of one primitive word (the exact clone-merger rule);
  - the source pullback `h_P(z_n) + n = v_2(s - rho_(w,P))`;
  - an actual `+7` refill at macro 86 of the root source;
  - unbounded refills at a fixed shape pair (no shape-weighted fuel is a rank);
  - the Frobenius defect stated as *not* the same object.
- **session66-aliquot / session66-root**, [aliquot REPORT](../session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md)
  and [divisor bridge](../session66_root27_aliquot_graph_interfaces_20261026/session66-root/DIVISOR_F1_BRIDGE.md):
  - the owner's list is the complete fibres `s^-1(m)`;
  - `n <= (m-1)^2`, a core-prime / squarefull certificate, Booker's hyperbola trick;
  - a certified Chen–Zhao instantiation, `liminf U(X)/X > 0.0604529`;
  - `U(X) = E(X) + o(X)`;
  - the monoid carrier `A_a(p) = s(n)`, `A_a(1) = tau(n) - 1`, with the no-go `4 -> 3` versus `9 -> 4`;
  - the 14 interior divisors of 210 with complement defect `7` and the pair `(14, 15)`.

Also inherited:
- my 2026-10-05 [stitched-sheets note](../collatz_spectrum_equality_and_stitched_sheets_20261005.md), Props. 1–3;
- the [root-27 packet](../root27_machine_bridges_20261026/anchor/README.md);
- THM-1370, THM-4741 and HYP-9352 (presentation counts);
- THM-4744 (Moon/Busch minimum);
- THM-2450 sec. 2 (stack multiplicativity).

Corrected near misses in this session:
- A first non-semiprime source enumerator stopped its prime loop at the first `k > N`, which is wrong at
  `p = spf(m)`. There `s(pm)` can be *smaller*, as for `n = 7^2 * 229 -> 1889`. Brute force to `6000` caught it, and
  the guard is in the script.
- The pre-audit draft made four mistakes (MISTAKE-672):
  - it tested PP's assumption with a sufficient proxy (`l` not dividing `sigma(n)`) instead of the assumption itself;
  - it attached thin-set exponents `(log k)^(-1/(l-1))` to spreads whose fitted decay is near `(log k)^(-0.75)`;
  - it presented the Moree classification and the kill/serve mechanism as new;
  - it missed the single-source numbers 3 and 7.

Canonical hostiles used:
- the deficient end of PP's per-number law;
- `q + 1` for prime powers `q` (sec. 4);
- source 15 and the halting sources 5, 9, 13 (sec. 7).

## 2. The share of untouchable numbers

### 2.1 Exact census

[`untouchable_sieve_20261026.py`](../../../04-computation/untouchables_s57/untouchable_sieve_20261026.py) records the
multiplicity `|s^-1(k)|` of every even `k <= N`. The parity lemma gives the preimages of even `k`: even `m <= 2k-2`,
odd prime squares `q^2` (`k = q + 1`), and odd composite squares with `q <= N^(2/3)`. Odd `k`: only 5 below
`4·10^18` (Goldbach, sec. 3). Counts to `10^10` equal A057978 at every power of ten: `2, 5, 89, ..., 165,826,606,
1,681,871,718` (707 s for `10^10`). The census reproduces known data and is not a record (Chum et al. reach `2^40`).

### 2.2 The Pollack–Pomerance constant is 0.17865689

PP's Conjecture 1.4 (sec. 3.1) sets `Delta = lim (1/log y) sum_(a <= y, 2|a) e^(-a/s(a))/a` and proves the limit
exists.
- **What PP say.** At `y = 2·10^10` "the expression inside the limit has value `0.171822`", and their counts to
  `10^10` are "consistent with a limiting density of approximately 0.17". In sec. 3.2 they add that `Delta` is
  effectively computable, that their printed numbers are intermediate values of (3.4), and that they did not
  evaluate the limit.
- **Quoted value.** The value quoted downstream as "about 0.17" (OEIS A005114, erdosproblems #418, Booker) comes
  from those intermediate values. The limit itself is not stated in PP, CGJM, OEIS, erdosproblems #418 or Booker,
  per the search here and the audit's.

Since `sigma(a)/a` has a limiting distribution, `Delta = (1/2) E[e^(-1/h)]` with `h = s(a)/a` and `a` a random even
integer. Evaluations ([`pp_delta_20261026.py`](../../../04-computation/untouchables_s57/pp_delta_20261026.py)):
- **window averages** of `(1/2) e^(-a/s(a))` over even `a` in `[10^k, 10^(k+1))`: `0.17865689` at `k = 5, 6, 7, 8`
  (agreeing to `1e-8`), and `0.1786578` at `k = 3`;
- **the audit's deterministic convolution** of the independent-valuation (Erdős–Wintner) model, at two prime cutoffs
  and two grids: `0.17865689139`;
- **Monte Carlo** of the same model (primes `<= 2000`, `4·10^7` samples): `0.178644 +- 0.000011`, plus `2.07e-5` for
  the primes above 2000;
- **partial sums** `E(y) = S(y)/log y`:
  - at `y = 10^6 .. 10^9` they are `0.1669232, 0.1685994, 0.1698566, 0.1708344`, equal to
    `Delta - 0.1621075/log y` within `7e-8` (within `3e-8` for `y >= 10^7`);
  - the formula gives `0.1718224` at `y = 2·10^10`, which is PP's printed value;
  - it gives `0.171617, 0.172257, 0.172790, 0.172810` at `10^10, 10^11, 10^12, 2^40`, which is CGJM's Table 2.

So the model predicts `Delta = 0.17865689`, and its partial sums converge only like `1/log y`. The exact finite-`y` class weights `phi(A_y/a)/A_y` (summing to `1/2`) give the truncations
`Delta_y = 0.0677, 0.1064, 0.1365, 0.1499, 0.1568, 0.1628, 0.1655, 0.1683, 0.1703, ..., 0.1772` at
`y = 2, 3, 4, 5, 7, 8, 9, 11, 13, ..., 53`. They increase to `Delta` with gap about `0.08/y`
([`pp_truncated_20261026.py`](../../../04-computation/untouchables_s57/pp_truncated_20261026.py)).

### 2.3 The data are below the model and approaching slowly

Local untouchable density (all integers) per decade:

| decade | 10^5–10^6 | 10^6–10^7 | 10^7–10^8 | 10^8–10^9 | 10^9–10^10 | 10^10–10^11 | 10^11–10^12 |
|---|---|---|---|---|---|---|---|
| density | 0.1515 | 0.1583 | 0.1630 | 0.1662 | 0.1684 | 0.1701 | 0.1713 |

The last three columns come from A057978. Per dyadic window the density goes `0.16515, 0.16601, 0.16679, 0.16750,
0.16812, 0.16870` at `W_28..W_33` (`W_e = (2^(e-1), 2^e]`), gaining about `0.0006–0.0009` per octave and
decelerating. The gap to `Delta` shrinks from `0.027` to `0.0074` between the first and last decades.
- A `1/log^2 N` fit of the last decades extrapolates to `0.1773–0.1778`.
- A `1/log N` fit extrapolates to `0.1845–0.1866`.
- `gap = 5.35/log^2 N` with limit `Delta` fits the last three decades within `0.0005`.

CGJM compare `U(X)/X` with the partial sum at `y = X` and see the two converge together, which is consistent with
this. The data are compatible with `Delta` as the limit and do not decide it. Matching each window to a truncation gives
`y_eff = 8.1, 9.5, 11.2, 12.8, 15.2` at `10^8 .. 10^12`, growing a little faster than `log N`.

### 2.4 Two ways the real map departs from PP's random map at computable sizes

PP's heuristic assumes that **`s` maps `T_a = {n : gcd(n, A_y) = a}` into `T_a`**, "asymptotically true if
`n > e^(e^y)`", and is otherwise a random map that multiplies by `s(a)/a`. The target's class then fixes the Poisson
mean `a/s(a)`. There are two exact measures of the distance from that model
([`pp_class_preservation_20261026.py`](../../../04-computation/untouchables_s57/pp_class_preservation_20261026.py),
[`pp_assumption_failure_20261026.py`](../../../04-computation/untouchables_s57/pp_assumption_failure_20261026.py)).

**(a) The literal assumption.** `gcd(s(n), A_y) = gcd(n, A_y)` fails at the prime `l` iff
`min(v_l(s(n)), e) != min(v_l(n), e)`, where `l^e || A_y`. For `l` not dividing `n` this means `sigma(n) = n (mod l)`.
On `n` in `(2^26, 2^27]`:
- per prime (uncapped), the failure share is `0.104, 0.177, 0.241, 0.175, 0.150, 0.124` for `l = 2, 3, 5, 7, 11, 13`;
- **jointly**, `gcd(s(n), A_y) != gcd(n, A_y)` for **62% of `n` at `y = 13` and 82% at `y = 53`**;
- for even `n` the joint shares are 64% and 86%.

So at these sizes PP's class assumption fails for most inputs, as their `n > e^(e^y)` qualifier anticipates. The
audit corrected the pre-audit figures, which had used criterion (b) below in place of (a).

**(b) Residue neutrality.** A random map into `T_a` treats the nonzero residues of `s(n) (mod l)` alike. The real
map does so on inputs with `l | sigma(n)`, where `s(n) = -n (mod l)`. The complementary share is
`f_l = share(l not dividing sigma(n)) ~ C_l (log X)^(-1/(l-1))` (Ramanujan, Watson 1935; Rankin 1961; Scourfield
1964). Measured on `(X/2, X]` to `X = 2^31`:

| `l` | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|
| `f_l`, `n` near `2^31` | 0.135 | 0.521 | 0.569 | 0.816 | 0.797 |
| fitted exponent of `log X` (`2^10` to `2^31`) | -0.525 | -0.258 | -0.167 | -0.101 | -0.083 |
| `-1/(l-1)` | -0.500 | -0.250 | -0.167 | -0.100 | -0.083 |

**Neither (a) nor (b) is small at any computable size, so the data cannot test PP's limit directly.**

Departure (b) has an exact signature in the target's residues `k (mod l)`. PP's class model sees `k` only through
`gcd(k, A_y)`, so it gives every residue `r != 0 (mod l)` the same share. The identity `s(aP) = s(a)P + sigma(a)` for
a prime `P` not dividing `a` gives:

**Kill/serve rule (elementary; classical mechanism).** This is the residue exclusion PP use in sec. 2 (proof of
Thm. 1.3, with `s(2^k m q r) = r s(2^k m q) + sigma(2^k m q)`), and Pomerance–Yang use in 2014 Prop. 2.3.
- If `l` does not divide `s(a)`, the core `a` (preimages `n = aP`) cannot serve `k = sigma(a) (mod l)` once
  `k > l s(a) + sigma(a)`, because then `l | P = (k - sigma(a))/s(a)`. It serves the other residues with local factor
  `l/(l-1)`.
- If `l | s(a)`, the core serves only `k = sigma(a) (mod l)`.
- Cores with `l | sigma(a)` treat all nonzero residues alike. Only the cores with `l` not dividing `sigma(a)`
  distinguish them.

[`core_decomposition_20261026.py`](../../../04-computation/untouchables_s57/core_decomposition_20261026.py) splits
every preimage of every even `k` in `[10^7, 10^8]` by core `a = n/P(n)`. The zeros and exclusive services appear
exactly as the rule says: core 2 (`k - 3` prime) misses `k = 3 (mod 5, 7)`, core 14 serves only `k = 4 (mod 5)`, and
core 8 serves only `k = 1 (mod 7)`. The audit checked the rule on every preimage in `[1e5, 1e6]` with `P > 13` and
found no violation.

Mean multiplicity by residue against the untouchable shares at `W_26`:
- **mod 5:** multiplicities `1.18, 1.05, 1.14, 1.15, 1.68`, shares `0.32, 0.38, 0.34, 0.35, 0.25`. Inverse order,
  with one near-tie.
- **mod 7:** multiplicities `1.15, 1.49, 1.38, 1.06, 1.24, 1.18, 1.17`, shares `0.34, 0.31, 0.28, 0.38, 0.31, 0.33,
  0.33`. Broadly inverse: rank correlation `0.77` over the nonzero residues, with three discordant pairs. Residue 1
  has the highest multiplicity but not the lowest share, since the share also depends on the dispersion.

The cores `a = 2p` alone kill `k = 1, 2, 0 (mod 5)` for `p = 1, 3, 4 (mod 5)` and serve only `k = 4` for
`p = 2 (mod 5)`. Spreads of the share over the nonzero residues at `W_33`:
- `1.453` (mod 5; argmax 1, argmin 4);
- `1.290` (mod 7);
- `1.281` (mod 11);
- `1.201` (mod 13);
- `1.092` (mod 3).

For `l = 5, 7, 11, 13` they fall from `1.92`, `1.73`, `1.90` and `1.43` at `W_14`. Over `W_24..W_33` the fitted
decay of `(spread - 1)` is like `(log k)^(-0.71..-1.09)`. That is small-core-dominated, about `1/log k`; the thin-core
exponent `-1/(l-1)` that the rule suggests for the eventual regime is not visible in the data. Mod 3 is not monotone:
the spread falls to `1.0095` at `W_18`, where the argmax flips from residue 1 to residue 2, then rises every window to
`1.092`.

The extremes mod 15 have the same cause. `k = 6 (mod 15)` has share `0.539`, because `5 | k-1` kills `(k-1)^2`,
`3 | k-3` kills `2(k-3)`, and core 4 is absent; `k = 4 (mod 15)` has share `0.174` (U-P3).

The **per-number law** is a separate discrepancy, not explained here. PP's model, read per target, predicts share
`e^(-x)` for `x = k/s(k)`. Measured shares and PP's `e^(-x)`:

| `x` bin | 0.2–0.4 | 0.4–0.6 | 0.6–0.8 | 0.8–1.0 | 1.0–1.2 | 1.2–1.4 | 1.4–1.6 | 1.6–1.8 | 1.8–2.0 |
|---|---|---|---|---|---|---|---|---|---|
| `W_16` | 0.392 | 0.369 | 0.333 | 0.304 | 0.269 | 0.239 | 0.183 | 0.239 | 0.212 |
| `W_24` | 0.462 | 0.432 | 0.389 | 0.354 | 0.298 | 0.278 | 0.219 | 0.220 | 0.226 |
| `W_34` | 0.506 | 0.468 | 0.418 | 0.374 | 0.315 | 0.288 | 0.220 | 0.220 | 0.224 |
| PP `e^(-x)` | 0.690 | 0.593 | 0.497 | 0.400 | 0.334 | 0.271 | 0.216 | 0.184 | 0.141 |

`W_34` is the partial window `(2^33, 10^10]`. The abundant end rises toward PP. For `x >= 1.4` the share is **flat at
0.220–0.226 from `W_24` to `W_34`**, where PP's value falls to `0.141`. The residue rule does not produce this excess. A
candidate is departure (a): the class of the preimage changes for most inputs, so a deficient target's preimages do
not come mainly from its own class, but this is not shown.

Class by class at `10^9` (`a = gcd(k, 2520)`), the weighted mean multiplicity is `1.227` against PP's `1.215`, and
the zero share is `0.3324` against `0.3310`. The averages nearly agree while the classes do not (`a = 2`: mean `1.73`
versus `2.00`, share `0.224` versus `0.135`). This is the statement to keep: **PP's constant is a reasonable
average at `10^10`, but its class mechanism is not what the data show.**

### 2.5 Other measured facts (pre-registered)

- Overdispersion (U-P5): at `10^9`, even `k >= 8` have mean multiplicity `1.2286`, variance `1.7466`, and zero share
  `0.3317` against `e^-mean = 0.2927`.
- The two easy preimages are not independent of the rest: the generic-zero share is `0.412` with neither `k-1` nor
  `k-3` prime, `0.347` with only `k-3` prime, `0.461` with only `k-1` prime, and `0.397` with both.
- `210 = 15 * 14` (U-P2): the shares mod 210 factor as `share(mod 15) * share(mod 14)/share`, with median relative
  error `0.059`. Splitting 15 as `3 * 5` gives `0.105` and the full split `3 * 5 * 7` gives `0.122`.
  - For even `k`, "mod 14" carries the same information as mod 7. The pass says the mod-7 profile is nearly
    independent of the mod-15 profile.
  - The kill/serve rule explains why 3 and 5 interact: the extremes mod 15 are coincident kills.
  - It is a residue statement, not a reason to identify 14 and 15 with anything else (MISTAKE-529; THM-4792).

### 2.6 What can be said about the exact share

Nothing proved here moves the bounds:
- lower density at least `0.0604529` (Chen–Zhao, certified by session66);
- upper density below `1/2` (Luca–Pomerance, session66);
- existence of the density open.

The exact content is:
1. PP's model, the only heuristic for the value, predicts `0.17865689`.
2. The observed local density is `0.1713` at `10^12` and rising, consistent with that limit at a `1/log^2`-type rate.
3. At every computable scale the model's class assumption fails jointly for 62–86% of inputs. Its residue
   neutrality fails for a share `~(log X)^(-1/(l-1))` at each small prime. The kill/serve rule accounts for the
   residue structure; the abundancy (per-number) discrepancy is unexplained.

**HYP-9400** states the falsifiable predictions; see it for the ranges.

## 3. A second odd untouchable

Odd `k >= 9` is `s(pq) = 1 + p + q` whenever `k - 1` is a sum of two distinct primes. So a second odd untouchable is
a distinct-Goldbach exception and exceeds `4·10^18` (classical; kind-pasteur Theorem C; session66; Booker's footnote
on distinctness). Additional measured facts
([`odd_goldbach_dependence_20261026.py`](../../../04-computation/untouchables_s57/odd_goldbach_dependence_20261026.py),
exact against brute force to 6000, multiplicities included):

| `N` | odd `k <= N` touched by a non-semiprime source | Goldbach-dependent (only semiprime sources) | share |
|---|---:|---:|---:|
| `10^6` | 443,732 | 56,266 | 11.25% |
| `10^7` | 4,537,267 | 462,731 | 9.26% |
| `10^8` | 45,882,159 | 4,117,839 | 8.24% |

At `10^6` this agrees with kind-pasteur, who count 3 (source 4) as Goldbach-dependent. The audit reproduced the
`10^6` counts and found `8,250` Goldbach-dependent odd numbers at `10^5`, by an independent brute force. Local
decade shares are `9.0%`, then `8.1%`.

**Below `10^8` exactly five odd numbers have a single source of any kind: `3 = s(4)`, `7 = s(8)`, `9 = s(15)`,
`11 = s(21)` and `39 = s(217)`.** Brute force over all `n <= (k-1)^2` confirms this to `2000`. 39 has a single source
because `38 = 7 + 31` and `19 + 19` is not distinct. Every other odd `k != 5` has at least two sources: two
distinct-prime partitions, or a partition plus another source. Verdicts on the tricks:
- **Search:** pointless below `4·10^18`.
- **Rule out:** non-semiprime coverage is not uniform. Its failure set has a share that is decreasing but still 8%,
  and a PP-type bounded mean would leave it positive. No covering argument can replace Goldbach, though failure would
  need both shapes to fail at once.
- **The structural reason 5 is alone:** the only odd `k` whose `k - 1` has no distinct-prime partition are 3, 5 and 7.
  Two of these are Mersenne, `3 = s(4)` and `7 = s(8)`, and 5 is not.

## 4. The F1 dictionary: what is exact

Session66 proved the monoid statement: the divisor enumerator `A_a(X)` evaluates to `s(n)` at the prime weights and
to `tau(n) - 1` at `X = 1`. Two exact additions:
- **Projective reading.** `prod (1 + p + ... + p^a_p) = |prod P^(a_p)(F_p)|` and `n = |prod A^(a_p)(F_p)|`, so `s(n)`
  counts the complement of the big cell in a product of projective spaces, each over its own prime field. In
  particular `s(p^a) = |P^(a-1)(F_p)|`, the hyperplane at infinity. The owner's list then reads:
  - `s(4) = 3 = |P^1(F_2)|`, `s(9) = 4 = |P^1(F_3)|`, `s(25) = 6 = |P^1(F_5)|`;
  - `s(8) = 7`, the Fano plane `P^2(F_2)`;
  - `s(49) = 8`, `s(121) = 12`, and `s(27) = 13 = |P^2(F_3)|`;
  - `s(6) = 3*4 - 2*3`, `s(10) = 8`, `s(15) = 9`, `s(14) = 10`, `s(21) = 11` and `s(35) = 13` are boundaries of
    `P^1 x P^1`.
- **Where F1 and F_p agree.** `s^-1(0) = (tau-1)^-1(0) = {1}` and `s^-1(1) = (tau-1)^-1(1) = {primes}`, exactly the
  owner's first two entries. They differ for every `k >= 2`: `(tau-1)^-1(k)` contains every `p^k` and is infinite,
  while `s^-1(k)` is finite (session66 Lemma A) and empty at `k = 2, 5`. **"2 with nothing" is the first place where
  the prime weights matter**, and the `F_1` shadow `tau - 1` has no untouchables at all.

Hostile, which refutes a tempting extension: `2 = |P^1(F_1)|` and `5 = |P^1(F_4)|`, but `q + 1` for odd prime powers
`q = p^a` (`a >= 2`) is untouchable for **785 of 3,661** such `q` up to `10^9`, the first at `q = 17^2`
(`290` untouchable). So "projective lines over non-prime fields" is not a law; the two small cases are numerology.
Also exact: `sigma(14) = sigma(15) = 24`, so `|P^1(F_2) x P^1(F_7)| = |P^1(F_3) x P^1(F_5)|`, which makes
`{2,7} | {3,5}` the unique such split of the four primes.

## 5. Tournament classes as graph pairs, and the fibres 7 and 21 (THM-4790)

A marked Hamiltonian path `v_1 -> ... -> v_n` of `T` defines a labelled graph `G` on `[n-1]`: `{i, j}` is an edge iff
`v_(j+1) -> v_i`. The pair `(G, complement of G)` of backward and forward non-path arcs is the owner's pair of not
necessarily connected graphs. Every graph occurs, and two paths give the same `G` iff an automorphism relates them.
So a class receives exactly `tc(T) = H/|Aut|` pairs (THM-4741), and `sum_T tc(T) = 2^binom(n-1,2)`. THM-1370 says
`H` is never 7 or 21. HYP-9352 (2) says the fibres of this graph-pair chart are never 7 or 21; the THM-4741 census
covered `n <= 9`. **THM-4790 proves it for `n <= 11`.**
1. A minimal counterexample is strong (multiplicativity).
2. Moon/Busch's floor forces `|Aut| >= f(n)/21`, which is `3, 3, 5, 7, 11` at `n = 7..11`.
3. An odd-order group that large contains an element of a listed cycle type. This uses Cauchy, groups of order 9 and
   25, `Z_3 wr Z_3`, and Sylow's and Jordan's theorems for the order-15 case on 11 points.
4. Every strong tournament invariant under each type was enumerated: up to `2^19` per type, 134,922 `|Aut|`
   computations, 256 s. None has `tc` in `{7, 21}`. Non-strong tournaments follow from step 1.

The audit's independent re-implementation reproduced every count, and found no hit without the candidate filter. Its
runs that bypass the lemma also found no hit:
- `n = 9`, type `(3)`, `2^22` invariant tournaments;
- `n = 11`, types `(5)` (`2^23`) and `(3)(3)` (`2^25`).
So `n <= 11` also holds using only Cauchy and the order-9 argument.

The method stops at `n = 12`, where the 3-group case needs `2^24` invariant tournaments.

## 6. The tables mod odd `n`, merged (THM-4792)

The 2026-10-05 note found the centre of every odd multiplication table to be `+-1/4`, swapped by `3x -+ 1`. In
THM-4792:
- the centre *labels* `k = -1/2`, `k + 1 = 1/2` are the fixed points of `3x+1` and `3x-1`, and the reflection centres
  of `ceil(3x/2)` and `floor(3x/2)`;
- negation swaps the two maps, and the unit shift `z = x - 1`, which is the root-27 section change
  (`8 -> 12 -> 18 -> 27` becomes `7 -> 11 -> 17 -> 26`), conjugates them;
- the halving shift between the tables exists iff `-1` is in `<2>` mod `m`, i.e. iff `m` divides some `2^t + 1`;
- that holds iff all `p | m` have `ord_p(2)` even with equal 2-adic valuation. This is **Moree's theorem** (Acta
  Arith. 80 (1997), Thm. 1, with `a = 2, b = 1`), and the moduli are OEIS A014657. kind-pasteur's open exercise is this
  known result, and the proof in THM-4792 reproduces Moree's.

The valuation classes are classical (second supplementary law):
- `p = 3 (mod 8)`: valuation 1;
- `p = 5 (mod 8)`: valuation 2;
- `p = 7 (mod 8)`: valuation 0, never alternate.

So 3 and 5 alternate separately, while `15 = 3 * 5`, the least odd number with two distinct prime factors, does not,
and `7` is the least prime that never alternates. Checks:
- on 1000 odd moduli here, with the literal table shift to 301;
- on all odd `m <= 100,001` in the audit, with the table shift to 1201.

## 7. Root-27: the drawdown law, and what the concurrent lanes settled

**R-P1 (pre-registered): PASS.** On the root-27 ray (`T(z) = (3z+c)/2`, credit `q' = q + 3c - 1`, source
`(q, z) = (1, 7)`), the largest deficit excursion over the first `10^5` V-returns is **20**. The window was
`[15, 35]` and the Haar scale is `ln J/ln phi = 23.9`. It is `10` at `10^3`, `13` at `10^4` and `27` at `3·10^5`
returns. Over `2·10^6` exact symbols (z of 1.17 million bits) the share of times whose future credit drop is `>= d`
matches the Haar law `phi^-d` within 5% for `d <= 16`, then fluctuates (`d = 1..6` within 0.4%). Controls reproduce the
packet's halting trees exactly, with `q0 = 1`:
- source 15: `1111000010000010100000` after 21 macros;
- 5: `1000`;
- 9: `1010000`;
- 13: `1001001010101000000`.
The carry word of 27 runs past 4096 symbols. This tests the root source against the random-source model of
kind-pasteur B.4 and klein; it is not a bound.

**Source-preserving rank.** It is settled as far as exact statements go, by the lanes in sec. 1:
- the rank of any tree shape is `|w| - LCP(w, kappa(7))` (kind-pasteur B);
- shape switches are exact 2-adic cancellations at the original source (session66 (2)–(6));
- no shape-weighted fuel is monotone (session66 Thm 3);
- the unique exact clone merger is "powers of one primitive word" (session66 Thm 2).

One corollary is immediate from B.4. Every carry cylinder of positive credit `q` has Haar halting mass `phi^-q > 0`, so
no rank depending on finitely many 2-adic digits of the source certifies survival. A proof must use that 7 is a
rational integer, not a 2-adic one. **Frobenius defect:** I agree with session66 that the switch defect
`Delta_(P,Q)` is not the orientation defect of session 65, and no map identifying them is known. The only valid
appearance of that machinery in this dispatch is session66's divisor-complement defect `7` on the 14 interior
divisors of 210. THM-4791 is released.

## 8. Pre-registration outcomes

| id | prediction | outcome |
|---|---|---|
| U-P2 | median error of the `15 x 14` factorization `<= 10%` | `5.9%`: **pass** |
| U-P3 | max mod 15 at `6` (`> 0.45`); min at `4` or `14` (`< 0.20`) | `0.539` at 6, `0.174` at 4: **pass** |
| U-P4 | share at `x` in `[0.4,0.6)` `>= 1.5x` share at `[1.8,2.0)` | ratio `2.0`: **pass** |
| U-P5 | variance/mean `> 1.2` | `1.42`: **pass** |
| T-P1 | no strong 10-tournament with `|Aut| > 1` has `tc` in `{7,21}` | none: **pass** (extended to `n <= 11`) |
| R-P1 | deficit excursion at `J = 10^5` in `[15,35]` | `20`: **pass** |
| U-C | count `<= 10^10` = A057978 | `1,681,871,718`: **pass** |
| U-P6 | local density in `W_33` in `[0.1685, 0.1700]` | `0.168696`: **pass** |
| U-P7 | share at `x` in `[1.8,2.0)` in `W_33` in `[0.215,0.235]` | `0.2241`: **pass** |
| U-P8 | share at `x` in `[0.2,0.4)` in `W_33` in `[0.495,0.520]` | `0.5039`: **pass** |
| U-P9 | mod-5 nonzero spread in `[1.42,1.48]`, argmax 1, argmin 4 | `1.4533`, 1, 4: **pass** |
| U-P10 | mod-13 nonzero spread in `[1.17,1.215]` | `1.2014`: **pass** |
| U-P11 | share(2 mod 3)/share(1 mod 3) in `[1.085,1.11]` | `1.0923`: **pass** |

The 1e10 predictions U-P6..U-P11 were trend extrapolations, so their passing confirms smooth drift; it does not
test the mechanism. U-P9 and U-P10 would have failed had the spreads stopped decreasing. The audit reproduced
independently:
- the `W_33` values of U-P6..U-P11;
- the mod-15 extremes;
- the U-P5 moments;
- the arrays against brute force for every even `k <= 2·10^6`.

## 9. Open and next

1. HYP-9400: a `2^40` census with residue and abundancy tallies (Chum et al.'s range) tests every clause. The cheapest
   strong test is the deficient-end share, predicted to stay in `[0.215, 0.230]`.
2. HYP-9352 (2) at `n = 12`: the types `(3)(3)(3)` and `(3)(3)(3)(3)` on 12 points (`2^24`, `2^22`) in C.
3. A target-side heuristic that keeps the core residues, replacing "`s` maps `T_a` to `T_a`" by the kill/serve rule,
   and predicts the residue spreads quantitatively. The weights are known; the thin-core sums are the open part.
4. The Goldbach-dependent share of odd numbers (8.24% at `10^8`): its limit, zero or positive, is undecided by the data.

## 10. Reproduction

All in `04-computation/untouchables_s57/`; big arrays go to `$S57_BIGDIR`.

```bash
python untouchable_sieve_20261026.py 1e10            # 707 s, 5 GB array; = A057978
python pp_delta_20261026.py 1e9 $S57_BIGDIR          # Delta windows, MC, partial sums
python pp_truncated_20261026.py $S57_BIGDIR 53       # Delta_y, y_eff
python pp_residue_drift_20261026.py 1e10 $S57_BIGDIR 28
python pp_xbin_drift_20261026.py $S57_BIGDIR 1e10 31 34 1
python pp_assumption_failure_20261026.py 31          # f_l = share(l not | sigma(n)), 158 s
python pp_class_preservation_20261026.py 27          # literal gcd-class failure per prime and jointly, 18 s
python core_decomposition_20261026.py 1e8            # kill/serve tables, 12 s
python odd_goldbach_dependence_20261026.py 1e8       # 25 s
python tc_fibres_upto11_20261026.py 11               # THM-4790, 256 s
python sheet_merge_20261026.py 2001                  # THM-4792 controls
python root27_drawdown_20261026.py 2e6               # R-P1, drop law, halting controls
```
