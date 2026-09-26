# The landing multiplicity `L` of the thin-divergence recursion: what it costs, why no orbit-blind count removes it, what would, and what the strip entropy already buys

**Status: REASSESSMENT (stepping back, as asked) with PROVED small
statements, FINITE-EXACT probes, and one new hypothesis (HYP-9161). No
Collatz claim. Session `collatz-landing-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/collatz_landing_20260926_probe.py`
-> `.out` (strip entropies, landing multiplicities on actual segments);
`collatz_zeckendorf_20260926_interlock_probe.py` -> `.out` (the owner's
Zeckendorf / odd-prime seed, section 7).

## 0. The bottleneck in one formula

THM-4476's recursion, for a positive one-signed orbit segment, `X`,
`L = log_2 X`, `k = floor(L)`, `theta`:

```text
N(X) <= k + k N(X^(1-theta)) + #F_b(X, theta),          (R)
```

where `#F_b(X, theta) <= C X^(h*) L^(-3/2) 2^(lambda* theta L) e^(s theta L)`
(THM-4499, Lemma 1.4c). The factor `k` in the middle term is the **landing
multiplicity**: a dipper `i` (an index whose orbit falls below
`y_i X^(-theta)` within `k` steps) is charged to its first landing point
`j`, and each landing point can serve at most `k` dippers. Suppose instead
every landing point served at most `M(L)` dippers, `M(L) = L^beta` with
`0 <= beta <= 1`. Then the bootstrap of THM-4499 (section 3) runs with
`c_1 h* = beta + eta` in place of `1 + eta`, and the exponent becomes

```text
a*(beta) = beta lambda*/h* - 3/2 :   a*(1) = -0.986 (THM-4499),   a*(1/2) = -1.243,   a*(0) = -3/2.
```

So the whole prize of the multiplicity attack is `(log X)^(0.514)`: from
`(log X)^(-0.986)` to `(log X)^(-3/2)`, the ballot floor of the no-dip
count itself (THM-4495). Nothing in this direction changes the exponent
`h*`, and nothing here excludes divergent orbits.

## 1. How large the multiplicity can be (corrected after the audit)

**Lemma (one-bit band; found by the audit).** Let `j` be a landing point
with `theta L >= 1` and `y > |b|` along the segment. The step into `j` is a
halving (`y_j < y_(j-1)`), so `y_j = y_(j-1)/2`, and a dipper `i` of `j` has
`y_i 2^(-theta L) in (y_j, y_(j-1)]`, i.e.

```text
y_i in ( 2^(theta L) y_j,  2^(theta L + 1) y_j ] :   every dipper of j lies in one half-open 1-bit band.
```

At most two of any three consecutive orbit values lie in a band `(a, 2a]`:
a halving leaves it, and two odd steps multiply by `(9y + 5b)/(4y) > 2`
(for `b = -1` once `y > 5`). Hence the multiplicity of every landing point
is at most `2 ceil(k/3) <= 2k/3 + 4/3`, for every orbit. So THM-4476's
`k` can be replaced by `2k/3 + 4/3`; this changes constants, not `beta`.

**What residue classes realise.** A first draft of this section claimed
that a climb followed by a drop, or a hover followed by a drop, gives one
landing point `~ L` dippers. The audit refuted both as stated: in a climb
the thresholds `y_i 2^(-theta L)` are `0.585` bits apart while the drop
crosses one bit per halving, so each halving lands about `1.7` climb
indices and a landing point receives at most `2` of them (word
`1^12 0^11`, `L = 34`: the `12` climb indices land on `7` points); a hover
in `[-W, 0]` spreads its indices over several landing points, and the
exhaustive search at `L = 24`, `theta = 0.15` over all hovers of length `12`
in `[-2, 0]` finds at most `9` dippers on one point. The largest
multiplicity found is `13 = 0.54 L` at `L = 24` (a `20`-step hover in
`[-3, 0]` followed by four halvings, `y_0 = 6654596`, all `13` dippers in
one 1-bit band). So multiplicity a constant fraction of `L` is realised
(`beta = 1`), but by hovers in a narrow band, not by climbs, and never
above the lemma's `2k/3 + O(1)`.

On actual orbit segments the multiplicities are small: with `L = 20` and
`theta = 0.05, 0.1, 0.2`, the mean multiplicity is `2.0`-`2.8` and the
largest `5`-`8` on the 5n+1 orbit of 7, the 3n+1 orbits of 27, of
`2^18 - 1` and of `2^19 - 1`, and even on a synthetic hover-then-drop
integer (probe (B)). The random-walk heuristic is that a landing point's
dippers are the indices spent above `y_j 2^(theta L)` before the final
descent, of order `theta L/0.2075 = 5 theta L`, i.e. `O(log L)` at the
`theta_X` of the bootstrap: polylogarithmic, not `L`.

**Scaling probe** (`collatz_landing_20260926_probe2.py` -> `.out`): at the
bootstrap's own `theta_X = 1.05 log_2 L/L` and `L = 20, 30, 40, 60, 80`, on
3n+1 segments from `2^L - 1`, `2^(L-1) + 1` and random odd starts, and on
5n+1 orbits, the mean multiplicity is `0.4`-`0.8` times `theta_X L`
(i.e. `2`-`5`), while the maximum is `0.04`-`0.35` times `L`. The audit
located the maximal landing points: they are not produced by the initial
climb of `2^L - 1` (which leaves `[1, X]` at step `1`) but by later
stretches in which many window points sit in one 1-bit band, exactly the
lemma's configuration, and random starts show the same `max/L`. So the
*average* is `O(theta L) = O(log L)` as the heuristic says, and the
*maximum* is a constant fraction of `L`: HYP-9161 must be stated for the
average, and it is the average that the recursion uses.

```text
L  theta   thetaL | segment                         total    ND     D   landing  maxmult  mean  mean/(thetaL)  max/L
20  0.227   4.54 | 3n+1 from 2^L-1                    43     6    37      18       4   2.06      0.45   0.200
20  0.227   4.54 | 3n+1 from 2^(L-1)+1                51    30    21       7       6   3.00      0.66   0.300
20  0.227   4.54 | 3n+1 from random odd ~2^L         114    61    53      18       6   2.94      0.65   0.300
20  0.227   4.54 | 3n+1 from random odd ~2^(L-3)     121    74    47      18       5   2.61      0.58   0.250
20  0.227   4.54 | 5n+1 from 7 (climbs)               53    53     0       0       0   0.00      0.00   0.000
20  0.227   4.54 | 5n+1 from random odd ~2^(L/2)      94    77    17       9       3   1.89      0.42   0.150
30  0.172   5.15 | 3n+1 from 2^L-1                   134    58    76      26       7   2.92      0.57   0.233
30  0.172   5.15 | 3n+1 from 2^(L-1)+1                97    26    71      23       9   3.09      0.60   0.300
30  0.172   5.15 | 3n+1 from random odd ~2^L         117    55    62      27       9   2.30      0.45   0.300
30  0.172   5.15 | 3n+1 from random odd ~2^(L-3)     115    54    61      28       7   2.18      0.42   0.233
30  0.172   5.15 | 5n+1 from 7 (climbs)              130   125     5       3       2   1.67      0.32   0.067
30  0.172   5.15 | 5n+1 from random odd ~2^(L/2)      80    75     5       2       4   2.50      0.49   0.133
40  0.140   5.59 | 3n+1 from 2^L-1                   184    47   137      31      11   4.42      0.79   0.275
40  0.140   5.59 | 3n+1 from 2^(L-1)+1                70     0    70      26       9   2.69      0.48   0.225
40  0.140   5.59 | 3n+1 from random odd ~2^L         137    28   109      37      11   2.95      0.53   0.275
40  0.140   5.59 | 3n+1 from random odd ~2^(L-3)     137    25   112      33       8   3.39      0.61   0.200
40  0.140   5.59 | 5n+1 from 7 (climbs)              181   177     4       3       2   1.33      0.24   0.050
40  0.140   5.59 | 5n+1 from random odd ~2^(L/2)     100    96     4       3       2   1.33      0.24   0.050
60  0.103   6.20 | 3n+1 from 2^L-1                   287    52   235      65      16   3.62      0.58   0.267
60  0.103   6.20 | 3n+1 from 2^(L-1)+1               303    37   266      66      21   4.03      0.65   0.350
60  0.103   6.20 | 3n+1 from random odd ~2^L         163     0   163      52      12   3.13      0.51   0.200
60  0.103   6.20 | 3n+1 from random odd ~2^(L-3)     326    43   283      65      14   4.35      0.70   0.233
60  0.103   6.20 | 5n+1 from 7 (climbs)              508   443    65      20       7   3.25      0.52   0.117
60  0.103   6.20 | 5n+1 from random odd ~2^(L/2)     230   198    32      14       4   2.29      0.37   0.067
80  0.083   6.64 | 3n+1 from 2^L-1                   294     1   293      85      18   3.45      0.52   0.225
80  0.083   6.64 | 3n+1 from 2^(L-1)+1               303     1   302      93      11   3.25      0.49   0.138
80  0.083   6.64 | 3n+1 from random odd ~2^L         427     0   427      90      16   4.74      0.71   0.200
80  0.083   6.64 | 3n+1 from random odd ~2^(L-3)     353     4   349      86      18   4.06      0.61   0.225
80  0.083   6.64 | 5n+1 from 7 (climbs)              563   468    95      25      11   3.80      0.57   0.138
80  0.083   6.64 | 5n+1 from random odd ~2^(L/2)     191   188     3       1       3   3.00      0.45   0.037
```

## 2. Why the shift structure does not turn into a bound by counting

Consecutive windows are shifts of one word, so the `m` dippers of a landing
point are determined by one residue class modulo `2^(L+m)`. This costs
nothing: the class count of "hover-then-drop" words of length `L` is
`2^(h_W m) 2^(L-m) ...`, about `X 2^(-(1-h_W) m)` integers, which is below
`X^(h*)` only for `m >= (1-h*) L/(1-h_W)`. So long hovers are rare *among
integers* (section 3), but hovers of length a constant fraction of `L`
are not, and a hypothetical orbit may contain many of them as far as
counting can tell. Charging dippers to landing points, to window minima,
to leaders (future minima of the orbit) or to record lows all reproduce
the same `k`. Splitting the orbit into visits below `X` or into
band-visits (maximal runs in a dyadic band) gives the identity
`N(X) = sum over band-visits of their lengths` and the bound
"band-visits of length `>= m` are at most `X 2^(-(1-h_W) m)`", but the
number of band-visits is itself of the order of the number of halvings,
i.e. of `N(X)`: the estimate closes on itself. What is missing is an
**oscillation bound**: control of the number of times a divergent orbit
crosses a level, or of the number of its dyadic band-visits, in terms of
something smaller than `N(X)`. That is exactly the "consecutive windows"
information, and it is not a counting statement.

## 3. What the strip entropy buys (PROVED numerically, exact DP)

Words of length `m` whose partial sums stay in a strip `[-W, 0]` (a
hovering segment in a band of `W` bits) are exponentially rarer than
no-dip words. Probe (A), growth rates between `m = 150` and `300`:

```text
W (bits):   1      2       3       4       6       8       12      16     one-sided
h_W:        --   0.2600  0.6165  0.7520  0.8542  0.8933  0.9244  0.9347   0.94996
1 - h_W:    --   0.7400  0.3835  0.2480  0.1458  0.1067  0.0756  0.0653   0.05004
```

(`W = 1` admits no word of length `>= 3`: the two steps `-1` and `+0.585`
span `1.585 > 1`.) Consequences, all elementary: the number of integers
below `X` whose orbit hovers in a `W`-bit band for `m` steps is at most
`X 2^(-(1-h_W) m + o(m))`, so hovers longer than `(1-h*) L/(1-h_W)`
(`0.20 L` at `W = 4`, `0.47 L` at `W = 8`, `0.77 L` at `W = 16`) contribute
at most `X^(h*)` elements in total over any orbit. This is the one genuine
structural gain of the two-sided constraint; it does not reach hovers
shorter than that, which is where the multiplicity lives.

## 4. The hypothesis (HYP-9161) and what it would give

**HYP-9161 (landing multiplicity is polylogarithmic).** For every `T_b`
orbit, every `X` and every `theta` with `theta log_2 X >= 1`, the sum of
the landing multiplicities is at most `C (log_2 X)^beta` times the number
of landing points, with some `beta < 1` (conjecturally, `C (theta log_2 X)`
suffices, i.e. `beta = 0` at the bootstrap's `theta_X`).

If HYP-9161 holds with exponent `beta`, THM-4499's proof gives
`N(X) <= K X^(h*) (log_2 X)^a` for every `a > beta lambda*/h* - 3/2`; with
`beta = 0`, `a > -3/2`, the ballot floor. The hypothesis is about the
orbit's oscillation, not about words: it fails for no residue class
(section 1) and holds on every actual segment probed (section 1).

## 5. Reassessment of the angle of approach

The one-window method (THM-4476, THM-4487, THM-4495, THM-4498, THM-4499)
is now closed up to the polynomial factor `(log X)^(0.514)`: exponent
`h*`, ballot factor `-3/2`, and the multiplicity `L` are its three
ingredients, and the first two are sharp for the method. The third is a
statement about a single orbit's oscillation; every attempt to prove it by
counting residue classes runs into the same wall as Korec's threshold
(Theorem 4, remark (iii) of the dip-spectrum note): beyond one window the
words are no longer free. So the multiplicity attack is, in substance,
the same problem as lengthening the window, and the honest next angle is
not a sharper count but an *orbit-coupled* quantity: something like the
2-adic-plus-height rank that the procgen reflection calls for, applied to
crossings rather than to descent. I record this rather than press the
count further.

## 6. The owner's snippets, placed

* **`139 = 3^7 - 2^11` and the sporadic cycle.** THM-4484's `-17` cycle has
  shape `(11, 7)`: `2^11 - 3^7 = -139` is its denominator. The pair
  `(2^11, 3^7)` is the atlas's `2187/2048`, and `1093 = (3^7 - 1)/2`; three
  atlas rows are one pair of prime powers.
* **Families with long initial growth and their frequency.** The
  crossroads-family lane's `k + 11` no-descent family has frequency
  decreasing geometrically in `k`; THM-4495 gives the exact rate for the
  full set: the residues with `k` initial no-descent steps have density
  `W_k/2^k = Theta(2^(-(1-h*)k) k^(-3/2))`, and the lane's fractal note
  already uses `W_k` as the exact cylinder count.
* **The `27` comb.** The odd predecessors of `41` are `(41 * 2^(2r+1) - 1)/3`,
  `n -> 4n + 1`, and every third one is divisible by `3` (no odd
  preimage): the lane's section 2 has the exact counts and densities; the
  family shares a tail, not a long initial growth.
* **`{11, 2, 3}` and "a multiple less than the next odd number squared".**
  Not decoded. Tested readings (the `p`-rough multiples `kp` below the
  square of the next odd number, or of the next odd prime; prime powers
  below it; multiples strictly between `p^2` and `(p+2)^2`) do not single
  out `{2, 3, 11}`; the one that fits `2` (multiples `4, 6, 8` below `9`,
  one a square) gives `3` for `3` and `2` for `5, 7, 11, 13`. Left open
  for the owner to restate.

## 7. Zeckendorf three-colourings and the odd-prime decomposition (probe)

The repo's Zeckendorf work (HYP-3739, HYP-4078, HYP-3000; klein's
"shallow diagonal sums of Pascal's triangle") does not contain a
three-colouring along diagonals of the Zeckendorf arrangement, so the
probe defines three natural ones (lowest index mod 3; row plus column of
the (number of terms, top index) arrangement mod 3; number of terms mod 3)
and, for odd numbers, the additive count of odd-prime copies `Omega(n) mod 3`
and the colour `n mod 3` that the map actually sees, and measures the
mutual information with the Collatz-relevant targets (parity, the number
of `2`s factored out of `3n + 1`, first descent within `20` steps, total
stopping time mod `3`):

```text
N = 200000; mutual information I(colour; target) in bits [target entropy in brackets]
                                    parity    v2(3n+1) mod 3 (odd n)               descent<=20       stopping time mod 3
cZ_low                     0.00000 [1.000]           0.00000 [1.379]           0.00001 [0.173]           0.00001 [1.585]
cZ_diag                    0.00001 [1.000]           0.00004 [1.379]           0.00000 [0.173]           0.00001 [1.585]
cZ_len                     0.00001 [1.000]           0.00027 [1.379]           0.00001 [0.173]           0.00000 [1.585]
cP_omega(odd)              0.00000 [-0.000]           0.00000 [1.379]           0.00000 [0.294]           0.00001 [1.585]
n mod 3                    0.00000 [1.000]           0.00000 [1.379]           0.00000 [0.173]           0.00001 [1.585]
baseline pseudo-random colour vs descent<=20: I = 0.000005  (noise level ~ 5.0e-06)
cZ_diag transition matrix under T (rows: colour of n; cols: colour of T(n)); row-normalised:
   colour 0: 0.335  0.328  0.337   (n=55691)
   colour 1: 0.330  0.339  0.330   (n=55507)
   colour 2: 0.332  0.332  0.336   (n=55467)
   marginal colour frequencies: 0.334  0.333  0.333
cZ_diag(n) -> cZ_diag(3n+1) for odd n, row-normalised:
   colour 0: 0.320  0.355  0.325   (n=11005)
   colour 1: 0.341  0.315  0.344   (n=11265)
   colour 2: 0.344  0.329  0.327   (n=11062)
```

Reading (corrected after the audit): every mutual information is below `0.0003` bits, i.e. below `0.02` percent of the target entropy, and the colour transition matrices under `T` and under `n -> 3n+1` are uniform to within `2.2` points. Two of these deviations are statistically significant at `N = 200000` (the number of Zeckendorf terms mod 3 against `v_2(3n+1) mod 3`: `0.00027` bits, about nine times the null expectation, `p = 1.4e-7`; the `n -> 3n+1` matrix, `p = 2e-8`), so the dependence is detectable but negligible in magnitude: there is no usable interlock between the Zeckendorf colourings, the odd-prime count and the Collatz step at this resolution. The only three-colouring the map sees structurally is `n mod 3` (multiples of `3` have no odd preimage), and even that carries no information about descent or stopping time.

The parity algebra the owner describes is the map's own shape: `3n` keeps
`n`'s parity, `+1` flips it, and the halvings remove the even part that
multiplication would otherwise keep; the Collatz step is one additive
parity flip between two multiplicative operations, and the quantity that
decides everything, `v_2(3n+1)`, is the multiplicative side's memory of
the additive flip. Zeckendorf digits are an additive normal form with no
compatibility with `x -> 3x + 1` (the golden base is not a ring), which is
what the probe measures.

## 8. What is provable cleanly: records are ballot-thin (PROVED)

**Proposition.** Let `(y_i)` be a positive one-signed `T_b`-orbit segment
with distinct terms, `X >= 2`, `L = log_2 X`. Call `i` a *future-minimum
record* (leader) if `y_l > y_i` for all `l > i` in the segment, and a
*running-maximum record* (peak) if `y_l < y_i` for all `l < i`. Then the
numbers of leaders and of peaks with `y_i <= X` are each
`O(X^(h*) L^(-3/2))`, with the constant of THM-4499's Lemma M at `y = 1.6`
(plus `2|b| X^(0.585)` small elements).

*Proof.* Leaders: for `y_i >= Y_0 = 2|b| X^(log_2(3/2))` and `1 <= s <= k`,
`y_(i+s) > y_i` and the carry bound `|beta_s| <= |b|(3/2)^s <= y_i/2` give
`M_s > 1/2`, i.e. the `k`-word of `y_i` has all partial sums `> -1`; such
classes number `M_k(1)`, each with at most two representatives below `X`,
and `M_k(1) <= D_s 2^(hk) k^(-3/2) 2^(lambda*) e^(s)` (Lemma M). Peaks: for
a peak `i >= k` with `y_i >= Y_0`, the element `z = y_(i-k)` has, in its
`k`-word, `M_k z + beta_k = y_i > y_(i-k+s) = M_s z + beta_s` for all
`s < k`. Since `|beta_s|, |beta_k| <= |b|(3/2)^k <= Y_0/2 <= y_i/2`,
`M_k z = y_i - beta_k >= y_i/2` and `M_s z = y_(i-k+s) - beta_s < 3 y_i/2`, so
`M_k/M_s > 1/3 > 2^(-1.6)` for every `0 <= s < k`, with no condition on `z`
(a first draft divided by `z`, which is not bounded below; the audit
supplied this form). So `S_k - S_s > -1.585` for all `s < k`: the reversed
word (THM-4495, Step 2) has all partial sums `> -1.6`, and is counted by
`M_k(1.6)`; `z` is an integer `<= X` in one of those classes. Leaders and
peaks with fewer than `k` successors or predecessors in the segment number
at most `k` each. Peaks with
`i < k` number at most `k`. ∎

So the extremal elements of any orbit (its records in either direction)
sit exactly on the ballot floor `X^(h*) (log X)^(-3/2)` of THM-4495,
whatever the orbit does; THM-4499's `(log X)^(0.514)` excess concerns only
the elements that are neither, i.e. the interiors of excursions between
consecutive leaders, and HYP-9161 is the statement that those interiors
are not much longer than the number of their leaders times a polylog. The
same argument bounds the *y-approximate* records (within `y` bits of a
running extremum) by `O(X^(h*) L^(-3/2) 2^(lambda* y) e^(sy))`.

## 9. Independent audit (2026-09-26)

Auditor subagent: `04-computation/experiments/collatz_landing_20260926_audit.py`
-> `05-knowledge/results/collatz_landing_20260926_audit.out` (61 checks, 0
failures; exact partial-sum comparisons). CONFIRMED: section 0's `a*(beta)`
and its arithmetic; the single-step factor; the strip entropies (to four
decimals; the `150 -> 300` rates are finite-size estimates, within `0.003`
of the `300 -> 600` rates); the statements of the Proposition (zero
violations on the orbits of `2^40 - 1`, random 40- and 80-bit starts, `27`,
and the 5n+1 orbit of `7`, at `L = 30, 40, 64`); the discrepancy corollary
(exact identity checks on three orbits). REFUTED and corrected: the first
draft of section 1 (climb-then-drop gives at most `2` dippers per landing
point; hovers give at most about `2L/3`, best found `0.54 L`); the
attribution of probe2's maxima to the initial climbs; the peak proof's
division by `z`; the display "`- 1.038 log_2 log_2 L`" in the corollary (the
coefficient is every value below `1.038`); the phrase "at the noise floor"
for the Zeckendorf probe. ADDED by the audit: the one-bit band lemma
(multiplicity `<= 2 ceil(k/3)` for every orbit). Cosmetic: the synthetic
hover of probe (A) is built in `[-1.49, 0]`, not `[-1, 0]`; for 5n+1 the
Proposition's `Y_0` exceeds `X`, so it asserts nothing there. Verdicts:
Proposition HAS GAPS (repaired above); corollary SOUND (display fixed).
