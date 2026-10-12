# The first half-band exit cannot be paid at coefficient three: least witness 56,924,955, the threshold 4.64, and source 187

2026-10-12, opus-2026-10-12-S24. **REFUTED, FINITE-EXACT:** the paid
dyadic-band hypothesis (3) of [the three-problems packet, section 3](three_problems_marked_defects_20261012a.md)
at coefficient 3 and at coefficient 4. **PROVED, scoped:** the two all-height
rules for source 187 and the mechanism statement under the random model.
**HEURISTIC + NUMERICAL:** the threshold coefficient `c* = 4.6407`
([HYP-9248](../hypotheses/HYP-9248-first-band-exit-charge-threshold-coefficient.md)).
No universal Collatz claim. Self-audited; no independent peer audit.

## Inheritance and concept board

The incoming [source 127 packet](collatz_paid_exit_127_20261012.md) compiles
all-height exit rules under a linear invoice `3(E(n)-E(h))` and a quadratic
invoice `3(E(n)^2-E(h)^2)`, telescoping to `o(n) <= 3E(n)+16` or
`o(n) <= 3(E(n)^2-1)` on accepted sources, and names 187 as the next
unresolved grammar source. [The three-problems packet](three_problems_marked_defects_20261012a.md)
states the universal target as the first-exit inequality (3) and records that
coefficient 3 survives a finite probe through `10^6` while coefficient 1 fails
at 27. [Charge receipts](collatz_charge_receipts_20261011.md) fix the charge as
the number of odd valuation letters.

Closest proved mechanism: the no-descent rate `3^(h*-1) = 0.9465` per odd
step ([Pascal tower](collatz_pascal_tower_lyapunov_20261005.md); the Beatty
structure of no-descent words in THM-4514). Canonical hostile: the glide
record holders, whose prefixes stay above the source for the longest time.
Corrected near miss: reading a finite probe through `10^6` as support for a
budget whose first failure is a large-deviation event at `5.7 * 10^7`.
Least-used sidecar: the odd-letter fraction of a no-descent word, which the
rate alone does not fix.

Board: first-exit charge / no-descent rate / odd fraction / band budget /
potential increment / own rule. Anchor: the hypothesis (3). Niche: the
six-colour core ([companion](six_colour_core_doubling_20261012.md)).
Wildcard: the two number identities, which the companion reads.

## 1. Coefficient three fails, and the failures are what the rate predicts

Definitions as inherited: `U(x) = (3x+1)/2^a`; the first half-band exit of
odd `n > 1` is the first state `h` with `2(h+1) <= n+1`; `charge` counts the
odd letters of that prefix; `E(n) = floor(log2(n+1))`.

**Least counterexample.** An exhaustive scan of all odd `n < 2^27` (int64,
no orbit left the guard) finds exactly nine sources with
`charge > 3(E(n)^2-E(h)^2)`, the least being

    n = 56,924,955:  198 odd steps to h = 25,168,801,
    E(n) = 25, E(h) = 24, charge 151 > 3(625-576) = 147.

The word is saved in the output. The same scan to `2^33` (the 24 orbits that
left int64 replayed with exact integers, none failing) finds exactly 37
counterexamples and the largest ratio

    charge/(E(n)^2-E(h)^2) = 176/51 = 3.451  at n = 96,883,183.

So the probe through `10^6` was simply below the first large deviation.

**Glide records.** The 34 glide-record holders to `2.6 * 10^18` (Roosendaal's
table, every glide recomputed exactly here) are the natural hostile
population: their prefixes stay above `n` longest. Nine of them fail
coefficient 3, all with drop one (`E(h) = E(n)-1`):

| n | odd steps to the exit | charge | E(n) | budget `3(2E-1)` | ratio |
|---:|---:|---:|---:|---:|---:|
| 56,924,955 | 198 | 151 | 25 | 147 | 3.08 |
| 63,728,127 | 239 | 156 | 25 | 147 | 3.18 |
| 898,696,369,947 | 347 | 243 | 39 | 231 | 3.16 |
| 2,081,751,768,559 | 388 | 272 | 40 | 237 | 3.44 |
| 31,835,572,457,967 | 449 | 321 | 44 | 261 | 3.69 |
| 1,008,932,249,296,231 | 560 | 415 | 49 | 291 | 4.28 |
| 118,303,688,851,791,519 | 569 | 423 | 56 | 333 | 3.81 |
| 180,352,746,940,718,527 | 610 | 438 | 57 | 339 | 3.88 |
| 2,602,714,556,700,227,743 | 634 | 478 | 61 | 363 | 3.95 |

`n = 1,008,932,249,296,231` also refutes coefficient 4 (`415/97 = 4.278`).
No tested source exceeds coefficient 5. Records whose exit has drop two or
more pass easily: the budget then starts at `3(4E-4)`.

The refutation does not touch the compiled finite certificates of the
127 packet (they are exact statements about residue classes), nor the
telescoped total bound on *accepted* sources. It kills the universal task
as posed: a total tester that stops at the first half-band exit with
coefficient 3 rejects 56,924,955.

## 2. Why: the exit charge and the band budget are both linear in E

Random model (letters i.i.d. with `P(a) = 2^-a`): about `N exp(-k KL(nu||P))`
of the odd starts below `N` carry a prefix of length `k` with empirical letter
law `nu` that stays above `n/2`; the constraint is `mean(nu) <= log2 3`. The
longest such prefix has `k = ln N / KL`, and its charge is `f(nu) k` with `f`
the odd-letter mass. Maximizing `f/KL` over admissible laws gives

    c* = (ln 2 / 2) max f/KL = 4.6407,
    optimum f = 0.7420, KL = 0.05541, law (0.638, 0.222, 0.090, 0.031, ...),

against the pure no-descent tilt (`f = 0.7304`, `KL = 0.05498`,
`exp(-KL) = 0.9465 = 3^(h*-1)`, 12.6 odd steps per bit). Thus the largest
first-exit charge below `N` grows like `9.28 E(n)`, while the budget of a
drop-one exit is `3(2E-1) ~ 6E`. Both are linear in `E`; the ratio tends to
`4.64`, and coefficient 3 had to fail. The observed drop-one ratios along the
records rise from 3.08 at `E = 25` to 4.28 at `E = 49`, approaching the model
from below (the longest prefix below `N` is `ln N/KL` only up to the
first-passage correction `-(3/2) ln ln N`), with the record prefixes' odd
fraction 0.72-0.75 matching the tilted law.

Consequence for the proof design: any potential `Phi(E)` that pays first exits
must satisfy `Phi(E)-Phi(E-1) >= 9.28 E`, so `Phi(E) >= 4.64 E^2`. The
quadratic shape of the invoice is forced, its coefficient is forced to be at
least `c*`, and the total it then returns, `o(n) <= 4.64 E^2`, is far weaker
than the expected linear truth. That is the price of paying at the *first*
exit; paying at a deeper exit (`E(h) <= E(n)-s` with `s` chosen by the
source) keeps the telescoping and spends the quadratic budget only where the
large deviation demands it, which is what the 127 compiler's variable drop
already does.

## 3. Source 187 has two paid all-height rules; the linear bank stalls at 235

The actual first-band word of 187 is `(1,2,1,3,1,1,3,4)`, exit 19, charge 6,
`E: 7 -> 4`. The inherited compiler accepts, for every `n = 187 mod 2^16`,
`n >= 187`, with `a = v_2` of the last numerator `>= 4`:

    linear invoice, drop 3:    8(h+1) <= n+1, charge <= 7 <= 9;
    quadratic invoice, drop 1: 2(h+1) <= n+1, charge <= 7 <= 3(49-36).

Both own the source (`first = 187`); the linear invoice rejects drops one and
two (`7 > 3, 6`). With the linear 187 rule the frozen depth-eight grammar
through 8,191 accepts 1,668 sources (first pending 159, whose folded rule is
quadratic); with the quadratic 187 rule, 1,702 (first pending 199).

An auto-compiler that gives each first pending source its own all-height rule
(shortest actual prefix and least drop that the invoice accepts, with the
source inside the compiled height range) separates the two invoices:

| invoice | own rules added | accepted through 8,191 | stops at |
|---|---:|---:|---|
| linear | 159, 187, 199, 219, 231 | 1,802 | 235: no compilable own rule |
| quadratic | 187, 199, 219, 231, 235, 295, 303, 315, 327, 363, 367, 379 | 2,048 | none within twelve rules (pending 411) |

Source 235 has total odd charge 30 against `3E(235) = 21`: its orbit exceeds
the linear budget at every exit (first exit 95 with charge 4 at drop one,
then 91, 103, 61, 5, 1), so under the linear invoice it can only be accepted
through another source's rule and the `+16` seed slack. The quadratic
invoice pays it at the first exit (`4 <= 3(49-36)`).

## 4. Connection ledger

| Source object | Map and preserved predicate | Lost data / sidecar |
|---|---|---|
| Hypothesis (3) at coefficient 3 | exact replay of first exits; counterexample census | nothing lost: explicit words retained |
| Glide records (Roosendaal) | candidate list only; every glide recomputed | record status not needed for the refutation |
| No-descent rate `3^(h*-1)` | `exp(-KL)` at the pure tilt | the odd fraction, supplied by the second multiplier |
| Band budget | `Phi(E) = cE^2`, increment `c(2E-1)` at drop one | the drop distribution at the exit (about half drop one) |
| 127 compiler | unchanged; rules are residue-class certificates | its universal tester at coefficient 3 is dead |

## 5. Validation and reproduction

[Script](../../04-computation/experiments/collatz_band_budget_hostile_20261012.py),
[output](collatz_band_budget_hostile_20261012.json) (default scan to `2^27`)
and [the `2^33` run](collatz_band_budget_hostile_20261012_2p33.json)
(`python3 -B 04-computation/experiments/collatz_band_budget_hostile_20261012.py 8589934593`).
Ordinary and `-O` runs agree byte for byte. The script recomputes every
fetched glide, replays the least witness exactly, checks the ratio record,
re-verifies the 187 rules at 65 members of the class and one member above
`2^90`, rejects the two unpaid linear plans, reproduces the ledger counts of
the 127 packet with the new rule, and performs the tilt optimization with
rational-free bisection (no SciPy). The exhaustive scan uses numba with an
int64 guard; every orbit that leaves the guard is replayed with exact
integers. The glide table was fetched from `www.ericr.nl/wondrous/glidrecs.html`
on 2026-10-12 and is used only as a list of candidates.

Next obligations: extend the census beyond `2^33` with 128-bit arithmetic or
the longer record lists; decide whether coefficient 5 is safe
(HYP-9248); and reformulate the universal tester with a source-chosen drop.
