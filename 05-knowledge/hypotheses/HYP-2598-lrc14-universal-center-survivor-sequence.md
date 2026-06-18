# HYP-2598: LRC(14) Universal-Center Survivor Sequence

**Status:** proved bounded-spread skeleton / residual classifier.  This does
not prove LRC(14), but it turns the integer-sequence lane into a clean proof
split.

Script: `04-computation/lrc14_integer_sequence_carriers_codex.py`  
Output: `05-knowledge/results/lrc14_integer_sequence_carriers_codex.out`

HYP-2597 proves that the pure cluster good set has universal centers

`0, 1/2, 1/3, 2/3`,

because at a rational `a/b` every integer offset lands on the `1/b` grid, so a
gap of size at least `1/b` exists, and `1/b > 2/7` exactly for `b < 7/2`.

This hypothesis asks when those universal centers survive the small-speed safe
set `G_P`.

## Survivor Lemma

Let `P subset {1,...,13}` and let `E` be a cluster-offset shape with
`max(E)=R`.

For `R >= 4`:

- If every speed in `P` is odd, then the center `1/2` contributes a safe
  reservoir of measure at least `3/(14R)`.
- If no speed in `P` is divisible by `3`, then the centers `1/3` and `2/3`
  contribute a safe reservoir of total measure at least `2/(21R)`.
- If both conditions hold, these add to `13/(42R)`.

The proof is just the denominator-cap window plus small-speed margin.  Around
`1/2`, HYP-2597/HYP-2590 gives a cluster-good window of half-width
`3/(28R)`.  Odd small speeds sit at distance `1/2` from the origin, so for
`R>=4` the entire conservative window remains level-`1/14` safe.  Around
`1/3` and `2/3`, the conservative half-width is `1/(42R)`; small speeds not
divisible by `3` sit at distance `1/3`, and for `R>=2` the whole window remains
safe.  The `R>=4` statement covers both at once.

## Integer Sequence

For `s=|P|`, the number of small parts whose universal-center skeleton survives
is

`survivor(s) = binom(7,s) + binom(9,s) - binom(5,s)`.

Reason:

- there are `7` odd speeds in `{1,...,13}`;
- there are `9` speeds not divisible by `3`;
- their intersection has size `5`.

Thus

`survivor(s) = 1, 11, 47, 109, 156, 146, 91, 37, 9, 1, 0, 0, 0, 0`

for `s=0..13`.  The mixed residual is the complement

`0, 2, 31, 177, 559, 1141, 1625, 1679, 1278, 714, 286, 78, 13, 1`.

Equivalently, the universal skeleton dies asymptotically exactly when `P`
contains at least one even speed and at least one multiple of `3`.

## LRC(14) Readout

Writing `k=|E|` and `s=13-k`, the center-surviving / mixed split is:

| `k` | center-surviving `P` | mixed residual `P` | total `P` |
|---:|---:|---:|---:|
| 3 | 0 | 286 | 286 |
| 4 | 1 | 714 | 715 |
| 5 | 9 | 1278 | 1287 |
| 6 | 37 | 1679 | 1716 |
| 7 | 91 | 1625 | 1716 |
| 8 | 146 | 1141 | 1287 |
| 9 | 156 | 559 | 715 |
| 10 | 109 | 177 | 286 |
| 11 | 47 | 31 | 78 |
| 12 | 11 | 2 | 13 |
| 13 | 1 | 0 | 1 |

This is a useful simplification: the universal denominator-cap skeleton proves
a `c/R` continuous reservoir for all-odd or 3-free small parts.  The true
Diophantine tail is now isolated as the mixed parity/triadic residual.

Named hard rows classify exactly as expected:

- the HYP-2595 constant-100 failure `P=(1,2,11)` is 3-free, so it has
  coefficient `2/21`;
- `quarter_min`, `near_via_min`, `via_zero_k7`, and the old `84m` hard-core row
  are mixed residuals with coefficient `0`;
- `P=empty` has coefficient `13/42`.

## Conservative-Window Stabilization

The exact computation intersects the conservative universal windows with
`G_P`.  As `R` grows, the positive skeleton counts converge to the survivor
sequence above.  At `R=160`, selected rows match exactly:

- `s=1`: `11` positive, `11` survivors;
- `s=2`: `47` positive, `47` survivors;
- `s=3`: `109` positive, `109` survivors;
- `s=4`: `156` positive, `156` survivors;
- `s=5`: `146` positive, `146` survivors;
- `s=8`: `9` positive, `9` survivors;
- `s=9`: `1` positive, `1` survivor.

The transient mixed-positive rows at smaller `R` are bounded-spread accidents:
their universal windows are wide enough to leak around small-speed danger
teeth, but they vanish as the spread grows.  This confirms the conceptual
split: fixed universal centers solve a bounded-spread skeleton; recurring
non-universal intervals are needed for the large-spread floor.

## Proof Target

The remaining LRC(14) S3 route should now be organized as:

1. Center-surviving small parts: use the explicit `c/R` reservoir plus a finite
   denominator-placement lemma.
2. Mixed residual small parts: prove colored resonance / spread recurrence
   prevents all non-universal recurring intervals from being blocked.
3. Low finite endpoint: classify q-grid obstruction families and exact-check
   only structurally reduced cores, as in HYP-2596.

See also HYP-2597, HYP-2596, HYP-2595, HYP-2593, HYP-2590, and OPEN-Q-108.
