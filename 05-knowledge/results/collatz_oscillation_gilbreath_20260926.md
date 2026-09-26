# The oscillation lemma is a shell-revisit statement; Gilbreath's difference triangle and the base-6 Collatz automaton compared as cellular automata; why neither conjecture proves the other

**Status: REASSESSMENT with PROVED small statements (the shell-visit
reformulation, the base-6 locality, the Gilbreath persistence mechanism),
FINITE-EXACT probes (orbit segments; the Gilbreath triangle of the primes
below `200000`), SPECULATION marked, and INDEPENDENTLY AUDITED (HAS GAPS ->
repaired, section 7). No Collatz claim, no Gilbreath claim. Session `collatz-oscillation-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/collatz_oscillation_20260926_shells.py`,
`gilbreath_20260926_ca.py`, `collatz_base6_ca_20260926.py`; outputs `.out`
alongside.

## 1. The oscillation lemma as a shell-revisit statement

**Proposition 1 (what the averaged hypothesis says).** In THM-4476's
recursion with depth `D` (dippers fall by `2^D` within `k = floor(log_2 X)`
steps), THM-4506 proves that all dippers of a landing point `j` lie in the
dyadic shell `S_j = (2^D y_j, 2^(D+1) y_j]` with an odd letter between
consecutive ones. Hence the multiplicity `m(j)` is exactly the number of
indices `i` in the `k` steps before `j` with `y_i in S_j` whose first
`2^D`-drop is at `j` (audit: verified on every landing point of the probe;
any two orbit indices inside one dyadic shell are automatically separated
by an odd letter, and a maximal run inside the shell can hold two of
them), and HYP-9161 in its averaged form (`#Dip(X, D) <= C L^mu N(X 2^(-D)) + C L^2`) is equivalent to:

```text
(shell form)  the orbit has, within a window of length k = log_2 X and before leaving a dyadic shell
              downward by 2^D, at most C L^mu indices in that shell on average over the elements below X 2^(-D).
```

Equivalently (regularity form), since every element below `X` is a
dipper, an ND element or one of the last `k`,
`N(X) <= (C L^mu + 1) N(X 2^(-D)) + #ND(X, D) + C L^2 + k`: **the count of orbit
elements cannot jump by more than `C L^mu` between the scales `X 2^(-D)`
and `X`, except through no-dip elements**, which are ballot-thin. With
`D = 1.05 log_2 L`, `X 2^(-D) = X/L^1.05`. *Proof.* Both directions are the
definitions up to constants: `#Dip = sum_j m(j)` over landing points
`j < X 2^(-D)`, and `N(X) = #Dip + #ND + #(last k)`; the converse returns the
averaged form with `C` replaced by `C + 1`. ∎

**Probe (FINITE-EXACT).** On the record orbit of `63728127`, the orbits of
`2^L - 1` and of random odd starts, and `5n+1` orbits, with `L = 24..64`
and `D = ceil(1.05 log_2 L) = 5..7`: the mean multiplicity is `2.3`-`5.8` on
the `3n+1` rows and `1.0`-`4.0` on the `5n+1` rows (growing about like `D`),
the heaviest landing points have `6`-`18` dippers whose window visits to
the shell fall into `4`-`13` maximal runs (the script's run count, which is
not the dipper count: a run can hold two dippers and some visits are cut
off by an earlier drop), and the printed scale ratio (both counts
restricted to full-window indices) is `1.0`-`2.4`; for the record orbit at
`L >= 48` every orbit value lies below `X 2^(-D)`, so its ratio `1.00` is
vacuous. These are orbits that reach `1`: they illustrate the shell
form (small multiplicities, heavy landing points = repeated visits to one
shell, as THM-4506's lemma says) but do not test it on a divergent orbit;
section 6 gives the random-address calibration.

**What would prove it.** A residue class realises any finite pattern of
returns, so the shell form is a property of the orbit as a whole
(THM-4506; crossroads-poset integer note, section 5). The two-place
identity gives no handle: the clock `Q_l = 2^(d_l) m_l/3^l` moves by
`O(1/m)` per step (crossings note, Proposition 2). The needed statement is
that a divergent orbit cannot spend a constant fraction of its time
below `X` oscillating inside single dyadic shells; the strip entropy
(`h_1` does not exist, `h_2 = 0.26`: a one-bit band admits no long stay,
a two-bit band only `2^(0.26 m)` words) makes each such stretch rare among
integers, but an orbit is one integer's address, and its stretches at
successive scales are not independent samples. I record the shell form
as the sharpest statement of the obstacle and do not claim progress
beyond THM-4506's saturation theorem.

## 2. Gilbreath's difference triangle as a cellular automaton (FINITE-EXACT)

Row `0` is the primes; row `r+1` is the absolute differences of row `r`.
Gilbreath's conjecture says every positive-depth row starts with `1`; row `0` starts with `2`. Three facts make
it a cellular-automaton statement (all classical; Odlyzko 1993):

* On the sublattice `{0, 2}` the rule `|a - b|` is `a XOR b` in units of
  `2`: the difference table of a `0/2` row is Pascal's triangle mod `2` in
  the difference orientation (`new[i] = a[i] XOR a[i+1]`, elementary rule
  `102`/`60`, an additive automaton).
* `|1 - a| = 1` for `a in {0, 2}`: once a row is `(1, then only 0s and 2s)`,
  every later row is too, and the conjecture holds from there on.
* An entry `d >= 4` ("defect") at column `p` produces `new[p-1] = d - row[p-1]`
  and `new[p] = d - row[p+1]`: a front moving one cell to the left per row
  and shrinking by `2` whenever the cell to its left is a `2`, plus a
  stationary copy that keeps re-emitting fronts (speeds `1` and `0`). A
  front of size `2j` starting at column `F` meets `F - 1` cells before the
  edge and must keep `2j - 2t >= 4`, so in the random model (each cell a
  `2` with probability `1/2`) it survives with probability
  `sum_(t <= j-2) C(F-1, t) 2^(-(F-1))` (audit; a first draft wrote
  `sum_(i<j) C(F, i) 2^(-F)`, off by one in `j` and `F`), a one-sided binomial
  tail at density `t/F -> 0`; the stationary copy raises the probability
  that the leading `1` is ever lost by a few per cent for `j >= 3`. The
  decay in `F` is exponential with a polynomial factor `F^(j-1)/(j-1)!`.

For the primes below `200000` (`17983` rows): the leading entry is `1` in
every positive-depth row (row `0` starts with `2`); the frontier `F(r)` (first entry `>= 4`) is `3, 8, 25, 59` at rows
`1, 2, 5, 10`, `2763` at row `50`, and from row `65` on the row after its leading `1` is entirely `0/2` (`F = row
length`, maximum entry `2`; rows `1`-`64` each contain an entry `>= 4`, and
none reappears later, by the closure of `{1} x {0, 2}`), so the triangle to the right of its leading `1` is
pure Pascal-mod-2 from row `65` and the conjecture is trivially true for
this range beyond that row. Fresh fronts at the frontier die at once: of
nine fronts in the first `3000` rows (columns `< 160`), eight of size `4`
travel `0, 0, 0, 0, 0, 0, 1, 1` rows and one of size `8` travels `3` (a first
draft sampled each surviving front once per row and reported a mean of
`0.18`). The zero regions of the `0/2` sea: the script's statistic is a
per-row histogram of maximal zero runs (`118983, 59228, 29388, ...` in the
kept block, ratio `2`), which does not distinguish the primes' rows from a
single seed; the audit counted true inverted-triangle sides (a triangle
of side `t` is topped by a maximal run of `t + 1` twos) over the whole sea:
`10031027, 5013216, 2507563, ...` for sides `1, 2, 3, ...` up to `26`, ratio
`2.00`, exactly what a random `0/2` row gives. **So the sides are geometric,
not restricted to `2^k - 1`**: the Sierpinski sides `1, 3, 7, 15, 31, 63` are
the single-seed diagram's; in the primes' triangle `3` and `7` (the cyclic
triangle and the Paley heptagon of the tournament thread) occur, but so
do `2, 4, 5, 6`, each about twice as often as the next. The owner's "zeros
of tournament size edged by 2s" is a feature of the single-seed
automaton, not of Gilbreath's triangle.

## 3. The Collatz map as a cellular automaton (PROVED, brute-force checked)

In base `6` the map `T(x) = x/2` (even) or `(3x + 1)/2` (odd) is a radius-1
cellular automaton on the digit string with one global bit (the parity of
the lowest digit) selecting the branch: the carry of `3 x_i + c` never
crosses a multiple of `6` for `c <= 2`, so the outgoing carry depends on
`x_i` alone, and halving reads one digit up. Brute force over `200000`
random `x` with `11` digits: radius `0` gives `1.3e6` inconsistencies,
radius `1` and `2` give `0`; the odd-branch rule depends on both neighbours
(`216` of `216` neighbourhoods seen); the audit wrote the explicit rule and
checked it exhaustively for all `x < 6^7`. Attribution: the cellular-
automaton reading of `3x+1` goes back to Cloney, Goles and Vichniac ("The
3x+1 problem: a quasi cellular automaton", Complex Systems 1 (1987)
349-360), who work in base two; the base-6 radius-1 form is standard in
later work (the audit cites Korec 1992 and Kari 2012), and a first draft
attributed it to the 1987 paper. So both
conjectures are **edge statements about one-dimensional automata with
arithmetic initial data**: Gilbreath's edge is the leftmost cell of the
difference automaton started from the primes; Collatz's edge is the
lowest base-6 digit (the parity vector) of the digit automaton started
from any integer, and the conjecture is that every finite string reaches
the two-cell cycle `1 <-> 2`.

## 4. Why neither proves the other

* **Encoding and its resource cost (CORRECTED).** Additivity alone does
  not obstruct an edge encoding. Every binary trace is the edge of a
  unique additive-difference seed by the self-inverse Pascal transform
  modulo two. Every nonnegative integer trace is the edge of an absolute-
  difference seed by the positive binomial transform. Applied to an actual
  Collatz orbit, these maps preserve the trace while placing its computation
  in an infinite seed. They do not produce consecutive primes or a
  finite-input reduction. A useful reduction must constrain seed arithmetic
  and resources and preserve the target magnitude or stopping predicate.
  The proofs and actual-source witnesses are in the
  [signed-difference audit](forest_20260926_gilbreath.md). The former phrase
  "which additivity forbids" was false and has been removed.
* **Same combinatorics, different thresholds.** Both heuristics are
  one-sided tails of the binomial distribution, i.e. values of the same
  entropy function `h(p)` at different densities: a Gilbreath front of
  size `2j` survives distance `F` with probability
  `sum_(t <= j-2) C(F-1, t) 2^(-(F-1))` (a lower tail at density `t/F -> 0`,
  where `h -> 0`, so fronts die like `2^(-F)` times a polynomial); a Collatz
  no-descent word of length `k` has probability `W_k/2^k = Theta(2^(-(1-h*)k) k^(-3/2))`
  (an upper tail at density `log_3 2`, where `h* = 0.95`, so bad words die
  only like `2^(-0.05 k)`). THM-4487's curve `E(gamma) = h(max(1/2, gamma/log_2 3))`
  covers only the Collatz range `[h*, 1]`; Gilbreath's tail sits at the
  `p -> 0` end of `h(p)` itself, off that curve (a first draft placed it on
  the curve's `theta -> 0` end, which is the Collatz end `gamma = 1`, an
  error caught by the audit and, independently, by the forest lane). These
  are heuristic comparisons of binomial-tail shapes, the same combinatorics
  at different densities: no map identifying the prime initial law with the
  arithmetic source law, or equating the two conjectures, has been
  established.
* **Rédei (scope corrected by the forest lane).** Gilbreath requires the
  exact leading value `1` at every positive depth. Oddness alone follows
  automatically from the initial `2` followed by odd numbers. Rédei's
  theorem gives an odd number of tournament Hamiltonian paths, not exactly
  one, so a parity analogy cannot supply Gilbreath's missing magnitude
  bound: forest's family `23 + 32s` preserves every row-comparison
  tournament and parity while its final difference grows as `1 + 2s`
  ([forest audit](forest_20260926_gilbreath.md)). My first draft recorded
  the resemblance as thematic; it is weaker than that.
* **What transfers.** Odlyzko's persistence argument (a defect at
  distance `F` needs `F` rows to reach the edge) is exactly the
  carry-bound step of the Terras window (a residue class fixes `k`
  steps): a statement true for `F` more rows is the same as a statement
  true for `k` more steps. Both methods are "one free window", and both
  stop where the window ends: Gilbreath needs the frontier to keep
  growing (an initial-condition statement about the primes), Collatz
  needs the oscillation lemma (an orbit statement). Neither window
  argument reaches the other's obstacle.

## 5. Speculation (marked)

If one wanted a single framework, it is the class of one-dimensional
automata whose edge cell is conjectured to be eventually constant: the
frontier growth of Gilbreath (`F(r) -> row length` by row `100` here) is a
"defect extinction" law, and the Collatz analogue would be a
"shell-revisit extinction" law (section 1): an orbit's stretches of
oscillation inside a shell should die out with scale as Gilbreath's
defects die out with row. In both cases the empirical extinction is
immediate (defects travel `0`-`3` rows; multiplicities are `2`-`6`) and the
missing proof is a statement about the specific initial data (the
primes; one integer's address), not about the automaton. This is the
same conclusion as the crossings note: the one-window method is
complete, and what remains is arithmetic of the seed.

## 6. Calibration in the random-address model: the natural conjecture is `mu = 1/2`, not `0` (FINITE-EXACT)

The script `collatz_oscillation_20260926_randommodel.py` runs the
shell-revisit statistics on the address of no integer: heights in bits
with i.i.d. valuations, `P(v = 1) = p`, `P(v = j) = (1-p) 2^(-(j-1))` for
`j >= 2`, drift `log_2 3 - (3 - 2p)` per odd step (divergent iff `p > 0.7075`),
twenty trajectories of `10^5` steps for each `(p, L)`, `D = ceil(1.05 log_2 L)`:

```text
random-address model: p = P(v=1); drift per odd step = log_2 3 - (3 - 2p) bits; 20 trajectories of 100000 T-steps each
 p     drift   L    D   mean tot<=L  mean ND  mean multiplicity  max multiplicity (over trajectories)  mean landing pts
0.72  +0.025   32   6      5082.3    4018.4     3.01              13                             342.9
0.72  +0.025   64   7      6906.5    4935.6     4.53              21                             441.1
0.72  +0.025  128   8     10456.7    6554.6     6.53              37                             597.6
0.72  +0.025  256   9     14957.7    8639.9     8.76              56                             717.0
0.75  +0.085   32   6       720.5     600.6     2.86              14                              38.8
0.75  +0.085   64   7      1532.3    1186.2     4.63              22                              74.8
0.75  +0.085  128   8      2080.9    1621.7     5.81              28                              76.3
0.75  +0.085  256   9      4970.1    3621.5     8.05              47                             165.2
0.80  +0.185   32   6       216.4     199.5     2.07              10                               5.6
0.80  +0.185   64   7       491.6     443.4     3.49              15                              12.3
0.80  +0.185  128   8       948.2     869.8     4.78              26                              16.4
0.80  +0.185  256   9      1894.6    1754.3     5.20              33                              24.6
0.90  +0.385   32   6       103.9     101.5     1.11               5                               0.9
0.90  +0.385   64   7       196.1     194.7     0.59               6                               0.6
0.90  +0.385  128   8       411.9     407.5     1.03              12                               1.4
0.90  +0.385  256   9       780.7     774.1     1.22               6                               2.4
```

Reading. For the near-critical drifts (`p = 0.72`, `0.75`; the regime of a
hypothetical divergent orbit that is thin, since such an orbit has
`Delta_l -> -infinity` slowly, i.e. vanishing drift) the mean multiplicity
grows like `L^(1/2)` (`3.0, 4.5, 6.5, 8.8` at `L = 32, 64, 128, 256`; doubling
exponents `0.59, 0.53, 0.42`): a zero-drift walk returns to a unit shell
about `sqrt L` times in `L` steps, and each return before the `2^D`-drop is
one more dipper of the same landing point. For a fixed positive drift
`delta` the returns saturate at about `1/delta` (`p = 0.80`: `2.1, 3.5, 4.8, 5.2`;
`p = 0.90`: about `1`). So the random-address heuristic for the shell form is

```text
mean multiplicity  ~  min( sqrt L,  C/delta ),      hence  mu = 1/2 at criticality,
```

and the S7 statement "the random-walk heuristic gives `O(theta L)`" was the
positive-drift case, not the critical one. The actual segments of
section 1 (means `2.4`-`5.8` at `L = 24`-`64`) sit on the `p = 0.72`-`0.75` rows.
Consequences: the honest conjectural exponent of the one-window method
for the thinnest orbits is `a*(1/2) = lambda*/(2h*) - 3/2 = -1.243`, not the
ballot floor `-3/2`; `mu = 0` should be expected only for orbits with
positive drift, whose element count is `O(log X)` anyway. HYP-9161 is
restated accordingly below and in its file. The Gilbreath side has no
such correction: its defects live at threshold `0` where the survival is
`2^(-F)` with no square-root effect (they are killed by the first `2`, not
by a return statistic).

## 7. Independent audit (2026-09-26)

Auditor subagent: `04-computation/experiments/collatz_oscillation_gilbreath_20260926_audit.py`
-> `05-knowledge/results/collatz_oscillation_gilbreath_20260926_audit.out`
(exact integers, numpy, Monte Carlo). CONFIRMED: Proposition 1 with the
precise reading `m(j) = #{shell indices in the window whose first
2^D-drop is j}` (every landing point of the probe), the equivalences up to
constants, all twelve probe cells, the Gilbreath tables (frontier,
run histogram, defect samples), the base-6 radius-1 rule (explicit rule,
exhaustive to `6^7`), the absorbing closure from row `65`. CORRECTED
(applied above): the sentence placing Gilbreath on the `theta -> 0` end of
THM-4487's curve (that end is the Collatz end; the entropy-`0` tail lies at
`p -> 0` of `h(p)`); the attribution of the base-6 automaton to the 1987
paper (base two, "quasi cellular automaton"); the first all-`0/2` row (`65`,
not "about `100`"); the defect statistics (each surviving front was
sampled once per row) and the survival tail (off by one in `j` and `F`,
and the seed's stationary copy re-emits fronts); the zero-triangle
statistic (a per-row zero-run histogram that also shows every size for a
single seed; true triangle sides, topped by runs of twos, are geometric
with ratio `2.00` over the whole sea); the probe's "separate returns"
(maximal runs, not dippers) and its vacuous scale ratio for the record
orbit at `L >= 48`; "Rule 90" for the difference rule (`102`/`60`). Verdict:
HAS GAPS, repaired; the bottom line (no claim on either conjecture, the
duality, no reduction either way) is unaffected.
