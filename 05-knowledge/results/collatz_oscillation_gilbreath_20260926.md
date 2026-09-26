# The oscillation lemma is a shell-revisit statement; Gilbreath's difference triangle and the base-6 Collatz automaton compared as cellular automata; why neither conjecture proves the other

**Status: REASSESSMENT with PROVED small statements (the shell-visit
reformulation, the base-6 locality, the Gilbreath persistence mechanism),
FINITE-EXACT probes (orbit segments; the Gilbreath triangle of the primes
below `200000`), and SPECULATION marked. No Collatz claim, no Gilbreath
claim. Session `collatz-oscillation-20260926` (opus), 2026-09-26.**

Scripts: `04-computation/experiments/collatz_oscillation_20260926_shells.py`,
`gilbreath_20260926_ca.py`, `collatz_base6_ca_20260926.py`; outputs `.out`
alongside.

## 1. The oscillation lemma as a shell-revisit statement

**Proposition 1 (what the averaged hypothesis says).** In THM-4476's
recursion with depth `D` (dippers fall by `2^D` within `k = floor(log_2 X)`
steps), THM-4506 proves that all dippers of a landing point `j` lie in the
dyadic shell `S_j = (2^D y_j, 2^(D+1) y_j]` with an odd letter between
consecutive ones. Hence the multiplicity `m(j)` is the number of
odd-separated visits of the orbit to `S_j` in the `k` steps before `j`
that are not cut off by an earlier drop, and HYP-9161 in its averaged
form (`#Dip(X, D) <= C L^mu N(X 2^(-D)) + C L^2`) is equivalent to:

```text
(shell form)  the orbit revisits one dyadic shell, within a window of length k = log_2 X and before
              leaving it downward by 2^D, at most C L^mu times on average over the elements below X 2^(-D).
```

Equivalently (regularity form), since every element below `X` is a
dipper, an ND element or one of the last `k`,
`N(X) <= (C L^mu + 1) N(X 2^(-D)) + #ND(X, D) + C L^2`: **the count of orbit
elements cannot jump by more than `C L^mu` between the scales `X 2^(-D)`
and `X`, except through no-dip elements**, which are ballot-thin. With
`D = 1.05 log_2 L`, `X 2^(-D) = X/L^1.05`. *Proof.* Both directions are the
definitions: `#Dip = sum_j m(j)` over landing points `j < X 2^(-D)`, and
`N(X) = #Dip + #ND + #(last k)`. ∎

**Probe (FINITE-EXACT).** On the record orbit of `63728127`, the orbits of
`2^L - 1` and of random odd starts, and `5n+1` orbits, with `L = 24..64`
and `D = ceil(1.05 log_2 L) = 5..7`: the mean multiplicity is `2.4`-`5.8`
(growing about like `D`), the heaviest landing points have `6`-`18`
dippers made of `4`-`13` separate returns to one shell, and the scale
ratio `N(X)/N(X 2^(-D))` is `1.0`-`2.4`. So on actual segments the shell
form holds with `mu = 0` and a small constant, and the heaviest landing
points are exactly repeated returns to a shell (oscillation), as
THM-4506's lemma says.

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
Gilbreath's conjecture says every row starts with `1`. Three facts make
it a cellular-automaton statement (all classical; Odlyzko 1993):

* On the sublattice `{0, 2}` the rule `|a - b|` is `a XOR b` in units of
  `2`: the difference table of a `0/2` row is Pascal's triangle mod `2`
  (an additive automaton).
* `|1 - a| = 1` for `a in {0, 2}`: once a row is `(1, then only 0s and 2s)`,
  every later row is too, and the conjecture holds from there on.
* An entry `>= 4` ("defect") moves one cell to the left per row and
  shrinks by `2` whenever the cell to its left is a `2`; it can reach the
  edge only by travelling its whole distance `F` through the `0/2` sea
  meeting fewer than `size/2` twos: the probability in the random model is
  `sum_(i < size/2) C(F, i) 2^(-F)`, a one-sided binomial tail whose
  threshold `i/F` tends to `0`.

For the primes below `200000` (`17983` rows): the leading entry is `1` in
every row; the frontier `F(r)` (first entry `>= 4`) is `3, 8, 25, 59` at rows
`1, 2, 5, 10`, `2763` at row `50`, and from about row `100` the entire row is
`0/2` (`F = row length`, maximum entry `2`), so the triangle is pure
Pascal-mod-2 below row `100` and the conjecture is trivially true for this
range beyond that row. Defects at the frontier die at once: size `4`
defects travel on average `0.18` rows (they meet a `2` immediately), the
largest observed travel is `3` rows (size `8`). The maximal inverted zero
triangles in the `0/2` region have sizes `1, 2, 3, ..., 15` with counts
halving at each step (`118983, 59228, 29388, ...`): a zero triangle of side
`t` is a run of `t + 1` equal entries in a row, which a random `0/2` row has
with probability `2^(-t)`. **So the sizes are geometric, not restricted to
`2^k - 1`**: the Sierpinski sizes `1, 3, 7, 15` arise only from a single
seed, not from the primes' rows; `3` and `7` (the cyclic triangle and the
Paley heptagon of the tournament thread) occur, but so do `2, 4, 5, 6`,
each about twice as often as the next. The owner's "zeros of tournament
size edged by 2s" is therefore a feature of the single-seed automaton,
not of Gilbreath's triangle.

## 3. The Collatz map as a cellular automaton (PROVED, brute-force checked)

In base `6` the map `T(x) = x/2` (even) or `(3x + 1)/2` (odd) is a radius-1
cellular automaton on the digit string with one global bit (the parity of
the lowest digit) selecting the branch (Cloney, Goles, Vichniac 1987): the
carry of `3 x_i + c` never crosses a multiple of `6` for `c <= 2`, so the
outgoing carry depends on `x_i` alone, and halving reads one digit up.
Brute force over `200000` random `x` with `11` digits: radius `0` gives
`1.3e6` inconsistencies, radius `1` and `2` give `0`; the odd-branch rule
depends on both neighbours (`216` of `216` neighbourhoods seen). So both
conjectures are **edge statements about one-dimensional automata with
arithmetic initial data**: Gilbreath's edge is the leftmost cell of the
difference automaton started from the primes; Collatz's edge is the
lowest base-6 digit (the parity vector) of the digit automaton started
from any integer, and the conjecture is that every finite string reaches
the two-cell cycle `1 <-> 2`.

## 4. Why neither proves the other

* **Duality of difficulty.** Gilbreath's automaton is additive on the
  bulk (`XOR`), with all the difficulty in the initial condition (the
  primes' gaps supply the defects) and a trivial absorbing rule at the
  edge. The Collatz automaton has a trivial initial condition (any
  integer) and all the difficulty in the rule, which is not additive in
  any base (the carry couples digits and the branch bit couples the whole
  string). A reduction of Collatz to Gilbreath would have to encode a
  non-additive automaton's global behaviour into an additive one's edge,
  which additivity forbids; a reduction of Gilbreath to Collatz would
  have to produce the primes' difference rows from a Collatz orbit, and
  no orbit quantity in this thread (parity words, valuations, carries
  `S_l`, the two-place clock) has prime-gap structure.
* **Same combinatorics, different thresholds.** Both heuristics are
  one-sided binomial tails: a Gilbreath defect survives distance `F` with
  probability `sum_(i < j) C(F, i) 2^(-F)` (threshold `j/F -> 0`, entropy
  `h -> 0`, so defects die exponentially fast in `F`: `2^(-F)`); a Collatz
  no-descent word of length `k` has probability `W_k/2^k = Theta(2^(-(1-h*)k) k^(-3/2))`
  (threshold `log_3 2`, entropy `h* = 0.95`, so bad words die only like
  `2^(-0.05 k)`). Gilbreath is the `theta -> 0` end of the entropy curve
  of THM-4487, Collatz the `gamma = 1` end. That is the precise sense in
  which they are the same problem at different temperatures, and why one
  is heuristically easy and the other hard.
* **Rédei.** The "leading `1`" of Gilbreath and Rédei's odd number of
  Hamiltonian paths are both parity statements, and Pascal's triangle
  mod `2` is the additive automaton behind both the `0/2` sea and the
  parity of binomial coefficients; but no map from difference rows to
  tournaments is in sight, and I record the resemblance as thematic.
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
