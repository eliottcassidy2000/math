# Insufficient precision is the no-descent residual: the classical residue sieve, the one-member theorem (THM-4512), why "every unresolved orbit enters a certified region" is the conjecture itself, and the Langlands and loom readings

**Session:** opus, `gilbreath6-collatz-precision-20260926`, 2026-09-26.
**Owner's directive:** "Collatz remains open. The next target is controlling
the insufficient-precision cases or proving that every unresolved orbit
enters a certified region. Think about past work regarding the Langlands
program and the loom for connections with those broad themes."
**Inherits:** the reset/swaplift lane
([board](reset_20260926_board.md), [swaplift](reset_20260926_swaplift.md):
hard swap `r = 1`, recreated precision `R`, the finite certificate bank of
`65` residue classes of density `0.1893` among all integers, target 2
"characterize or charge the insufficient-precision residual"); Terras 1976
and Lagarias 1985 (coefficient stopping time; the stopping-time classes);
[THM-4495](../../01-canon/theorems/THM-4495-collatz-no-descent-exact-order-spitzer.md)
(exact order of the no-descent counts),
[THM-4499](../../01-canon/theorems/THM-4499-collatz-thin-divergence-little-o.md),
[THM-4506](../../01-canon/theorems/THM-4506-collatz-landing-multiplicity-shell.md);
the S8 crossings note (the two-place identity, the earlier Langlands remark).

**Status: PROVED (THM-4512, elementary) + FINITE-EXACT (residual densities,
thresholds, stopping-time comparison to `10^7`) + REASSESSMENT of the target
+ SPECULATION marked. Collatz OPEN; nothing here changes that.** Scripts:
`04-computation/experiments/collatz_precision_residual_20260926.py`,
`collatz_coefficient_stopping_20260926.py`, with `.out` files.

## 0. The answer in one paragraph

"Precision" is the number of 2-adic bits of `n` that the orbit has consumed:
after `j` Syracuse steps with valuations `v_1, ..., v_j` it is `A = v_1 + ... +
v_j`, and the odd integers sharing that word form one residue class modulo
`2^A`. A class is *certified* when its word has a coefficient descent
(`3^j < 2^A`), and then, by THM-4512, every member above an explicit
threshold descends, and the threshold excludes at most the class
representative (none at all for `j <= 14` except `n = 1`, and no odd
`3 <= n <= 10^7` has actual and coefficient stopping times differing). So the
"insufficient-precision residual" is exactly the set of *no-descent words*,
whose density among odd integers after `k` Syracuse steps is `D(k)`:
`D(41) = 0.00265` (the classical sieve certifies `99.73` per cent of odd
integers within `41` steps, against the swaplift bank's `37.87` per cent),
`D(60) = 0.00062`; in `T`-coding the residual counts are THM-4495's `W_k`,
of exact order `2^(h* k) k^(-3/2)`, and per Syracuse step the residual
decays with exponent `(1 - h*)/rho* = 0.0793` (`D(k)^(1/k) -> 2^(-0.0793) =
0.9465`), since a no-descent word of `k` odd steps has about `k/rho*` map
steps at the critical density. Controlling this residual *in density*
is therefore done; controlling it *pointwise* is the conjecture: an orbit is
unresolved at precision `A` iff its `A`-bit word has no descent, "resolving"
it means extending the word until a descent appears, and "every unresolved
orbit enters a certified region" is, verbatim, "every `n` has finite
stopping time", which is Terras's form of the conjecture. The density form
of that statement is trivial (`2^(-L)`); the pointwise form is the problem.

## 1. Dictionary (exact)

| lane / owner vocabulary | classical object | in this repository |
|---|---|---|
| precision `R`, recreated precision at a swap | 2-adic bits consumed, `A = v_1 + ... + v_j` | Terras bijection `n mod 2^A <-> word` (THM-4476 setting) |
| hard swap `r = 1` | valuation `1` step (`n = 3 mod 4`) | the letter `v = 1` |
| certified class (descent within `J + 1` steps for all admissible `v`) | residue class mod `2^A` whose word has a coefficient descent at some `j <= J + 1` | THM-4512: certified up to at most its representative |
| insufficient-precision residual | words without coefficient descent up to the available `A` | the no-descent set: `W_k` (`T`-coding, THM-4495), `D(k)` (Syracuse coding, below) |
| "every unresolved orbit enters a certified region" | every `n` has finite stopping time | the conjecture (Terras) |

The swaplift bank's `65` classes descend within `J + 1 <= 41` Syracuse steps,
and an actual descent forces a coefficient descent (THM-4512, part 1), so
the bank is contained in the classical certified set `C_41` of precision at
most `65` bits. Its density `0.1893` among all integers (`0.3787` among odd)
should be read against `1 - D(41) = 0.9973` among odd integers.

## 2. The residual density (FINITE-EXACT)

`D(k)` = sum over valuation words of length `k` with no coefficient descent at
any `j <= k` of `2^(-A_k)` = the density, among odd integers, of classes not
certified within `k` Syracuse steps (dynamic programme over `(j, A)`, exact
rationals):

```text
 k :  1     2     3     4      5       6     7         8          12        16        20        24
D(k): 1/2   3/8   1/4   13/64  19/128  1/8   113/1024  367/4096   0.05212   0.03092   0.01925   0.01314
 k :  28       32       36       40       41       48       52       56       60
D(k): 0.00875  0.00593  0.00427  0.00299  0.00265  0.00156  0.00112  0.00082  0.00062
```

`D(k)^(1/k)` rises from `0.865` at `k = 40` to `0.884` at `k = 60`; the limit is
`2^(-(1-h*)/rho*) = 2^(-0.0793) = 0.9465` (a no-descent word with `k` odd
steps has `A ~ k/rho` map steps and weight `2^(-A)`, and `(1 - h(rho))/rho` is
minimal at `rho = rho* = log_3 2`, where it equals `0.0793`); the polynomial
factor of THM-4495 (`A^(-3/2)` in `T`-coding) is what makes the approach
slow, and its exact form in Syracuse coding is not derived here. In `T`-coding the corresponding counts are
`W_k` (`W_41 = 12805670000`, `W_60 = 2216134944775156`, density `W_k/2^k`
among all integers). The two codings index precision differently (`k`
Syracuse steps versus `k` map steps) and are not to be compared row by row.

## 3. THM-4512: certified classes have at most one uncertified member (PROVED)

`U^j(n) = (3^j n + S_j)/2^A` with `S_j = sum_(t<j) 3^(j-1-t) 2^(A_t) > 0`.
Actual descent at `j` forces `3^j < 2^A`; conversely `3^j < 2^A` gives
`U^j(n) < n` for every `n > N(w) = S_j/(2^A - 3^j)`. Since `A_t <= A - (j - t)`,
`S_j <= 2^A ((3/2)^j - 1)`, so `N(w)/2^A <= ((3/2)^j - 1)/(2^A - 3^j)`, which
is `< 1` at the minimal admissible `A = ceil(j log_2 3)` for every
`j <= 5000` (worst `0.507` at `(j, A) = (5, 8)`; beyond `5000` any effective
irrationality measure for `log_2 3` does it). Hence a coefficient-descent
class contains at most one member not certified by its coefficient descent,
its representative `rho_w = -S_j 3^(-j) mod 2^A`, and only when
`rho_w <= N(w)`. Enumerating all `606746` first-descent classes with
`j <= 14` (valuations to `40` beyond critical; beyond that
`S_j/2^A <= 2^(-v) 2((3/2)^j - 1) < 1`): the only uncertified representative
is `rho = 1` (word `(2)`, `U(1) = 1`). Directly, `sigma(n) = sigma_inf(n)` for
every odd `3 <= n <= 10^7` (maximal `155`). Canon:
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md).
Terras's conjecture (`sigma = sigma_inf` for all `n >= 2`) stays open beyond
these ranges; the theorem reduces it to representatives.

The largest thresholds sit at the near-convergent pairs `(j, A)` of
`log_2 3`: `N(w)/2^A = 0.25` at `(1, 2)`, `0.144` at `(3, 5)` (word
`(1, 2, 2)`), `0.096` at `(5, 8)` (word `(1, 2, 1, 2, 2)`), then below `0.05`.
This is the precise form of the lane's "coupling between the core's division
budget and the precision supplied by the orbit": the coupling is
`2^A - 3^j` against `S_j`, and it costs at most one member per class.

## 4. The target, honestly

* **Density form.** Odd `n` whose first `L` Syracuse iterates all avoid the
  `k`-step certified region must have `v = 1` at each of them (a valuation
  `>= 2` is a one-step descent), so the density is at most `2^(-L)` for every
  `k >= 1`, exactly `2^(-L)` for `k = 1`. "Orbits enter certified regions
  quickly, in density" is trivial.
* **Pointwise form.** `n` is unresolved at precision `A` iff its `A`-bit word
  is a no-descent word; it becomes resolved exactly when the word acquires a
  descent, i.e. at its stopping time; "every unresolved orbit eventually
  enters a certified region" is "every `n >= 2` has finite stopping time",
  Terras's form of the conjecture (the strong induction from stopping times
  to reaching `1` is classical). No intermediate statement of the form
  "unresolved orbits are controlled" is available beyond the density
  results, which the repository already holds in exact order: the
  no-descent counts (THM-4495), the polynomial orders of the dip spectrum
  (THM-4498), and the thin-divergence count `o(X^(h*))` (THM-4499). The
  shell-revisit lemma (HYP-9161, S9/S10) is the one pointwise-flavoured
  statement in this thread, and it remains open.
* **What the lane's board can use.** Its target 2 asks to characterise the
  insufficient-precision residual: it is the no-descent set, with the
  exact order above; its members are the words whose partial valuation
  sums stay below `j log_2 3`, the width-two poset / ballot structure of
  THM-4503 and THM-4495. Its target 1 (extend the bank by structural
  families) is the classical residue sieve, which already covers `99.7`
  per cent at `41` steps; an all-height certificate per class is exactly the
  coefficient descent, complete up to the representative (THM-4512).

## 5. Langlands and the loom (SPECULATION, marked)

The owner asked for the broad-theme connections. The honest map:

* **Two sides.** The 2-adic side of Collatz is completely understood:
  Terras's bijection `n mod 2^k <-> parity word`, Lagarias's conjugacy of
  `T` on `Z_2` to the shift, Haar measure, ergodicity, the densities of
  section 2. The integer side is where the conjecture lives: `Z^+` is a
  measure-zero subset of `Z_2`, and the conjecture says its points behave
  like Haar-generic points. In Langlands-program language the closest
  honest analogue is not functoriality but the *local-global* gap: local
  (2-adic and 3-adic) information is complete and global integrality is the
  obstruction, as with the difference between analytic density statements
  (Dirichlet, `L`-functions) and pointwise ones (Linnik-type) for primes.
  The 3-adic thread is the repository's two-place identity
  `2^(d_l) m_l = 3^l n + S_(l-1)`, `m_l = 2^(-d_l) S_(l-1) mod 3^l` (S7/S8): the
  carry `S` is the only place where the two threads meet.
* **The loom.** If the loom is the weaving of the 2-adic and 3-adic threads,
  the weft is the carry `S_l`, and the repository's exact results about it
  (the carry bound of THM-4476, the reset lane's carry-defect potential)
  are the only global statements. The braids/seams notes of September
  (`arithmetic_braids2_*`, `arithmetic_seams_*`) are about floor-reciprocity
  and Fano-code identities; they do not touch the Collatz carry.
* **Gilbreath.** Its sea is the `F_2`-Frobenius tower (S10), and the Fermat
  primes enter through cyclotomic 2-extensions (Gauss), i.e. abelian class
  field theory, the `GL_1` corner of Langlands; the conjecture itself has no
  automorphic side that anyone has found. The tournament `L`-function idea of
  the concept map (item 34) and the S8 remark (`11` as the conductor of
  `X_0(11)`, the theta series) are unrelated to either conjecture.
* **Verdict.** No reduction, no bridge; the useful content of the analogy is
  the discipline it imposes: every density statement is a statement on the
  analytic side, and only a global (integral) mechanism, not more precision,
  can move the pointwise problem.

## 6. Next

The pointwise target has one candidate lever in this thread: the
shell-revisit lemma (HYP-9161 shell form) with `mu = 1/2` as the calibrated
target. Everything about precision and certified regions is settled in
density and reduces to the conjecture pointwise.
