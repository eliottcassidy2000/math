# The Fibonacci line with its zero removed: three consecutive 1s, Zeckendorf and negaFibonacci as two copies, the real Binet function; and the Collatz reading: growth families are 2-adic shadows of negative rational cycles

**Session:** opus, `fibonacci-two-copies-20260926`, 2026-09-26.
**Owner's directive:** "a new smooth nice function that is 1 at 1, -1 and 0
and allows us to set up our Zeckendorf tricolor decomposition along with
another copy of it such that there is only 1 extra 1, in between the
positive and negative numbers which each use each other's copies of 1 as
their 3rd 1 for their coloring. Consider how this relates with Binet's
formula and what information it captures more precisely, and synthesize
new related ideas that help progress our Collatz proof attempts." The
owner also summarised the reset lane's state: recognizable classes with
arbitrarily long growth certified by three parameterized rules, a separate
closed family through `27` with a decreasing recursion parameter, and the
gap: a finite certificate checker may reject, and visiting a locally
decreasing value does not repay the original source.
**Inherits:** the S8 Wythoff/Zeckendorf tricolor
([`collatz_crossings_20260926_potential_and_seeds.md`](collatz_crossings_20260926_potential_and_seeds.md)),
[THM-4507](../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)
(rational shadows), the reset lane's families
([swaplift](reset_20260926_swaplift.md) sections 3, 6),
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
and the S11 precision note, Knuth's negaFibonacci representation, Lagarias
1985 (periodic parity vectors and rational cycles).

**Status: PROVED elementary identities (the real Binet function's twist and
Cassini identities, the interpolant family; the shadow proposition) +
FINITE-EXACT (negaFibonacci uniqueness and class densities to `10^5`; the
three integer negative cycles as the only integer fixed points among the
`120` primitive no-descent words with `p <= 8`; the two families' 2-adic
limits) + SPECULATION marked (the analogy). Collatz OPEN.** Scripts:
`04-computation/experiments/fibonacci_two_copies_20260926.py`,
`collatz_cycle_shadows_20260926.py`, with `.out` files.

## 0. The answer in one paragraph

Delete the `0` from the two-sided Fibonacci sequence and the three units
`F_(-1) = F_1 = F_2 = 1` become consecutive: `..., -8, 5, -3, 2, -1, 1, 1, 1,
2, 3, 5, 8, ...`. The positive Zeckendorf system lives on the right of the
middle `1` (it uses `F_2, F_3, ...`), Knuth's negaFibonacci system on the
left (it uses `F_(-1) = 1, F_(-2) = -1, F_(-3) = 2, ...` and represents *every*
integer of either sign), and the middle `1 = F_1` is the extra unit used by
neither. The smooth function behind this is the real Binet function
`F(x) = (phi^x - cos(pi x) phi^(-x))/sqrt 5`: it passes through all Fibonacci
numbers of both signs, satisfies the recurrence for real `x`, and its
negative copy is the positive one with the sign twist,
`F(-x) = -cos(pi x) F(x) + phi^(-x) sin^2(pi x)/sqrt 5`; the information the
smooth function adds to the integers is exactly this phase `cos(pi x)`
(and nothing else that the integers can see: every smooth interpolant
obeying the recurrence differs from Binet's by `r phi^(-x) sin(pi x)`, and
Binet's is the one with minimal oscillation). No such function is `1` at
`0`: `F(0) = F(2) - F(1) = 0` is forced, so "1 at `-1`, `0`, `1`" is a property
of the zero-removed *sequence*, not of a recurrence-respecting function.
The tricolor extends: the negaFibonacci classes by lowest index have the
same densities for positive and for negative integers, `phi^(-2), phi^(-3),
phi^(-3), phi^(-4)` (four classes: lowest term `+1`, `-1`, other odd, other
even), and the two "unit classes" on the positives (Zeckendorf ending in
`F_2`, negaFibonacci ending in `F_(-1)`) are independent in density
(`phi^(-2) * phi^(-2) = phi^(-4) = 0.1459`, observed `0.1459`). The Collatz
reading that actually does work is the *other copy of the map*: the
negative integers, whose three odd cycles `-1`, `-5`, `-17` are exactly the
three parameterized growth rules, and more generally every no-descent
valuation word `w` has a negative rational cycle point `x_w` whose 2-adic
neighbourhoods are the recognizable growth classes (shadow proposition,
section 4); the `27` family converges 2-adically to `-5`, the
`(4^k + 17)/3` family to `17/3`, a positive point whose orbit reaches `1`,
which is why its large members descend at step `6`. Precision is 2-adic
distance to a cycle point, and "not repaying the source" is the factor
`(3^p/2^A)^m` of the shadowed periods.

## 1. The zero-removed line and the two systems (CLASSICAL, checked)

`g(n) = F_(n+1)` for `n >= 0`, `g(-n) = F_(-n)` for `n >= 1`:
`g(-6..6) = -8, 5, -3, 2, -1, 1, 1, 1, 2, 3, 5, 8, 13`, with `g(-n) = (-1)^(n+1)
g(n - 1)`. Zeckendorf: every positive integer is uniquely a sum of
non-consecutive `F_k`, `k >= 2`. NegaFibonacci (Knuth): every integer, of
either sign, is uniquely a sum of non-consecutive `F_(-k)`, `k >= 1`
(`F_(-k) = (-1)^(k+1) F_k`); checked by enumeration for `|N| <= 10^5`
(existence and uniqueness). Small cases: `1 = F_(-1)`, `-1 = F_(-2)`,
`2 = F_(-3)`, `-2 = F_(-1) + F_(-4)`, `3 = F_(-1) + F_(-3)`, `-3 = F_(-4)`,
`4 = F_(-2) + F_(-5)`, `-4 = F_(-2) + F_(-4)`, `5 = F_(-5)`, `8 = F_(-1) + F_(-3) +
F_(-5)`. So the two copies share the units as the owner says: the positive
system's unit is `F_2`, the negative system's units are `F_(-1) = +1` and
`F_(-2) = -1`, and `F_1` sits between, used by neither.

## 2. The smooth function and what it captures (PROVED identities)

`F(x) = (phi^x - cos(pi x) phi^(-x))/sqrt 5` is the standard real Binet
extension (`psi^x = (-1)^x phi^(-x)` read as `cos(pi x) phi^(-x)`). Checked
to `10^-12` on a grid, and proved by the same two-line computation each:

* `F(n) = F_n` for all integers `n` of both signs; `F(x + 2) = F(x + 1) + F(x)`
  for all real `x`.
* Twist: `F(-x) = -cos(pi x) F(x) + phi^(-x) sin^2(pi x)/sqrt 5`. At integers
  the correction vanishes and this is `F_(-n) = (-1)^(n+1) F_n`.
* Smooth Cassini: `F(x + 1)^2 - F(x) F(x + 2) = cos(pi x)` for all real `x`
  (using `phi^2 + phi^(-2) = 3`).
* All interpolants: a function satisfying the recurrence for real `x` and
  taking the values `F_n` at all integers is `F(x) + r(x) phi^(-x) sin(pi x)`
  with `r` `1`-periodic; for constant `r`, `F_r(x) = F(x) + r phi^(-x)
  sin(pi x)/sqrt 5` satisfies the recurrence (checked) and oscillates on the
  negative axis with amplitude `phi^|x| sqrt(1 + r^2)/sqrt 5`, minimal at
  `r = 0`: **Binet's formula is the unique minimal-oscillation smooth
  extension, and the integers cannot see `r`.** That is the precise sense
  in which the smooth function captures more: the phase `cos(pi x)` (the
  sign pattern of the negative copy) is forced, the quadrature component
  `sin(pi x)` is free and invisible.
* Zeros on the negative axis: `0, -0.1838, -1.5708, -2.4704, -3.5109,
  -4.4958, -5.5016, -6.4994, -7.5002, ...`, converging to the half-integers
  `-(k + 1/2)` (where `cos(pi x) = 0`): the negative copy is an oscillation
  of growing amplitude with zeros at the half-integers.
* `F(0) = 0` for every interpolant, since `F(0) = F(2) - F(1)`. A function
  that is `1` at `-1`, `0` and `1` therefore cannot respect the recurrence;
  the zero-removed sequence `g` is a pair of Binet branches, `F(x + 1)` on
  `x >= 0` and `F(x)` on `x <= -1`, glued across `(-1, 0)`. If a smooth marker
  of the three points is wanted, `w(x) = exp(-x^2 (x^2 - 1)^2)` is one (`1`
  exactly at `-1, 0, 1`, positive, rapidly decaying), but it carries no
  Fibonacci structure; I do not claim it is the owner's function.

## 3. The tricolor on both copies (FINITE-EXACT, `|N| <= 10^5`)

Zeckendorf classes on positives (lowest index `2` / odd `>= 3` / even `>= 4`):
`0.3820, 0.3820, 0.2361` = `phi^(-2), phi^(-2), phi^(-3)` (the S8 Wythoff
classes). NegaFibonacci classes by lowest index `k` (`1`: the term `+1`; `2`:
the term `-1`; odd `>= 3`; even `>= 4`):

```text
positive N:  k=1: 0.3820   k=2: 0.2361   odd>=3: 0.2361   even>=4: 0.1459
negative N:  k=1: 0.3820   k=2: 0.2361   odd>=3: 0.2360   even>=4: 0.1459
```

i.e. `phi^(-2), phi^(-3), phi^(-3), phi^(-4)`, identical for the two signs:
the negative copy is coloured like the positive one. On the positives the
two unit classes, "Zeckendorf ends in `F_2`" and "negaFibonacci ends in
`F_(-1)`", have density `phi^(-2)` each and intersect in density `0.1459 =
phi^(-4) = phi^(-2) phi^(-2)`: independent in density. (A first guess that
the negaFibonacci index set is the Zeckendorf set shifted by one is false:
it holds for no `N <= 10^5`.) These are exact-looking observations at
`10^5`; I have not proved the densities for the negative system.

## 4. The Collatz reading: growth families are 2-adic shadows (PROVED + FINITE-EXACT)

The "other copy" of the Collatz map is the map on negative integers (the
`3x - 1` map on positives). Its odd cycles with `|n| <= 10^6` are exactly
three: through `-1` (valuation word `(1)`, `A = 1`, factor `3/2`), through
`-5` (word `(1, 2)`, `A = 3`, factor `9/8`), through `-17` (word
`(1, 1, 1, 2, 1, 1, 4)`, `A = 11`, factor `2187/2048`; the owner's `139 = 3^7 -
2^11` is `3^p - 2^A` here). The general mechanism (Lagarias 1985: periodic
parity vectors are rational cycles; THM-4507's rational shadows) in exact
form:

**Proposition (shadow rule).** Let `w = (v_1, ..., v_p)` be a valuation word,
`A = v_1 + ... + v_p`, `S_w = sum_(t<p) 3^(p-1-t) 2^(A_t)`. The odd integers
whose first `p` valuations are `w` form one residue class modulo `2^(A+1)`
(the last valuation must be exact, which costs the extra bit). The affine
map of the word is `n -> (3^p n + S_w)/2^A`, with fixed point
`x_w = S_w/(2^A - 3^p)`, the unique 2-adic integer with word `w^infinity`. If
`n = x_w mod 2^(mA+1)` then `n` follows `w^m` and

```text
U^(mp)(n) = (3^p/2^A)^m (n - x_w) + x_w .
```

`x_w < 0` iff `3^p > 2^A` iff `w` is a growth word. *Proof.* The first
statement is Terras's correspondence between `n mod 2^(A+1)` and the parity
word of length `A + 1`; the rest is the affine composition and the fixed
point. ∎ Checked: the three integer cycles for `m = 1, 3, 10` on random
members of the classes (`n = 3 mod 4`, `15 mod 16`, `2047 mod 2^11`;
`n = 11 mod 16`, `1019 mod 2^10`, `2147483643 mod 2^31`; `n = 4079 mod 2^12`,
...), and every one of the `120` primitive no-descent words with `p <= 8`
(up to rotation) with `m = 2`.

Consequences for the lane's board:

1. **The three parameterized rules are the three integer cycles.** Among the
   `120` primitive no-descent words with `p <= 8`, exactly three have an
   integer fixed point: `(1) -> -1`, `(1, 2) -> -5`, `(1, 1, 1, 2, 1, 1, 4) ->
   -17`. Every other no-descent word has a negative *rational* cycle point
   and gives a recognizable growth class just the same (for example `(1, 1,
   2)` gives `x = -19/11`, factor `27/16` per three steps, classes `n = -19/11
   mod 2^(4m+1)`), so there are infinitely many "rules"; the integer ones are
   the three whose shadowed point belongs to the family itself.
2. **The `27` family.** `n_(k+1) = 8 n_k + 35` (`27, 251, 2043, 16379, ...`)
   has fixed point `-5`, and `n_k = -5 mod 2^(3k)`: it is the `-5` shadow,
   growing by `9/8` per period for `k` periods (the lane's `U^(2j)(n_k) = 4 *
   9^j 8^(k-j) - 5` is the shadow formula with `x = -5`).
3. **The `(4^k + 17)/3` family** (`7, 11, 27, 91, 347, 1371, ...`,
   `n_(k+1) = 4 n_k - 17`) converges 2-adically to `17/3`, a *positive*
   rational whose orbit under `U` is `17/3 -> 9 -> 7 -> 11 -> 17 -> 13 -> 5 ->
   1`; its members shadow a descending point, and indeed their first-descent
   times are `37, 28` (for `27, 91`) and then `6, 6, 6, 6, 6` for `k >= 5`. The
   "decreasing recursion parameter" family is a descent shadow, not a
   growth one.
4. **Precision is 2-adic distance.** A positive `n` shadows the cycle point
   `x` for exactly `floor((v_2(n - x) - 1)/A)` periods; the certificate
   length of a growth rule is the number of bits `n` shares with `x`. "A
   finite certificate checker may reject" is the statement that `n` is
   still inside the shadow at the checker's precision.
5. **Not repaying the source.** After `m` shadowed periods the orbit sits at
   `(3^p/2^A)^m (n - x) + x`; the debt is the factor, and it is repaid only
   if the continuation word has a coefficient descent *relative to `n`*
   (THM-4512), never by a locally smaller value. This is the exact form of
   the owner's sentence.
6. **Coverage.** At precision `K` bits the three integer rules cover density
   `3 * 2^(-K)` of the integers, the union of all rational-cycle shadows at
   that precision is the set of classes whose `K`-bit word is periodic, and
   the whole no-descent residual is `D(k) ~ 10^(-2)`-`10^(-3)` (S11):
   aperiodic no-descent words dominate the residual by an exponential
   margin. Growth *rules* explain what the residual's periodic sliver looks
   like, not the residual.

## 5. The analogy (SPECULATION, marked)

Fibonacci's two copies are joined by the sign twist `cos(pi x)`; Collatz's
two copies are joined by the 2-adic completion: the negative integers are
2-adic limits of positive growth families (`-5 = lim 27, 251, 2043, ...`), and
the twist is the sign of `3^p - 2^A`. The three consecutive units
`F_(-1), F_1, F_2` have a natural Collatz counterpart in the three trivial
cycle points `-1, 0, 1` of `T` (`0 -> 0`; `1 -> 2 -> 1`; `-1 -> -2 -> -1`): a
symmetric core of three units around `0`, with the two copies on either
side and the middle one (`0`, the fixed point) used by neither sign. The
three colours could be read as the three negative integer cycles, the
three integer growth rules. This is a pattern match: base `phi` numeration
and the bases `2`, `3` of Collatz share the "two copies with a twist" shape
and nothing mechanical, and I claim no more.

## 6. Next

The shadow proposition turns the lane's "growth rules" into a statement
about the periodic sliver of the residual; the residual itself is
aperiodic, and the pointwise question ("does every `n` leave every shadow
and every aperiodic no-descent stretch") is the conjecture in Terras's
form (S11). The one open exact question raised here that is not the
conjecture: prove the negaFibonacci class densities `phi^(-2), phi^(-3),
phi^(-3), phi^(-4)` and the independence of the two unit classes.
