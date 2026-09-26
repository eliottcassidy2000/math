# Gilbreath's sea as the Frobenius tower of the Fermat numbers; the five Fermat primes and the five Platonic solids; the wall theorem for size-4 defects

**Session:** opus, `gilbreath-fermat-platonic-20260926`, 2026-09-26.
**Owner's directive:** continue the S9 work (Gilbreath's difference triangle as an
automaton, the oscillation lemma), "especially the connection regarding
Gilbreath's conjecture and the 5 Fermat primes, and the 5 Platonic solids".
**Inherits:** S9 note
[`collatz_oscillation_gilbreath_20260926.md`](collatz_oscillation_gilbreath_20260926.md)
(the `0/2` sea is the XOR automaton; fronts move left at speed 1; the first
all-`0/2` row of the primes below `200000` is row `65`); forest's
[board](forest_20260926_board.md) section 3 (the von Staudt-Clausen bridge
`denom(B_(2^r)) = 2 prod_(2^j <= r, F_j prime) F_j`, the constructible-polygon
criterion) and [Gilbreath audit](forest_20260926_gilbreath.md) (every trace has
an additive seed; parity is weaker than the exact leading `1`);
[HYP-3772](../hypotheses/HYP-3772-platonic-tilings-schlafli-curvature-genus-covering-min.md)
and [HYP-3771](../hypotheses/HYP-3771-five-tilings-five-solids-apex-prime-geometry.md)
(the Schlafli count `(p-2)(q-2) < 4` and the `(2,3,n)` angle-defect spine);
[THM-871](../../01-canon/theorems/THM-871-fermat-rung-rigidity.md) (Fermat-rung
rigidity of rotational tournaments: at a Fermat prime `p`, `Aut = Z_p`).

**Status: PROVED small statements (the sea kernel theorem; the Mersenne-side
and `3^(K-1-m)` count of the single-seed zero triangles; the skew-Hadamard
doubling tower and the containment `F_21 <= Aut(T_k)`; THM-4511, the wall
theorem for a lone size-4 defect, with the exact extinction law `2^(1-F)`);
FINITE-EXACT (all checks on the primes below `200000`, the tower tournaments to
order `255`, the group identifications, the extinction table to `F = 12` over
all `2^(F-1+10)` contexts); CLASSICAL facts cited as such; SPECULATION marked.
No claim on Gilbreath's conjecture, none on Collatz.** Scripts:
`04-computation/experiments/gilbreath_fermat_tower_20260926.py`,
`..._tournaments.py`, `gilbreath_fermat_platonic_20260926_groups.py`,
`gilbreath_extinction_20260926.py`, with `.out` files.

## 0. The answer in one paragraph

The `0/2` sea of Gilbreath's triangle evolves by `new[i] = a[i] XOR a[i+1]`, and
`t` rows down this is the Frobenius power `(1 + x)^t` over `F_2`: the cell is
the parity of the sea bits over the Lucas support `{j : j subset t}`, whose
generating polynomial, evaluated at `x = 2`, is `prod_(i in bits t) F_i`, a
product of Fermat numbers `F_i = 2^(2^i) + 1`. So the **Fermat numbers are the
dyadic coarse-graining kernels of Gilbreath's sea**, exactly, for every `t`.
Their *primality* never enters: the automaton uses `1 + x^(2^i)` whether or not
`2^(2^i) + 1` is prime. The five known Fermat primes enter three other exact
places: Gauss-Wantzel (the first `32` rows of the single-seed Sierpinski
diagram, read in binary, are the `32` products of distinct known Fermat primes,
i.e. the odd orders of constructible polygons; row `32` is the composite `F_5`),
forest's Bernoulli denominators, and THM-871's Fermat rungs. The five Platonic
solids come from a different mechanism (the Schlafli inequality
`(p-2)(q-2) < 4`), and touch the Fermat primes only through the set
`{3, 4, 5} = {F_0, 2^2, F_1}`: these are the Schlafli numbers, the orders of
the polygons Euclid needs in Book XIII, *and the sizes of the fields whose
projective groups are the rotation groups of the solids*
(`PSL(2,3) = A_4`, `PGL(2,3) = S_4`, `PSL(2,5) = PSL(2,4) = A_5`). The two fives
are not the same five. What the automaton picture does deliver is a theorem
about Gilbreath's actual dynamics: **a lone defect of size `4` never crosses
the first sea `2` on its left, whatever lies to its right** (THM-4511), so its
extinction probability is exactly `2^(1-F)`, and the owner's "zeros of
tournament size edged by 2s" is made exact: in the single-seed diagram every
zero triangle has Mersenne side `2^m - 1`, each such order carries a doubly
regular tournament built by the same doubling, and the symmetry those
tournaments keep is precisely the Paley heptagon's Frobenius group of order
`21`.

## 1. The sea kernel theorem: Fermat numbers as coarse-graining kernels (PROVED, FINITE-EXACT)

Write the sea entries as bits `b = a/2`. The rule `|a - b|` on `{0, 2}` is
`2 (b_i XOR b_(i+1))`, i.e. multiplication by `1 + x` in `F_2[[x]]` (difference
orientation). Hence

```text
(sea kernel)   b_(r+t)(i) = XOR_(j subset t) b_r(i + j)      whenever a_r(i..i+t) is all in {0, 2},
```

because `(1 + x)^t = sum_j C(t, j) x^j` and `C(t, j)` is odd iff `j subset t`
(Lucas). The window condition is the light cone: the rule looks only to the
right, so the cell `(r + t, i)` depends only on `a_r(i), ..., a_r(i + t)`. The
kernel of step `t` is the row `t` of Pascal's triangle mod `2`, and read as a
binary number it is

```text
sum_(j subset t) 2^j  =  prod_(i in bits(t)) (1 + 2^(2^i))  =  prod_(i in bits(t)) F_i .
```

Checked: for `n < 128`, `C(n, j)` odd iff `j subset n`, and row value
`= prod F_i` (script part 1). On the actual triangle of the primes below
`200000`: rows `1`-`64` each contain an entry `>= 4`, row `65` is `1` followed
by `0/2`, and rows `65 + t` equal the kernel image of row `65` for
`t = 1..64` and `t = 127, 255, 511, 1023, 2047` (69 exact checks, leading `1`
throughout). Inside the sea left of the frontier of row `10` (`F(10) = 59`,
defect `4`), all `1653` in-cone cells over `t = 1..57` equal the kernel image of
the sea word. Frontier values `F(r) = 3, 8, 25, 59, 291, 870, 2770, 2763, 5942,
5940` at `r = 1, 2, 5, 10, 20, 30, 40, 50, 60, 64`.

So at step `t = 2^k` the sea is `b(i) XOR b(i + 2^k)` (the kernel `F_k`), and
at step `2^k - 1` it is the full window parity over `2^k` cells (the kernel
`prod_(i<k) F_i = 2^(2^k) - 1`). **The primality of `F_k` is invisible to the
automaton.**

## 2. Where the five Fermat primes do live (CLASSICAL + FINITE-EXACT)

* **Rows of the single-seed diagram.** Row `2^k` has value `F_k`; rows
  `2^k - 1` are all ones with value `2^(2^k) - 1 = prod_(i<k) F_i`
  (`3, 15, 255, 65535, 4294967295`). The `32` rows `0..31`, read in binary, are
  exactly the `32` products of distinct known Fermat primes
  (`1, 3, 5, 15, 17, 51, 85, 255, 257, ...`), i.e. `1` and the `31` odd orders
  of constructible regular polygons (Gauss-Wantzel); row `32` is
  `F_5 = 4294967297 = 641 * 6700417`, the first non-constructible row, and every
  row `n` with `32 <= n < 2^33` carries a factor `F_i`, `5 <= i <= 32`, all
  known composite. Checked exactly (script part 2). Forest's bridge is the
  same list read through von Staudt-Clausen: `denom(B_(2^r)) = 2 prod F_j` over
  the Fermat primes with `2^j <= r`.
* **Zero triangles of the single-seed diagram.** With the seed placed with
  room on both sides, every maximal run of `2^m` ones in a top row `n < 2^K`
  is bounded by zeros and tops an inverted zero triangle of side `2^m - 1`
  whose left edge is vertical and whose right edge is the diagonal of ones,
  and there are no other zero triangles; the number of side-`(2^m - 1)`
  triangles with top row `n < 2^K` is `sum_(n: exactly m trailing ones)
  2^(popcount(n) - m) = 3^(K-1-m)` for `m < K`, plus the single triangle of
  side `2^K - 1` under the all-ones row `2^K - 1`. Checked for `K = 5, 7, 9`:
  counts `27, 9, 3, 1, 1`; `243, 81, 27, 9, 3, 1, 1`; `2187, ..., 3, 1, 1`.
  The sides `1, 3, 7, 15, 31` in the first `32` rows are the Sierpinski sides
  of the five known Fermat primes' rows (the side-`(2^k - 1)` triangle opens
  under row `2^k`, whose binary value is `F_k`); every side `2^m - 1` with
  `m >= 2` is `3 mod 4`. In the primes' triangle (S9, audited) the sides are
  geometric with ratio `2.00`, all sizes: **the owner's tournament-size
  triangles are the seed's, not the primes'** (S9, unchanged).
* **THM-871.** At a Fermat prime `p` the multiplier group `Z_(p-1)` is a
  `2`-group, so every rotational tournament on `p` vertices has `Aut = Z_p`
  (klein). That is the tournament thread's own use of the same five primes.
* **The `F_k n + 1` maps.** Drift per odd step `log_2 F_k - 2`:
  `3: -0.415`, `5: +0.322`, `17: +2.09`, `257: +6.01`, `65537: +14.0`.
  Collatz is the only Fermat-prime Conway map with negative drift; the `5n+1`
  map of the S9 probes is `F_1`.
* **A non-result.** The constructible sequence `s_n = prod_(i in bits n) F_i`
  has no Gilbreath-type regularity: the leading column of its absolute
  difference triangle is `2, 0, 8, 8, 16, 0, 80, 64, 32, 32, 896, ...`.

## 3. The doubling tower as tournaments: Mersenne sides carry doubly regular tournaments (PROVED + FINITE-EXACT)

The Sierpinski kernel is the Kronecker power of `[[1,0],[1,1]]`; the same
doubling with signs gives skew Hadamard matrices:

```text
H_2 = [[1, 1], [-1, 1]],     H_(2n) = [[H, H], [-H^T, H^T]].
```

*Proof that `H_(2^k)` is skew Hadamard.* If `H = S + I` with `S^T = -S` then
`H_(2n) - I = [[S, S + I], [S - I, -S]]`, whose transpose is its negative; and
`H_(2n) H_(2n)^T = [[2 H H^T, 0], [0, 2 H^T H]] = 2n I`. Normalization (first
row `+1`, first column `-1` off the diagonal) is inherited. Deleting row and
column `0` gives a tournament `T_k` on `2^k - 1` vertices (`i -> j` iff
`H[i][j] = +1`), doubly regular by the classical skew-Hadamard correspondence
(every pair has `(n - 3)/4` common out-neighbours). Checked for `k <= 8`
(`n = 3, 7, 15, 31, 63, 127, 255`). The recursion is explicit and checked:

```text
T_(k+1) = T_k + {0'} + T_k'   (T_k' = arc-reversed copy),   i -> i',   i -> j' iff i -> j,
          i' -> j iff i -> j,   0' -> every i,   every i' -> 0'.
```

Every automorphism `s` of `T_k` extends to `T_(k+1)` as `(s, s)` fixing `0'`
(all arcs above are defined through `T_k`'s arcs and the copy map), so
`Aut(T_k) <= Aut(T_(k+1))`. Computed (individualization-refinement
backtracking, verified on `P_7` and `P_31`): `|Aut(T_2)| = 3` (cyclic
triangle), `|Aut(T_3)| = 21` and `T_3` is isomorphic to the Paley heptagon
`P_7`, and `|Aut(T_k)| = 21` for `k = 4, 5, 6` (`n = 15, 31, 63`; the run for
`k = 7, 8` is in the `.out` file). For `T_5`: element orders `{3: 14, 7: 6}`,
i.e. the Frobenius group `Z_7 x| Z_3` of the heptagon, with orbits of sizes
`7, 1, 7, 1, 7, 1, 7`. `T_5` is **not** isomorphic to the Paley tournament
`P_31` (`|Aut(P_31)| = 465`) although its 4-vertex census
(`4340` with a source, `13020` strong, `4340` with a sink, `9765` transitive)
equals `P_31`'s (the census did not separate them; the automorphism
groups do). So: the zero-triangle sides `2^k - 1` of the sea are exactly the orders
of a tower of doubly regular tournaments generated by the sea's own doubling,
the Paley heptagon is its `k = 3` level, and from there on the tower keeps the
heptagon's `21` symmetries and nothing else ([HYP-9162](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md)).
This is the exact content behind the owner's "zeros of some tournament size
edged by 2s": true of the seed diagram (sides Mersenne, tournaments doubly
regular), false of the primes' rows (sides geometric).

## 4. The five Platonic solids and the fields of size 3, 4, 5 (CLASSICAL, VERIFIED)

* `(p - 2)(q - 2) < 4` has exactly five solutions `{3,3}, {3,4}, {4,3}, {3,5},
  {5,3}` and `= 4` exactly three (`{3,6}, {4,4}, {6,3}`): HYP-3772/HYP-3771,
  re-checked. So the Schlafli numbers are `{3, 4, 5} = {F_0, 2^2, F_1}`, and
  `(17 - 2)(q - 2) >= 15` means the first Fermat prime without a solid is `17`.
* The rotation groups are the projective groups over the fields of those
  sizes, verified by building each group as a permutation group of the
  projective line `P^1(F_q)`: `PSL(2,3)` has order `12` on `4` points, all
  even, hence `A_4` (tetrahedron; the `4` points are its vertices);
  `PGL(2,3)` has order `24` on `4` points, hence `S_4` (cube/octahedron; the
  `4` points are the cube's diagonals); `SL(2,4)` has order `60` on the `5`
  points of `P^1(F_4)`, all even, hence `A_5`; `PSL(2,5)` has order `60` on
  `6` points (the six five-fold axes of the icosahedron), has `15` involutions
  and exactly `5` Klein four-subgroups, and its conjugation action on them is
  faithful with image of order `60` inside `A_5`, hence `PSL(2,5) = A_5` (the
  five objects are the five cubes in the dodecahedron). Also `PSL(2,2) = S_3`
  (the triangle), `PSL(2,7)` of order `168` (Klein quartic, HYP-3771's `n = 7`
  frontier) and `PSL(2,17)` of order `2448`, not a finite rotation group of
  `3`-space (Klein's list: cyclic, dihedral, `12`, `24`, `60`).
* Constructibility: the icosahedron `(0, +-1, +-phi)` and the dodecahedron
  `(+-1, +-1, +-1), (0, +-1/phi, +-phi)` are regular (checked: `12` resp. `20`
  vertices, equal edges, `5` resp. `3` nearest neighbours) with coordinates
  in `Q(sqrt 5)`; Euclid XIII's constructions use exactly the constructibility
  of the `3`-, `4`- and `5`-gon, i.e. the Fermat primes `3` and `5`. The
  binary polyhedral groups (orders `24, 48, 120`) sit in the unit quaternions,
  THM-871's rung `5 = F_1` of the Cayley-Dickson table; the icosians of
  THM-870 are the `120` case.
* **Verdict on "5 = 5".** The Platonic five counts integer solutions of a
  curvature inequality; the Fermat five counts the `k` for which the
  evaluation `2^(2^k) + 1` of the tower's kernel `1 + x^(2^k)` is prime, and
  whether it stays five is open. They share the set `{3, 5}` and nothing
  else: no map sends solids to Fermat primes. What is exact is the triangle
  `{3, 4, 5}` = Schlafli numbers = constructible polygons up to `5` = field
  sizes of the rotation groups, and `3 = F_0` is also the Collatz multiplier
  (the only Fermat multiplier with negative drift, section 2).

## 5. THM-4511: the wall theorem for a lone size-4 defect, and the exact extinction law (PROVED + FINITE-EXACT)

**Theorem.** Run `a_(r+1)(i) = |a_r(i) - a_r(i+1)|` from a row with
`a_0(0) = 1`, `a_0(F) = 4`, and `a_0(i) in {0, 2}` for every other `i`
(finitely or infinitely many columns to the right). Let `c < F` be the largest
column with `a_0(c) = 2` and `1 <= c`, if it exists. Then

1. if `c` exists, every column `i <= c` stays in `{0, 2}` for all rows and the
   leading entry stays `1`: **no entry `>= 4` ever crosses the wall `c`**,
   whatever the cells right of `F` are;
2. if `c` does not exist (`a_0(1..F-1) = 0`), then `a_(F-1)(1) = 4` and
   `a_F(0) = 3`: the leading `1` is destroyed at row `F`.

Hence a lone size-`4` defect at distance `F` destroys the leading `1` iff the
`F - 1` cells between them are all `0`; in a uniform random sea the
probability is exactly `2^(1-F)`, and re-emissions from the stationary copy
never add anything.

*Proof.* Column `i` depends only on columns `i` and `i + 1` (the rule looks
right). Columns `> F` start in `{0, 2}` and stay there. Column `F` starts at
`4`; while its right neighbour is `0` it stays `4`, the first time the
neighbour is `2` it becomes `|4 - 2| = 2`, and afterwards it stays in `{0, 2}`.
Call a column sequence *good* if it takes values in `{0, 4}` up to some row,
then possibly the value `2`, then values in `{0, 2}` for ever (or stays in
`{0, 4}` for ever). Column `F` is good. If column `i + 1` is good and column
`i` starts at `0` (`c < i < F`), column `i` is good: it stays in `{0, 4}` while
column `i + 1` does (`|x - y|` for `x, y in {0, 4}`), it becomes `|x - 2| = 2`
the row after column `i + 1` first shows a `2`, and afterwards both are in
`{0, 2}`. So columns `F - 1, ..., c + 1` are good. Column `c` starts at `2`:
`|2 - x| = 2` for `x in {0, 4}`, so it stays `2` while column `c + 1` is in
`{0, 4}`; it becomes `0` when column `c + 1` shows `2`; afterwards both are in
`{0, 2}`. So column `c` never leaves `{0, 2}`, hence (induction leftward)
neither do columns `< c`, and `|1 - a(1)| = 1`. For 2: on `{0, 4}` the rule is
XOR in units of `4`, so the cone of the defect over the zeros is Pascal mod
`2` and its left edge `(s, F - s)` is `4` for every `s`; the cone
`[F - s, F]` never reaches column `F + 1`, so the right side is irrelevant.
∎

*Verification.* Exhaustive enumeration over all `2^(F-1+R)` sea contexts
(`F - 1` cells left, `R = 10` right, tableau run for `F + R` rows, which is
exact for every loss caused while the copy is provably alive) for `F = 3..12`:
the minimum column ever holding a `4` equals `F - z` (`z` = trailing zeros of
the left sea) in every configuration, and `p_4(F) = 2^(1-F)` to all printed
digits; Monte Carlo with `4 * 10^6` samples: `p_4(14) = 1.267e-4` vs
`1.221e-4`, `p_4(16) = 2.875e-5` vs `3.052e-5` (Poisson noise).
Canon file: [THM-4511](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md).

*What the theorem says about Gilbreath.* Rows `1`-`64` of the primes below
`200000` contain `29` fresh fronts (rows where the frontier is not the previous
frontier minus one), `23` of them of size `4`, at distances `3, 8, 14, 14, 25,
59, 98, ...` (list in the `.out`). A size-`4` front is harmless unless its
whole zero-run reaches column `1`; the riskiest moment of the whole triangle
was row `1` (the `4` at distance `3`, with the sea `2, 2` before it: survival
probability `1/4` in the random model), then row `2` (distance `8`,
`1/128`). The front-only risk sum over fresh fronts is `0.258`, over all rows
`0.258`, with exact re-emission corrections `0.258`; Gilbreath's conjecture
survived its first two rows with random-model probability about `0.74` and
has faced only negligible risk since. **The theorem does not touch larger
sizes**: for a lone `6` the wall breaks to a `4` (`|6 - 2| = 4`) and the exact
extinction probability exceeds the front-only tail `F 2^(1-F)` by
`6`-`10` per cent (`F = 8`: `0.067993` vs `0.062500`; `F = 12`: `0.006338` vs
`0.005859`); for a lone `8` by `4`-`8` per cent (`F = 8`: `0.244003` vs
`0.226563`), the excess being the stationary copy's re-emitted fronts, which
for `d >= 6` can pass a wall the first front has already weakened. Nor does it
cover several defects: a second `4` arriving at column `F` after the first has
died can reopen the wall (column `c` may then be `0`), which is the one-window
limitation of every persistence argument (S9).

## 6. Bearing on the two conjectures (honest)

* **Gilbreath.** The Fermat tower is the structure of the *safe* region; the
  conjecture is a statement about the frontier, i.e. about defect extinction.
  THM-4511 settles extinction for lone size-`4` defects exactly (the wall);
  sizes `>= 6` and interacting defects are governed by the same left-looking
  mechanism but no closed law is proved. The five Fermat primes and the five
  solids do not enter the conjecture.
* **Collatz.** `3 = F_0` is the multiplier; the base-`6` automaton (S9) is
  nonlinear; the `F_k n + 1` maps for `k >= 1` have positive drift. No
  reduction in either direction (S9 and forest, unchanged).
* **Tournaments.** The Mersenne orders `2^k - 1` of the sea's zero triangles
  carry the doubling tower `T_k`, doubly regular with `Aut = F_21` as far as
  computed; the Paley heptagon is its third level; THM-871's Fermat-rung
  rigidity is the other face of the same five primes.

## 7. Speculation, marked

1. **An extinction theorem for size `2j`.** The proof of THM-4511 is a
   monotonicity of column sequences (`{0, 4}` then `2` then `{0, 2}`). For a
   lone `2j` the column values form a descending chain of levels
   `2j, 2j - 2, ...` with the same left-looking structure; a level-by-level
   version of goodness might give the exact law `p_(2j)(F)` in closed form,
   reproducing the `6`-`10` per cent excess computed here. This is the next
   exact target if the Gilbreath thread continues.
2. **Fermat and Collatz through `3 = F_0`.** The only structural fact found is
   the drift sign; the S8/S9 material on `{2, 3, 11}` and the base-`6` locality
   is unrelated to Fermat primality.
3. **The owner's fives.** If a sixth Fermat prime existed, Gauss-Wantzel would
   add `32` constructible rows and nothing would change in the sea, in the
   tower, or in the solids. The fives are a coincidence of counts; the exact
   common object is the set `{3, 4, 5}`.
