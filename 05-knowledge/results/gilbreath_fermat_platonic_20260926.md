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
and `3^(K-1-m)` count of the single-seed zero triangles, proof in section 2; the skew-Hadamard
doubling tower and the containment `F_21 <= Aut(T_k)`; THM-4511, the wall
theorem for any number of size-4 defects, with the exact extinction law `2^(1-F)` for a lone one);
FINITE-EXACT (all checks on the primes below `200000`, the tower tournaments to
order `255` including their automorphism groups, the group identifications, the extinction table to `F = 12` over all `2^(F-1+R)` contexts, `R = 10`, `R = 9` at
`F = 12`); CLASSICAL facts cited as such; SPECULATION marked; INDEPENDENTLY
AUDITED (HAS GAPS -> repaired, section 8). No claim on Gilbreath's conjecture,
none on Collatz.** Scripts:
`04-computation/experiments/gilbreath_fermat_tower_20260926.py`,
`..._tournaments.py`, `gilbreath_fermat_platonic_20260926_groups.py`,
`gilbreath_extinction_20260926.py`, `gilbreath_certificate_20260926.py`, with `.out` files.

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
are not the same five. What the automaton picture does deliver is a theorem about Gilbreath's actual
dynamics: **in a row whose entries are `0`, `2` and `4` only, no `4` ever
crosses the first sea `2` on the left of the first `4`, whatever `0/2/4`
pattern lies to the right, and the leading `1` is destroyed iff the first `4`
is preceded by zeros only** (THM-4511; the audit removed my restriction to a
single defect, which the proof never used). So a lone size-`4` defect's
extinction probability is exactly `2^(1-F)`, a finite Gilbreath triangle is
decided at its first all-`{0,2,4}` row, and all danger to the leading `1`
comes from entries `>= 6`. And the owner's "zeros of tournament size edged by
2s" is made exact: in the single-seed diagram every zero triangle has
Mersenne side `2^m - 1`, each such order carries a doubly regular tournament
built by the same doubling, and the symmetry those tournaments keep, as far
as computed (orders up to `511`), is the Paley heptagon's Frobenius group of
order `21` (HYP-9162).

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

## 2. Where the five Fermat primes do live (CLASSICAL + PROVED + FINITE-EXACT)

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
  counts `27, 9, 3, 1, 1`; `243, 81, 27, 9, 3, 1, 1`; `2187, ..., 3, 1, 1`; the
  audit also checked that every interior zero lies in exactly one such
  triangle and the closed form to `K = 16`. *Proof.* A maximal zero run of
  length `t` in row `r + 1` comes from `t + 1` equal cells of row `r` bounded
  by unequal neighbours, and those cells are ones unless the run continues a
  zero run above; so the zero region is exactly the union of the inverted
  triangles under the maximal runs of ones. In row `n` the ones sit at
  `p - j` for `j subset n` (Lucas); if `n` has exactly `m` trailing ones, the
  sets `{H, H + 1, ..., H + 2^m - 1}`, `H` ranging over the subsets of the
  higher bits of `n`, are the maximal intervals of such `j` (adding `1` to
  `H + 2^m - 1` carries into the zero bit `m` of `n`; `H - 1` has bit `m` set),
  so every run has length `2^m`, every side is `2^m - 1`, and the count over
  `n < 2^K` with exactly `m` trailing ones is `sum_H 2^(popcount H) =
  3^(K-1-m)`. ∎ The sides `1, 3, 7, 15, 31` in the first `32` rows are the
  Sierpinski sides of the five known Fermat primes' rows: the side-`(2^k - 1)`
  triangle opens under the all-ones row `2^k - 1` (value `prod_(i<k) F_i`) and
  its top zero row is row `2^k`, whose binary value is `F_k = 1 0...0 1` (a
  first draft said the triangle opens under row `2^k`); every side `2^m - 1`
  with `m >= 2` is `3 mod 4`. In the primes' triangle (S9, audited) the sides are
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
`P_7`, and `|Aut(T_k)| = 21` for `k = 4, 5, 6, 7, 8` (`n = 15, 31, 63, 127, 255`;
the order-`255` search took `648` s and `35140` nodes). For `T_5`: element orders `{3: 14, 7: 6}`,
i.e. the Frobenius group `Z_7 x| Z_3` of the heptagon, with orbits of sizes
`7, 1, 7, 1, 7, 1, 7`. `T_5` is **not** isomorphic to the Paley tournament
`P_31` (`|Aut(P_31)| = 465`), and `T_7` is not `P_127` (`|Aut(P_127)| = 8001`) although its 4-vertex census
(`4340` with a source, `13020` strong, `4340` with a sink, `9765` transitive)
equals `P_31`'s (the census did not separate them; the automorphism
groups do). So: the zero-triangle sides `2^k - 1` of the sea are exactly the orders
of a tower of doubly regular tournaments generated by the sea's own doubling,
the Paley heptagon is its `k = 3` level, and as far as computed (`k <= 8` here, `k = 9` in the audit) the tower
keeps the heptagon's `21` symmetries and nothing else; the audit also proved
that the stabilizer of `0'` in `Aut(T_(k+1))` consists exactly of the diagonal
extensions, so `Aut(T_k) = F_21` for all `k` is equivalent to `0'` being fixed
by every automorphism ([HYP-9162](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md)).
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

## 5. THM-4511: the wall theorem for size-4 defects, and the exact extinction law (PROVED + AUDITED + FINITE-EXACT)

**Theorem.** Run `a_(r+1)(i) = |a_r(i) - a_r(i+1)|` from a row with
`a_0(0) = 1` and `a_0(i) in {0, 2, 4}` for every `i >= 1` (finitely or
infinitely many columns), containing at least one `4`. Let `F` be the first
column with `a_0(F) = 4` and `c` the largest column with `1 <= c < F` and
`a_0(c) = 2`, if it exists. Then

1. if `c` exists, every column `1 <= i <= c` stays in `{0, 2}` for all rows
   and the leading entry stays `1`: **no entry `>= 4` ever crosses the wall
   `c`**, whatever `0/2/4` pattern lies to the right of `c`;
2. if `c` does not exist (`a_0(1..F-1) = 0`), then `a_(F-1)(1) = 4` and
   `a_F(0) = 3`: the leading `1` is destroyed at row `F`.

Hence a row of `0`s, `2`s and `4`s destroys the leading `1` iff its first `4`
is preceded by zeros only; for a lone size-`4` defect at distance `F` in a
uniform random sea the probability is exactly `2^(1-F)`, and re-emissions
from the stationary copy never add anything. (My first draft assumed a
single `4` and claimed that a second `4` could reopen the wall; the audit
observed that the proof never uses the assumption, and that mechanism does
not exist.)

*Proof.* Column `i` depends only on columns `i` and `i + 1` (the rule looks
right), and cell `(r, i)` only on `a_0(i..i+r)`, so for a statement about
rows `<= R` and columns `<= C` the row may be truncated at column `C + R`;
assume finitely many `4`s. Call a column sequence *good* if it takes values
in `{0, 4}` up to some row and values in `{0, 2}` from the next row on
(either phase may be empty). Every column right of the last `4` is in
`{0, 2}` for ever, hence good. If column `i + 1` is good, column `i` is good
whatever its initial value: while column `i + 1` is in its `{0, 4}` phase,
column `i` stays in `{0, 4}` if it started there and stays `2` if it started
at `2` (`|2 - 0| = |2 - 4| = 2`); from the first row where column `i + 1` is
in `{0, 2}`, a column `i` at `0` or `2` is in `{0, 2}` from then on, and a
column `i` at `4` stays `4` while its neighbour is `0`, becomes `2` at the
neighbour's first `2`, and is then in `{0, 2}`. So every column is good, by
induction from the right. Columns `c + 1, ..., F - 1` start at `0`, so their
`{0, 4}` phases begin at row `0`; column `c` starts at `2`, stays `2` while
column `c + 1` is in `{0, 4}`, becomes `0` at the first row where column
`c + 1` is `2`, and stays in `{0, 2}` for ever; hence (induction leftward) so
do columns `1, ..., c - 1`, and `|1 - a(1)| = 1`. For 2: on `{0, 4}` the rule
is XOR in units of `4`, so the cone of the first `4` over the zeros is
Pascal mod `2` and its left edge `(s, F - s)` is `4` for every `s <= F - 1`;
that cone is `[F - s, F]` and never reaches column `F + 1`, so everything to
the right is irrelevant. ∎

*Verification.* Exhaustive enumeration over all `2^(F-1+R)` sea contexts
(`F - 1` cells left, `R = 10` right, `R = 9` at `F = 12`, tableau run for
`F + R` rows) for `F = 3..12`:
the minimum column ever holding a `4` equals `F - z` (`z` = trailing zeros of
the left sea) in every configuration, and `p_4(F) = 2^(1-F)` to all printed
digits; Monte Carlo with `4 * 10^6` samples: `p_4(14) = 1.267e-4` vs
`1.221e-4`, `p_4(16) = 2.875e-5` vs `3.052e-5` (Poisson noise). The audit
verified the multi-`4` form exhaustively on all `3^L` rows of length `L <= 13`
and on `200000` random rows of length `40` with about eight `4`s each.
Canon file: [THM-4511](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md).

*What the theorem says about Gilbreath.* **Corollary (light cone).** If row
`r` has `a_r(1..t)` in `{0, 2, 4}`, the leading `1` survives at least until
row `r + t`, unless the first `4` among those cells is preceded by zeros
only, in which case it is destroyed at row `r + F` (truncate the row at
column `t`; the finite-row theorem applies). So **all danger to the leading
`1` comes from entries `>= 6`**, and a finite Gilbreath triangle is decided at
its first all-`{0,2,4}` row `r_4`: the leading `1` survives every later row
iff the first `4` of row `r_4` is preceded by a `2`. For the primes below
`200000`: `r_4 = 59` against the first all-`0/2` row `r_2 = 65`; below `10^6`:
`r_4 = r_2 = 95`; below `10^7`: `r_4 = 132`, `r_2 = 135`; in each case the
certificate holds, and no computed row has a first `4` preceded by zeros only
(`gilbreath_certificate_20260926.py`). The light-cone certificate per row is
modest in the first rows (`t = 7` at row `1`, `27` at row `2`, `119` at row
`9`, `1804` at row `30`, `3366` at row `50`, all certified) because entries
`>= 6` sit near the edge there. The random-model numbers of the S9 heuristic
(front-only tails: for the primes below `200000` the `29` fresh fronts of
rows `1`-`64`, `23` of size `4`, give a risk sum `0.258` dominated by the
row-`1` term `1/4`) are **not** consequences of the theorem, since rows
`1`-`64` contain thousands of entries `>= 6`; a first draft called the
row-`1` defect "the riskiest moment", a model artefact the audit removed.
**The theorem does not cover sizes `>= 6`**: for a lone `6` the wall breaks
to a `4` (`|6 - 2| = 4`) and the finite-context (`R = 10`) extinction
probability exceeds the front-only tail `F 2^(1-F)` by `6`-`10` per cent
(`F = 8`: `0.067993` vs `0.062500`; `F = 12`: `0.006338` vs `0.005859`; these
values are `R`-stable to seven digits); for a lone `8` by `4`-`8` per cent
(`F = 8`: `0.244003` at `R = 10`, still moving to `0.244017` at `R = 14`, vs
`0.226563`), the excess being the stationary copy's re-emitted fronts, which
for `d >= 6` can pass a wall the first front has already weakened.

## 6. Bearing on the two conjectures (honest)

* **Gilbreath.** The Fermat tower is the structure of the *safe* region; the
  conjecture is a statement about the frontier, i.e. about defect extinction.
  THM-4511 settles the size-`4` question exactly (any number of `4`s: the
  wall, the light-cone certificate); sizes `>= 6` are governed by the same
  left-looking mechanism but no closed law is proved. The five Fermat primes and the five
  solids do not enter the conjecture.
* **Collatz.** `3 = F_0` is the multiplier; the base-`6` automaton (S9) is
  nonlinear; the `F_k n + 1` maps for `k >= 1` have positive drift. No
  reduction in either direction (S9 and forest, unchanged).
* **Tournaments.** The Mersenne orders `2^k - 1` of the sea's zero triangles
  carry the doubling tower `T_k`, doubly regular with `Aut = F_21` for
  `k = 3..8` (`k = 9` in the audit); the Paley heptagon is its third level; THM-871's Fermat-rung
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

## 8. Independent audit (2026-09-26)

Auditor subagent: `04-computation/experiments/gilbreath_fermat_platonic_20260926_audit.py`
(sha256 `365f22327926b4974fde93431bc5a26c18b6af24bd516d65397ee8f886062a0a`) ->
`05-knowledge/results/gilbreath_fermat_platonic_20260926_audit.out`
(sha256 `a969b97e330299e19b4ed8cd761c5b0aa604d6fe389c5b3883800d3e89ca33e5`),
everything re-implemented from scratch. Verdict: HAS GAPS -> repaired.
CONFIRMED: every quoted number (frontier, kernel checks, triangle counts,
automorphism orders by two independent methods and `|Aut(T_9)| = 21`,
`Aut(P_7), Aut(P_11), Aut(P_31)`, the census, group orders and
identifications with sympy, Klein subgroups, the extinction table and the
seeded Monte Carlo, the fresh-front statistics, von Staudt-Clausen against
actual Bernoulli denominators to `2^9`), every written proof, every citation
(THM-871, THM-870, HYP-3771, HYP-3772, forest, S9), and the classical facts.
CORRECTED (applied above): (1) the scope of THM-4511, which holds for any
number of `4`s (the first draft's "a second `4` reopens the wall" cannot
happen), and with it the "interacting defects open" wording everywhere; (2)
the `3^(K-1-m)` count and the completeness of the triangle decomposition were
labelled PROVED without a proof on the page (now proved in section 2); (3)
"whatever lies to its right" now reads "whatever `0/2/4` pattern lies to the
right" (`1 2 4 10` is a counterexample to the literal reading); (4) the
random-model risk numbers for the primes were presented as consequences of
the theorem, which they are not (rows `1`-`64` contain entries `>= 6`);
replaced by the light-cone corollary and the `r_4` certificates; (5) the
size-`8` "exact" values are finite-context values still moving in the fifth
digit, and `F = 12` used `R = 9`; (6) the side-`(2^k - 1)` triangle opens under
row `2^k - 1`, with row `2^k` as its top zero row; (7) the summary stated the
`F_21` conjecture as fact; (8) column `0` excluded from "columns `<= c`"; (9)
`T_5 != P_31` is established, not conditional on HYP-9162; (10) the word
"once" in the definition of a good column. Novelty of THM-4511 relative to
Odlyzko 1993 is not established (the auditor recalls the light-cone and
absorption remarks there, not the all-time wall statement or the exact law).
