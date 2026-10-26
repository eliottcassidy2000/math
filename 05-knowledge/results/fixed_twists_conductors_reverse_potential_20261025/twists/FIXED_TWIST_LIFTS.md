# Fixed-twist rational lifts on the hourglass boundary: a conic on the central line, genus five on the outer lines

2026-10-25, klein. **PROVED** pointwise criterion and central-line theorem; **CITED** genus-five
components (Gale-cube / Session60) and Faltings; **FINITE-EXACT** Hilbert-symbol decisions against
brute force (`check_central_conic.py`, `check_boundary_lifts.py`; receipts byte-identical under `-O`).
This is the "residual squareclasses for prescribed nontrivial twists" target of
[normseven_boundary](../../normseven_boundary_20261025/README.md), made decidable where it can be.

## 1. Pointwise criterion (any field k, char k != 2)

Forms `L_0..L_6 = (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z)`, twist `d = (d_0..d_6) in (k^x)^7` modulo
common scaling, surface `X_d : s_i^2 = d_i L_i` (projective root surface). For `y in P^2(k)` let
`r(y)` be the number of vanishing forms. Then `y` lifts to `X_d(k)` iff the nonzero values
`d_i L_i(y)` lie in one class of `k^x/k^{x2}`, and then it has exactly `2^{6-r(y)}` lifts.
(This is the fibre count in THM-4765's proof; it is also BOUNDARY §5's residual torsor `V_0/W_0`
of dimension `r`, written out as `r` squareclass equalities.)

## 2. Central line theorem

Over `E = 0` with `x z (x+z) != 0` the six outer values are `d_1 x, -d_2 x, d_3 z, -d_4 z, d_5(x+z),
-d_6(x+z)`. Hence a lift exists iff

  (i)  `-d_1 d_2`, `-d_3 d_4`, `-d_5 d_6` are squares in `k`  (the three parallel-pair port
       characters are trivial on the twist; these are BOUNDARY §5a's necessary tests), and
  (ii) `d_1 d_3 x z` and `d_1 d_5 x(x+z)` are squares.

Writing `x = d_1 d_3 s^2`, `z = t^2` solves the first condition of (ii), and the second becomes the
conic

    C_d :  w^2 = d_1 d_5 d_3^2 s^2 + d_3 d_5 t^2.                                     (C)

**Theorem.** `X_d(k)` has a point over the open central line iff (i) holds and `C_d(k) != ∅`
(a conic with one rational point has infinitely many, and the finitely many degenerate ones with
`st(d_1d_3 s^2 + t^2) = 0` can be avoided). Over `Q` this is decided by Legendre/Hilbert symbols of
`(d_1 d_5 d_3^2, d_3 d_5)` at `2`, `oo` and the odd primes dividing the coefficients.
Consequences:
- untwisted `d = 1`: (i) fails (`-1` is not a square), recovering the central case of BOUNDARY 5a;
- Pythagorean twist `d = (1,1,-1,1,-1,1,-1)`: `C_d` is the unit circle `w^2 = s^2 + t^2`, so the
  rational boundary lifts over `E = 0` are parametrized by Pythagorean triples; `(3,4,5)` gives the
  base point `(0:9:16)` and the point `(0:3:3:4:4:5:5) in X_d(Q)`, which satisfies the three
  pair-sum quadrics `s_1^2/d_1 + s_2^2/d_2 = 2 s_0^2/d_0` etc. and the mixed relation
  `s_5^2/d_5 - s_1^2/d_1 - s_3^2/d_3 = -s_0^2/d_0` (`25 - 9 - 16 = 0`);
- `d = (1,1,-1,1,-1,7,-7)`: (i) holds, `C_d: w^2 = 7s^2 + 7t^2` has no rational point
  (`(7,7)_7 = -1`), no lifts; `d = (1,2,-2,3,-3,6,-6)`: `w^2 = 108 s^2 + 18 t^2` fails at `3`.
Receipt `conic-normal.txt`: the Hilbert decision agrees with brute-force existence (height ≤ 25) on
all 300 random twists with entries in `{±1,±2,±3,±5,±6}`, and on eight displayed twists.

## 3. Outer lines and points

On an outer line the lift locus is the union of two curves of genus five (the branch components of
[GENUS_FIVE_ELLIPTIC_FACTORS](../../six_seven_gale_cube_20261025/GENUS_FIVE_ELLIPTIC_FACTORS.md) and
[BRANCH_CUBE_AND_DOUBLE_COVER](../../six_seven_gale_cube_20261025/BRANCH_CUBE_AND_DOUBLE_COVER.md),
twisted by `d`; the twist changes fields of definition, not genus). By Faltings (**CITED**) each has
finitely many points over a number field. The twelve double points and three triple points are
finite squareclass tests: e.g. `(0:0:1)` lifts iff `d_3 ~ -d_4 ~ d_5 ~ -d_6`, so the twist
`(1,1,1,1,-1,1,-1)` lifts it while `d = 1` does not.

**Dichotomy.** For a fixed twist over a number field, the central line carries either no rational
lift or infinitely many (a conic), each outer line finitely many (genus five), and the fifteen
special points finitely many; the untwisted surface has none at all. Over `Q` with height ≤ 12 the
Pythagorean twist lifts exactly 4 of 6521 points of `P^2`, all on the boundary (receipt
`lifts-normal.txt`).

## 4. Ports and what is not claimed

Conditions (i) are the three parallel-pair characters `S_{{1,2}}, S_{{3,4}}, S_{{5,6}}` of THM-4765's
Walsh spectrum evaluated on the twist; (C) encodes the two remaining independent even characters
restricted to `E = 0`. Nothing is claimed about interior rational points, about the hourglass over
`F_q` beyond THM-4765, or about integer magic squares of squares; session63's Theorem B (fixed twists
over extensions are not fresh twists) is consistent with and independent of this note.
