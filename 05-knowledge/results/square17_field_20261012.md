# The seventeen-square packing constant generates a degree-18 field with Galois group S_18; its side formula is the mirror image of the eleven-square one

2026-10-12, opus-2026-10-12-S24 (third addendum). **CITED (trusted per the
owner, not re-verified):** the optimality theorem of
[sq17-lean](https://github.com/ofitzharding/sq17-lean) — seventeen unit squares
pack in a square of side `L` iff `L >= 4.675530093604551...`, the exact side of
Bidwell's 1998 packing (Lean 4.34.1, pinned Mathlib, no `sorryAx`, native
decision procedures, same trust model as the eleven-square formalization of
Ahmed and collaborators). **PROVED + FINITE-EXACT (this note):** everything
about the field of that constant. **CONNECTION:** with the repo's
[eleven-square octic note](square11_octic_langlands_collatz_20261006.md)
(opus-2026-10-06-S16). Self-audited; no independent peer audit.

## 1. The constant, exactly

The formalization defines the side from a point `(t, beta)` in a `1e-12` box
around `(0.36204389331432338..., 0.33094867472273991...)` that is the unique
root there of two integer polynomials (`Sq17Lean/EndpointRootData.lean`;
`polynomial0` of degree 6, `polynomial1` of degree 8, both transcribed and
validated by the root to `1e-58`), and

    side = (6 + 4 t) / (1 + 2 t - t^2).

With `t = tan(theta/2)`, `theta = 39.80495898...` degrees is the tilt of the
rotated squares, and

    s(17) = (3 + 3 cos theta + 2 sin theta) / (cos theta + sin theta).

Eliminating `beta` (resultant of degree 28 = three squared linear factors, a
squared quadratic, and one factor of degree 18) gives the minimal polynomial
of `t`,

    t^18 - 2t^17 - 11t^16 - 24t^14 - 112t^13 + 48t^12 - 16t^11 - 622t^10 - 308t^9
    + 194t^8 - 480t^7 - 224t^6 - 576t^5 + 24t^4 - 16t^3 + 101t^2 - 26t + 1,

irreducible over `Q`, and the side's minimal polynomial is the degree-18
polynomial `4775 s^18 - 190430 s^17 + ... + 26662976550 s - 2631254953`
(leading coefficient `5^2 * 191`; full coefficients in the output). So
`[Q(s(17)) : Q] = 18` and `Q(s(17)) = Q(t)`.

## 2. The Galois group is S_18 (PROVED)

Dedekind cycle types of the `t`-polynomial modulo primes: `(5, 13)` at `p = 3`
(a 13-cycle; `13` is prime, `13 > 18/2` makes the group primitive, and
`13 <= 18 - 3` lets Jordan's theorem give `A_18`), and `(18)` at `p = 7` (an
odd permutation). Hence `Gal = S_18`, as generic as possible, like the
eleven-square octic's `S_8`. Over 237 primes below 1,500 there are 98 distinct
cycle types; root counts cross-check the factorization at every prime. The
discriminant of the `t`-polynomial is
`-2^181 * 5 * 154943^2 * 388183 * 193319963`; the real roots are
`0.0470, 0.2489, 0.3620 (= t), 4.7521`, signature `(4, 7)`.

## 3. The relation with the eleven-square packing

| | eleven squares (Trump 1979, Lean: Ahmed et al.) | seventeen squares (Bidwell 1998, Lean: sq17-lean) |
|---|---|---|
| side | `3.8770835900...` | `4.6755300936...` |
| tilt angle | `40.1819373...` deg | `39.8049590...` deg |
| `tan(theta/2)` | `u = 0.3657693...`, root of a degree-8 polynomial | `t = 0.3620439...`, root of a degree-18 polynomial |
| side formula | `(2 + 2 cos + 3 sin)/(cos + sin) = (6u + 4)/(1 + 2u - u^2)` | `(3 + 3 cos + 2 sin)/(cos + sin) = (6 + 4t)/(1 + 2t - t^2)` |
| field degree | 8 | 18 |
| Galois group | `S_8` (the octic note) | `S_18` (this note) |
| discriminant primes | `2, 5, 31, 104053` (field) | `2, 5, 154943, 388183, 193319963` (polynomial) |

The two formulas are mirror images: the coefficients `(2, 3)` of the eleven
packing become `(3, 2)` for seventeen — the same quarter-turn frame
`cos theta + sin theta` in the denominator and a numerator `p(1 + cos) + q sin`
with `(p, q)` swapped. Both tilts sit near 40 degrees, less than a degree
apart, yet the fields share no prime of ramification beyond `2` and `5`,
have no common subfield (symmetric Galois groups leave no proper subfields in
`Q(u)` and `Q(t)`), and are of different degrees; the near-coincidence of
the angles is numerical, not algebraic. In the octic note's Langlands
ladder, `Q(s(11))` was a `GL(7)` Artin case (the standard representation of
`S_8`); `Q(s(17))` is a `GL(17)` Artin case, the standard representation of
`S_18`, with the same unconditional piece only at `GL(1)` (the sign character,
here the Kronecker character of `Q(sqrt(disc))`). The earlier proven optima
`s(5)`, `s(10)` live in `Q(sqrt 2)`: the pattern "a packing that needs a
tilted sub-block has a generic field whose degree grows with the number of
coupled contacts" now has two data points, degrees 8 and 18.

What is *not* claimed: any relation between the two Lean proofs beyond the
shared trust model (native decision procedures) that the sq17-lean README
itself cites from the eleven-square work; any arithmetic link between the two
fields; or an independent check of either optimality theorem.

## 4. Reproduction

    python3 -B 04-computation/experiments/square17_field_20261012.py

[Script](../../04-computation/experiments/square17_field_20261012.py),
[output](square17_field_20261012.json): transcription check of the two Lean
polynomials at the stated root, the resultants, both minimal polynomials, the
cycle-type census with root-count cross-checks, the discriminant, the
signature, the angle identities and the mirrored side formulas.
