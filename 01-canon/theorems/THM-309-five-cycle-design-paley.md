# THM-309: 5-Cycle Design in Paley Tournaments

**Status:** PROVED; SYMMETRY PROOF REPAIRED AND EXTENDED TO ALL CYCLE LENGTHS
**Verified computationally:** p = 7, 11, 19, 23 for k = 5; p = 7, 11
for every 3 <= k <= 7
**Found by:** opus-2026-04-05-S28
**Verified in:** `04-computation/five_cycle_designs_s28.py`, `cycle_design_formulas_s28.py`
**Repair companion:** `04-computation/paley_all_cycle_lengths_design_thm309.py`
with output `05-knowledge/results/paley_all_cycle_lengths_design_thm309.out`

## Statement

For the Paley tournament P_p (p ≡ 3 mod 4 prime), the directed 5-cycles form a **2-(p, 5, λ₅) design** with:

**λ₅ = (p+1)(p-2)(p-3)/8**

Every pair of vertices participates in exactly λ₅ directed 5-cycles.

The total count of directed 5-cycles is:

**c₅ = p(p²-1)(p-2)(p-3)/160**

Here cycles are counted modulo cyclic rotation, and their vertex sets are
blocks with multiplicity when one support carries several directed cycles.

More generally, for every `3<=k<=p`, the multiset of simple directed
`k`-cycles is a **2-(p,k,lambda_k) design**, where, if `c_k` is their total
number,

**lambda_k = c_k C(k,2)/C(p,2).**

## Proof Sketch

The required symmetry is on **unordered pairs**, not ordered pairs.  Let

```text
G={x -> ax+b : a is a nonzero square in F_p, b in F_p}.
```

Every element of `G` preserves the Paley orientation.  The group has order
`p(p-1)/2`, exactly the number of unordered vertex pairs.  Given source pair
`{x,y}` and target pair `{u,v}`, one of
`(v-u)/(y-x)` and its negative is a square, and exactly one is: `-1` is a
nonsquare because `p=3 mod 4`.  Choosing that ordering supplies the unique
square-affine map between the pairs.  Thus `G` acts sharply transitively on
unordered pairs.

The multiset of simple directed `k`-cycles is invariant under `G`.  Every
unordered pair therefore has the same cycle incidence, proving the all-`k`
2-design statement.  Double-counting `(cycle, contained pair)` gives
`lambda_k=c_k C(k,2)/C(p,2)` and hence the displayed `lambda_5`.

The formula for c₅ comes from the Jacobi sum evaluation: the number of directed 5-cycles through a fixed vertex is c₅ × 5/p, and this equals a sum of Legendre symbol products over 4-tuples of field elements, which evaluates to the given closed form via the Hasse-Davenport relation.

The earlier proof incorrectly identified the automorphism group with the full
`AGL(1,F_p)` and invoked 2-transitivity.  Nonsquare multipliers reverse the
tournament, so that premise is false.  The sharp unordered-pair action above
is the precise repair and is also why even cycle lengths require no new
character-sum identity once their total counts are known.

## Relation to 3-Cycle Design

The ratio λ₅/λ₃ = (p-2)(p-3)/2 = C(p-2, 2), where λ₃ = (p+1)/4 is the 3-cycle design parameter (THM-211). This means:

**Each pair's 5-cycle incidence is exactly C(p-2, 2) times its 3-cycle incidence.**

## Verified Cases

| p | c₃ | λ₃ | c₅ | λ₅ | λ₅/λ₃ |
|---|-----|-----|------|------|--------|
| 7 | 14 | 2 | 42 | 20 | 10 |
| 11 | 55 | 3 | 594 | 108 | 36 |
| 19 | 285 | 5 | 11628 | 680 | 136 |
| 23 | 506 | 6 | 31878 | 1260 | 210 |

## Special Case: P_7

At p = 7, the 5-cycles have extraordinary additional structure:

- c₅ = 42 = 2 × C(7,5): every 5-element subset supports exactly 2 directed 5-cycles
- Every induced 5-vertex subtournament of P_7 is REGULAR (scores (2,2,2,2,2))
- The 5-cycles form a **4-(7, 5, 6) design** (not just a 2-design)
- This equals 2 copies of the trivial complete design on 5-subsets of 7 points
- Each complement pair (the 2 vertices NOT in a 5-cycle) appears exactly 2 times

At p ≥ 11, the 5-cycles form only a 2-design (the 3-design and 4-design properties fail).

## Exact repair audit

The repair companion enumerates simple cycles modulo rotation without floating
point or third-party packages.  It checks the sharp unordered-pair orbit for
`p=7,11,19,23`, all lengths `3..7` for `p=7,11`, and the existing 5-cycle
formulas for all four primes.  Ordinary and optimized runs byte-match the
stored transcript.  Script/output/semantic SHA-256 are
`40182a73db81a1ce9991eaa2fb7d4fbb5593f4b79430be19e491c35cfe9f3730`,
`478d2d54446b2b70c0a208bb66b915a9515b759feec2ca56acccd4b80032c94e`,
and `81e1de44977eb1757a0bb6fe2a9e9cb75a010b7e53b493ee4e49147f117fc5c4`.

This is an unordered-pair multidesign statement.  It does not assert uniform
incidence on ordered pairs or uniform multiplicity on every `k`-vertex support.

## Non-Paley Comparison

The interval tournament {1,2,3} mod 7 (circulant with non-QR connection set) has:
- 3-cycles: NOT a 2-design (pair distribution {1:7, 2:7, 3:7})
- 5-cycles: NOT a 2-design (pair distribution {13:14, 14:7})

Only the QR structure of Paley creates the design uniformity.

## Related

- THM-211: 3-cycle BIBD for Paley
- THM-305: v_2(T(n)) = (n-1)/2
- The Barvinok-BIBD-QR triangle (07-reflections/barvinok-bibd-qr-triangle.md)
