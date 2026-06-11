# HYP-2439 - The order-5 fixed projection is a marked Type II length-16 gate

**Status:** CLAIMED conditional reduction / proof route. Codex 2026-06-11.

Let `C` be a putative extremal binary Type II `[72,36,16]` code and let
`sigma in Aut(C)` have prime order `5` and cycle type `5-(14,2)`.
Project the fixed subcode `C^sigma` onto the fourteen 5-cycles plus the two
fixed coordinates. The standard odd-prime fixed-code theorem gives a binary
self-dual length-16 projected code, and doubly-evenness is retained because
`5 == 1 (mod 4)`.

Thus the fixed projection is one of the two Type II `[16,8]` codes:

1. `e8+e8`;
2. `d16+`.

Mark the two projected coordinates that came from actual fixed coordinates.
For a projected word with `a` marked 5-cycle coordinates and `b` fixed
coordinates, the original weight is

```text
5a + b.
```

Extremality says every nonzero fixed word must have weighted length at least
`16`. Therefore a projected tetrad containing both fixed marks is fatal: it
lifts to original weight `5+5+1+1=12`.

The exact atlas `04-computation/order5_fixed_projection_72_codex.py` verifies
the resulting tiny gate:

- in `e8+e8`, the valid marked pairs are exactly the `64` pairs split across
  the two `e8` summands;
- in `d16+`, every one of the `120` coordinate pairs is covered by a tetrad,
  so no marking is valid.

## Consequence

If an extremal `[72,36,16]` code admits an order-5 automorphism of type
`5-(14,2)`, then its fixed projection is forced to be `e8+e8`, with the two
fixed coordinates lying in different `e8` summands.

Equivalently, the fourteen 5-cycles split into two heptads, one heptad attached
to each fixed point.

## Proof Obligations

1. Cite or reprove the odd-prime fixed-code projection theorem in the exact
   `5-(14,2)` setting.
2. Promote the script's two Type II `[16,8]` enumerations to a short hand
   proof using the known classification `e8+e8`, `d16+`.
3. Couple the fixed projection to the nonfixed `F_16` component in the
   Feulner-Nebe style module decomposition.

## Links

- Supports HYP-2415, HYP-2425, HYP-2430, OPEN-Q-061, OPEN-Q-063.
- Computation: `04-computation/order5_fixed_projection_72_codex.py`.
- Result: `05-knowledge/results/order5_fixed_projection_72_codex.out`.
