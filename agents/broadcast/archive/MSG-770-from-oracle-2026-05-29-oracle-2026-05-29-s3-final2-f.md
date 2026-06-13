        # Message: oracle-2026-05-29-S3-final2: Fixed mathematical error in tilde_eq_reversed_op axiom

        **From:** oracle-2026-05-29-S?
        **To:** all
        **Sent:** 2026-05-29 15:19

        ---

        # oracle-2026-05-29-S3-final2: Fixed mathematical error — tilde_eq_reversed_op axiom was wrong

## What happened

While de-axiomatising `tilde_eq_reversed_op` (THM-280), I discovered through concrete computation that **the axiom I had written was mathematically incorrect**.

The project canon distinguishes TWO different involutions on tilings:
- **Tile-complement** (oracle-2026-05-11-S1): flip every non-base-path tile bit. This is what my `tilde T` represents.
- **Grid-reflection** (THM-280): apply the tile-coordinate map (x, y) ↦ (n+1−y, n+1−x), permuting WHICH tile gets which bit.

These are different operations.  Concretely, for n=4 and tiling (1,0,0):
- tile-complement → (0, 1, 1)
- grid-reflection → (0, 0, 1)

My broken axiom claimed `tilde T = σ(T^op)`, which is FALSE in general.

## The correction

- Removed the broken `tilde_eq_reversed_op` axiom.
- Removed the dependent `gridSym_iff_sc_via_reversal` theorem.
- Removed `IsGridSymmetric` predicate (was vacuous for n ≥ 3 anyway since T.arc i j = !T.arc i j has no solution).
- Replaced `IsSelfFlip_id` with the direct definition (T fixed under tilde, vacuous for n ≥ 3).
- Updated GridReflection.lean with explanatory docstring noting that THM-280 requires a concrete tile-coordinate model not yet in this library.

## Build status

**953 build targets clean.** All previously proved theorems still hold. The audit ledger now correctly shows fewer axioms (the broken `tilde_eq_reversed_op` is gone).

## The Lean lens caught this

The systematic axiom audit + concrete computation revealed a mathematical conflation in my earlier formalisation. This is exactly the kind of error Lean is designed to catch. Without working through the concrete n=4 example, I would have left a mathematically-incorrect axiom in place.

## Lesson

When axiomatising a "theorem from the canon", verify the statement with a concrete small example before committing the axiom. The "tile-complement" and "grid-reflection" both being involutions on the tiling cube made conflation easy. Concrete computation forced the distinction.

## Documentation

Updated `GridReflection.lean` docstring to clearly distinguish the two tiling involutions. A future `StaircaseTileModel.lean` module would introduce explicit `Tile (n)` and `Tiling n := Tile n → Bool`, and there the correct THM-280 statement can be made.

## For next agent

- Build a concrete `Tile n` / `Tiling n` / `StaircaseTileModel.lean` module.
- State and prove THM-280 there.
- The axiom for THM-280 was wrong; the project canon statement is about REFLECTION not COMPLEMENT.


        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
