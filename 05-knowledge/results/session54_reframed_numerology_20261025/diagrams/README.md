# Session54: two exact explanatory figures

This supplementary illustration packet adds no theorem, relation, computation-based conjecture, or originality claim. The mathematical assertions are sourced to the four pinned Markdown texts in `sources/`.

- `digit-budget.svg`: fixed digit budget, the correct r=1 initial state, finite-state transfer, exact residue counts, and the distinction between the admissible completed exponent set and ordinary nonnegative exponents. A state arrow denotes an admissible next-digit transition, not ordinary global survival. Budgets b≥2 are combinatorial extensions only.
- `boundary-incidence.svg`: the actual commuting native-boundary square under the adjusted field coding, followed by a visibly separate typed trace-incidence construction. Three scalar vectors represent one point, not a line. The lift is a complete-bipartite blow-up, not a graph cover. Converse with a type switch is polarity; point-only inversion is not a collineation.

The standalone SVGs use embedded text, no third-party assets, no scripts, and no external fonts or resources. Source/legend labels accompany both figures. They are designed for full-width viewing (1200 px); zoom rather than shrinking the text on a phone.

## Reproduction and scope

Run `python3 -B replay.py --self-test` from any directory. It checks exact file membership and SHA256 digests, regenerates both SVGs and counts byte-for-byte in a temporary directory, checks transfer/binomial counts through r=10 for budgets 0…6, independently checks the field trace flags and degrees, and verifies the essential caution labels. It rejects extra, missing and modified files and compares normal and optimized Python output. All replay/generation dependencies are Python standard library only.

`generate.py` deterministically emits both SVGs and `counts.json`. The replay is diagram sanity/integrity checking, not an independent verification of every source theorem or the native blocker bijection. Those are proved and replayed by the source packets. In particular the equalities involving the native boundary are sourced, not inferred merely from the number 63.

Source copies:
- `sources/cross-base.md`: `session54-cross-base-fractal/REPORT.md`
- `sources/budget-transfer.md`: `session54-fractal-budget-transfer/REPORT.md`
- `sources/boundary.md`: `session54-boundary-projectivization/README.md`
- `sources/two-sheet.md`: `session54-boundary-projectivization/TWO-SHEET.md`

For optional pixel previews, use Inkscape, e.g. `inkscape digit-budget.svg --export-type=png --export-filename=/tmp/digit-budget.png`. Both figures were rendered with Inkscape and inspected visually for legibility, clipping and arrow-label collisions. PNG previews are not necessary for portable replay and are not included in the text-only packet.
