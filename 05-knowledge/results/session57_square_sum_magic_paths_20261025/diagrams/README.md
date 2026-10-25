# Session 57 exact diagrams

Two deterministic SVG illustrations. Run `python3 generate.py` and `python3 check.py`; `python3 -O check.py` performs the same checks. The checker uses explicit exceptions, exact integers/rationals, exhaustive six-vertex path enumeration, pinned source hashes, and a byte-identical SVG regeneration check. No third-party Python package is needed.

- `edge-label-repair.svg`: the same F59 magic set moves from entry vertices to additive AP carrier edge labels; color codes P-only, Q-only and shared labels. This is a cross-skeleton, not an asserted induced finite-field graph.
- `square-totals-not-magic.svg`: the integer square-total variant, explicitly not magic. The unrestricted magic-square-of-squares existence question remains unresolved by these results.

Sources are copied and pinned in SOURCE_PINS.json. Character and AP text come from frozen-second. The variant source was absent from that packet; its current source report and certificate were locally snapshotted and pinned before drawing. Independent arithmetic is in check.py.

Local previews were rendered with Inkscape and visually inspected. Publish the SVGs, Python sources, JSON results, README and pinned text sources only. Do not publish PNG previews or render.log. No image generation, remote writes, Library operations, or older research-file edits occurred.

Optional PNG QA command: `inkscape edge-label-repair.svg --export-type=png --export-filename=edge-label-repair-preview.png` (likewise for square-totals-not-magic.svg).
