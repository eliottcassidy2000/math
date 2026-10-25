# Session 56: exact query and cycle bridge figures

Two deterministic, accessible, self-contained SVG figures. No bitmap generation, embedded raster images, network assets, or third-party runtime dependencies. SVGs are the only published visual assets. QA PNGs are excluded and held in a separate temporary directory.

- `shared-choice-overlap.svg`: overlapping CRT blocks [2,5] and [5,3], original F5 state identity versus independent resampling, and the distinct fixed-mask/correlated-copy/independent-law notions of idempotence. The 128 fabricated all-distinct outcomes count independent ordered block-row pairs, not necessarily distinct resulting masks. The 112 true arrays have zero all-distinct members.
- `path-locked-score-fibre.svg`: every arc in each of the two full labelled tournaments is shown. Dashed arcs are the exact locked Hamilton path; teal arcs form the reversible residual quadrilateral; blue residual arcs are unchanged. Arrow direction is authoritative. Lines may cross without forming a vertex. H counts directed Hamilton paths as vertex permutations. The exact full ordered score is (3,2,1,2,2), with labels 0 through 4 fixed. No legal directed residual triangle occurs. All-size forest-lock connectivity is a theorem of the cited source packet, not inferred from this five-vertex audit. Native complement wiring admissibility is not established by that theorem.

## Portable replay

Run `python3 generate.py --verify` and `python3 -O generate.py --verify` from any working directory. Standard library only. The generator independently enumerates local modular magic-array states, both block-row sets, their independent pairing and shared-ID join. It exhausts all 64 fixed-path tournaments, finds the entire score fibre, enumerates all vertex permutations for H, checks the directed quadrilateral reversal and absence of residual triangles, compares the pinned certificates, and byte-checks both SVGs and `counts.json`. To regenerate, omit `--verify`.

Run `sha256sum -c SHA256SUMS` for file integrity. `SOURCE_PINS.json` gives SHA-256 hashes and original source locations; copies are bundled for portability. Original sources are Session 56 correlated-choice gluing and path-locked cycle-moves packets; the local enumerator is their pinned Session 55 helper.

Rendered at natural resolution with Inkscape and visually inspected on 25 October 2026: no clipped text, missing glyphs, overlapping labels or ambiguous arrowheads. PNG previews are QA-only and not in this directory. The optional tree-crossing inset was omitted to keep both primary figures readable.
