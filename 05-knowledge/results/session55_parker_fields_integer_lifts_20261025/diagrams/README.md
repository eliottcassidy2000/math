# Session 55: exact arithmetic diagrams

Two self-contained, deterministic SVG scientific figures; no generated bitmap art, external fonts, scripts, or remote assets.

- `crt-shared-cells.svg`: shared row-major cell labels preserve the actual collision relation. The modulo-5 and modulo-7 colours show equality locally, not root choices. The CRT product is the exact witness (0,9,1;16,15,14;29,21,30). Ordinary row sums are 10,45,80, not equal integer sums.
- `fields-rings-reentry.svg`: distinguishes extension fields from quotient rings at orders 27,81,243; gives the exact all-exponent characteristic-two contrast and smallest doubling re-entry 50 non-P → 100 P → 200 non-P.

Parker means absence of nine pairwise distinct square entries in a 3×3 magic grid. Zero is allowed. Non-P means existence; none of these diagrams claims an integer square solution.

## Replay

Run `python3 -B replay.py --self-test` and `python3 -O -B replay.py`. Only Python's standard library is needed. `generate.py` regenerates SVG and `countdata.json` byte for byte.

Replay independently checks all ring moduli 1 through 100 plus 128,200,243; the CRT witness, roots, masks and all eight modular line sums; and complete normalized field counts at 27,81,243 using independent polynomial arithmetic and irreducibility checks. This also checks minimality of the doubling re-entry, rather than merely showing three examples. Universal all-exponent claims rely on the pinned proofs, not on finite enumeration.

`SOURCE_PINS.json` records SHA-256 hashes of the exact source packets. `sources/` contains internal theorem copies and narrowly selected certificate rows. No external article is copied. Pins of entire source JSON files identify provenance; extracted data has its own packet integrity hash. This packet is portable and does not require its sibling directories.

## Visual QA

Both SVGs were rendered with Inkscape and visually inspected at native resolution. Cell labels, legends, row-sum warning, field/ring entries and re-entry arrows were checked; no clipping or overlap was observed. PNG previews are disposable QA outputs, not source artifacts. No Library or remote publication was performed.

## Verified affine-coordinate interpretation

Index rows and columns by i,j in {0,1,2}. The displayed modulo-5 equality partition A is the fibre partition of i−j modulo 3; the modulo-7 partition B is the fibre partition of i+j+1 modulo 3. Each has three classes, and all nine ordered coordinate pairs occur exactly once. The replay checks these identifications directly against the displayed residues. This is an affine nine-position coordinate interpretation, not the seven-point Fano plane or an asserted tournament map. Colour names may be relabelled; shared cell positions may not.
