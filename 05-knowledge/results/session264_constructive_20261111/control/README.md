# Control lane: physical modular targets, not just free source words

**New PROVED mechanisms:**

- A legal four-vertex sign flip and one-sink padding generate elementary
  residue shears. Every slope residue is constructible at either exact
  unit offset modulo every n coprime to30, using at most8n²-2 source vertices.
- The signed full-state observable z=r(M-4(h+r)+E) decreases exactly2
  under F5. The paired operation keeps r and decreases z by4.
- That hidden affine clock becomes a pure full-state translation modulo
  every prime p other than2,5,761 after p²-1 applications. Conjugates span
  all three state coordinates. Every target (M,h+r,E) and either exact
  offset sign are therefore reachable by positive physical F3/F5 words.
- A complete250-state finite certificate supplies the omitted p=5 case.
  Thus full-state modular controllability holds for every odd prime except
  possibly761. The prime761 remains an exception of this proof only.

The inherited Session199 construction already gives scalar congruence
density for all moduli in another architecture. The new contribution is
exact-offset moving-mark control, especially all three source coordinates,
not a new global scalar-density claim. No integer target-height bound,
prime-power lift, or Chan–Pak conclusion follows.

## Artifacts

- `PROOF.md`: complete new arguments, quantifiers and boundaries.
- `SOURCES.md`, `SOURCE_PINS.json`: diverse actual old-topic archaeology and
  typed transfers, including inherited profinite and binary-append hostiles.
- `controls.py`: exact tests and executable positive-word constructors.
- `RESULTS.json`: exact receipts, including literal small graph recounts.
- `MOD5_TREE.json`: all250 residue states and parent edges for the p=5 proof.

Run the new controls:

    python3 -B 05-knowledge/results/session264_constructive_20261111/control/controls.py

Construct a literal source with r=-1, M=3 mod7, q=1 mod7:

    python3 -B 05-knowledge/results/session264_constructive_20261111/control/controls.py --construct-slope 7 3 -1

Construct a literal source with r=+1 and M=q=E=0 mod3:

    python3 -B 05-knowledge/results/session264_constructive_20261111/control/controls.py --construct-full 3 1 0 0 0

Both commands output the actual outgoing adjacency lists, selected vertex
indices, positive operation word and exact integer source state. Vertices
use zero-based labels. A graph is supplied, not just residue matrices.
The full-controller constructor uses O(p³) residue states, so larger primes
are theorem-level targets rather than a claim that the implementation scales
comfortably to arbitrary primes.

The controls pass38 literal graph/closing checks,255 mixed-word clock checks,
24 quotient orbits,15 prime translation identities, and261,630 full residue
states across all odd moduli3 through31. The composite full-state cases are
FINITE-EXACT observations only.
