# Two-edge escapes and source-preserving Boolean inverses

**2026-11-02. PROVED scoped lemmas / FINITE-EXACT classifications /
CITED inputs. Same-agent audit with independent computational paths.**

The768 target17 traps split into480 with three two-edge escapes each
and288 with none. All768 have sharp three-step escapes of maximum squared
error36; no route has a lower maximum. A runnable six-vertex Boolean
decoder now finds H17 within three edge flips from every orientation.

The deeper obstruction is symmetry. A deterministic relabelling-equivariant
inverse must preserve the source stabilizer. The trap groups have orders3
and5, so an H17 inverse needs a distinguishing port or external choice.
Their smallest invariant target-supported laws have3 and5 atoms.
The actual81/27 hostile similarly has sharp minimum equivariant
target-support sizes3 and1 for H27.

At order seven the smallest observer for C plus the17/71/33 re-marking
has10,455 states; adding converse requires16,149. All named legal
single-vertex relocations force all32,768 orientations. Even legality
alone separates them. This full separation is verified through seven,
not proved at every size.

The G25-to-G26 birth calculation independently recovers9,184 active new
orientations from inactive old sources. It gives an exact two-port
distribution and a conditional-energy decoder that preserves all old bits
while creating exactly one square-sum-compatible Hamiltonian path.

Read:

- [Escapes, symmetry and carrier minima](ESCAPES_AND_CARRIERS.md).
- [Two-port birth formulas and actual constructed sources](BIRTH_AND_MOMENTS.md).
- [Magic squares, the nine-cycle, lambda cube, octonions and the pasted-source audit](ARITHMETIC_AND_SOURCE_AUDIT.md).
- [Validation and inheritance](AUDIT.md), [next decisive tests](HANDOFF.md).

## Reproduce

From the repository root:

    python -B 05-knowledge/results/two_edge_transport_20261102/scout.py
    python -B 05-knowledge/results/two_edge_transport_20261102/escape_depth.py
    python -B 05-knowledge/results/two_edge_transport_20261102/carrier_scout.py
    python -B 05-knowledge/results/two_edge_transport_20261102/birth_scout.py
    python -B 05-knowledge/results/two_edge_transport_20261102/verify.py
    python -O -B 05-knowledge/results/two_edge_transport_20261102/verify.py --output 05-knowledge/results/two_edge_transport_20261102/verification_optimized.json
    python -B 05-knowledge/results/two_edge_transport_20261102/seal.py --check

Construct an actual target17 tournament, retaining every intermediate step:

    python -B 05-knowledge/results/two_edge_transport_20261102/decode17.py 140
    python -B 05-knowledge/results/two_edge_transport_20261102/decode17.py 140 --max-error 36

All code uses the standard library. The exact transcript and normalized
source hashes accompany the files. The ordinary all-order spectrum,
all-order carrier injectivity, magic-square-of-squares existence and
the broader conjectures remain open. The pasted Hecke table is not a
defined eigenform calculation; its weight-two cuspidal interpretation
is ruled out by genus zero.

