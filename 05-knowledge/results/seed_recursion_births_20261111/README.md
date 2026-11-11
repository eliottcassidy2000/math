# From stable seed choices to recursive families and birth obstructions

**11 November 2026. PROVED IN NOTE / FINITE-EXACT with stated inherited
inputs. Universal odd-c seeds and Chan–Pak remain OPEN.**

The 95 strong-core stable-port presentations and 49 stable pairs extend
to 60 full commuting inflation families, with up to five independent
block sizes. All remain unit seeds for every positive size tuple.
The true marked-isomorphism quotient has 58 such core types.

Seven missing head families now have explicit Hamiltonian-path
constructors. In particular, for every r,s>=1, we construct tournaments
with **H=1+(3*2^r+1)*2^s** and **H=1+(21*2^r+1)*2^s**. This closes all
below-ray gaps on these two infinite coefficient families.
Equivalently, in the unique indexing H=1+2^s+d*2^(r+s), the entire
odd-d grids d=3 and d=21 are now filled.

Composition is still constrained. All 322 old-vertex transport tests
are classified at every positive inflation size: a single inflation
never activates a previously unstable old vertex. One specific
Hamilton-doubling interface fails after changing its source. A finite
fixed bank of these templates has only O((log X)^6) Hamilton values below
X, so universal coverage needs a justified change of template.

The marked seed birth graph has 4,7,43 deletion-root types at orders
5,6,7. Its full marked-orbit sequence is 1,2,4,9,22,59,194 through seven;
only 1,2,4,9,18,32,53 of these can descend to the singleton by unit-seed
deletions. Keeping the actual source distinguishes roots with c=1 from
the singleton with that same coefficient.

- [Proofs and exact scope](PROOFS.md): R1–R9, full tensors, the completed
  lines, global fixed-kernel exclusions using the inherited H<=375
  inventory, source-compatible recursion, and new marked sequences.
- [Arithmetic birth laws](ARITHMETIC.md): summand/multiplicand fibres,
  doubles/squares, exact strict generators, and a complete all-size
  square-product graph classification beside the square-sum graph.
- [Verification and reproduction](AUDIT.md): universes, inherited filters,
  independent counting checks, source pins, and commands.
- [Research board and next decisive steps](BOARD.md): unresolved heads
  and the operation needed to escape the fixed-bank boundary.

An explicit constructor, from the repository root:

```text
python3 05-knowledge/results/seed_recursion_births_20261111/realize.py --family three --parameter 3 --exponent 3
```

This emits an eight-vertex tournament with H=201. The supported complete
line families are `three` and `twentyone`. `mixed11` supports exponent 1
and all exponents >=3; it explicitly rejects the unresolved exponent 2.
The adjacency convention is bit j of row i means i -> j.

These results refine A000568's marked tournament structures and connect
with the existing A135244/A002496 inverse work; no new general enumeration
record or new term of those three sequences is asserted.
