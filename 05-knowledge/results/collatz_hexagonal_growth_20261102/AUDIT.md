# Verification and scope audit

Status: same-agent proof review and FINITE-EXACT replay. No independent agent
audit and no Lean build. General theorems use the written proofs, not a finite
search extrapolation. The verifier uses explicit exceptions, not assert.

agents/check_docs.py passed after integration of the concurrent checkpoints:
bounded startup, frontier statuses, correction sentinels and router smoke.

## Reproduction

From the repository root:

```text
python -X utf8 -B 05-knowledge/results/collatz_hexagonal_growth_20261102/verify.py
python -X utf8 -O -B 05-knowledge/results/collatz_hexagonal_growth_20261102/verify.py
python -X utf8 -B 05-knowledge/results/collatz_hexagonal_growth_20261102/seal.py --check
```

The first two stdout receipts are stored in verification.json and
verification.optimized.json and must be byte-identical. The checker also
rewrites the source-labelled TSV deterministically. seal.py checks LF-normalized
text hashes and literal binary hashes, including the imported hostile source.
The optional plot.py uses Matplotlib; it is not a verification dependency.
The PNG and SVG were rendered and the PNG visually checked. Matplotlib was
installed in an isolated temporary environment, leaving repo dependencies intact.

## Declared universes and independent paths

| Check | Universe | Expected evidence |
|---|---|---|
| Native inverse census | All ordinary C first-hit predecessors of5 through26 steps, no source filter |2263 labels; every edge and full forward path independently replayed |
| Lattice capacities | H and triangular balls of all radii0..105 | Coordinate BFS agrees with the separately proved formulas |
| SF5 | Selected two-branch odd-unit family through13 odd levels |16383 distinct labels; raw clocks replayed; no all-integer coverage claim |
| Joint native meetings | Every unordered pair of511 SF5 nodes through depth8 |130305 pairs; literal path intersection agrees with J4 |
| Weighted native union | All nodes/first-hit clocks in expanded SF5 depth0..8 | All coefficients agree with J5; depth8 has2001 states, not2811 |
| Gram controls | Six legal sibling exponents3,5,9,11,15,17;256 depth8 leaves | Literal tail-edge overlaps; exact rational elimination agrees with J7; private edges certify leaf rank256 |
| Periodic clock hostile | Sources13/853, fixed port53 | Same exponent mod6, different meetings40/160 |
| Residue precision hostile | k2..8, roots19 and19+2*3^k, exponent2 | Same legal branch and root residue, different successor residue |
| Minimal chart | All unit source residues modulo3^(h+1), h1..7 |2*3^h distinct lower-branch menu traces and weighted polynomials; totals lose information |
| Independent chart replay | All240 root classes at h1..4, with actual odd lifts>1 | Polynomial coefficients equal literal native unions; one merged STOP state |
| Prime multiplier | All odd n<=201 not divisible by5; all exponents3..22 | Exactly four increasing unit predecessors each; separate13/33/83 cycle replay |
| Palette transport | All48 maps x->Px+v, P in S3,v in F2^3 | Literal fixed frames versus kernel/image formula; histogram0:33,2:8,4:6,8:1 |
| Functional colouring | Every total map on each labelled set of size1..5,3413 maps | Literal3-colourings agree with tree/unicycle formula |
|81/27 hostile | Original twelve-vertex pair in regular_parent_symmetry_20261029/verify.py | Same67 cyclic triples; Aut81/27; regular binary masks1620/1728 |

The last check imports only hostile_checks from the pinned original verifier;
it does not claim to independently reproduce all tournament computations in
that packet. All other arithmetic and combinatorial checks are implemented
locally, with literal forward replay or exhaustive finite enumeration as
described. A duplicate computation path is not an independent human audit.

## Direction, boundary and type audit

- Lattice obstruction assumes bounded edge image length and bounded fibre
  size. It does not rule out arbitrary curves, shrinking scales, or an
  unseparated unit-distance drawing. The ordinary1/4/2 triangle is a
  separate honeycomb obstruction; the root5 certificate is acyclic.
- Radius21 and26 are the first failures of this root's ball-count test,
  not minimal graph obstructions over all roots and layouts.
- The dyadic drawing is of macro edges. Shared even tails must be glued;
  a planar macro drawing alone proves no claim about that particular
  drawing after identification. J1 describes the abstract native tree.
- J4 is an asynchronous first meeting, retaining both clocks. It is not
  a synchronous equivalence statement and does not conflate U with C.
- J5 counts every expanded source state once; children subtrees are
  disjoint by odd parts. Clocks are ordinary steps; exponent a costs a+1.
- General p branching assumes2 primitive modulo p^2. This is checked at5;
  no unproved assertion for every prime is imported.
- Chromatic number3 concerns the simple arithmetic arrow graph, not the
  Euclidean plane or a graph containing all unit-distance coincidences.
- Functional palette sections use a different graph carrier, retaining
  one cell per arrow, including loops and opposite arrows. This avoids
  deleting cycle holonomy while simplifying the colouring graph.
- Gram rank counts independent port responses; it is not a lower bound
  on every possible algorithm or a proof that scalar formulas cannot exist.
- Source7 and the5n+1 cycle block interpreting a selected inverse family,
  a decreasing rank on that family, or compatible transport as universal
  Collatz convergence. No statistical average is used to exclude integers.

## Corrections caught before release

An initial phase probe used parent7 and exponent2; its child9 is divisible
by3 and is not legal SF5. The explicit guard rejected it. The corrected
probe19/37 uses legal odd units on both sides. No false claim was promoted.

The initial macro-drawing prose did not state the shared-even-tail gluing
obligation. Reading the concurrent root85 proof exposed that gap;13/53
supplied the exact witness. JOINT_PORTS repairs it and the global mistakes
ledger preserves the failed inference. The antecedent SF5 theorem remains
valid; it concerns the odd inverse tree with clocks.

A later incoming dictionary repeated THM-4852's already retracted generic
one-primorial-digit selector claim. The table and CA comparison were repaired
to retain denominator valuations, with the existing2/8 modulo6 hostile.
Its README hash was refreshed; all simulation receipts were preserved.

## External claims

Primary sources inspected:

- Applegate--Lagarias, [Density Bounds for the3x+1 Problem I](https://websites.umich.edu/~lagarias/doc/applegateI.pdf),
  Theorem1.1. Classical quantitative inverse growth is cited for context;
  the geometric proof uses the repository's elementary SF5 construction.
- OpenAI, [The Euclidean plane is not five-colorable, pinned introduction](https://raw.githubusercontent.com/openai/math/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/build/sections/introduction.tex),
  Theorems1.1 and1.3. This is CITED external status. No complete independent
  proof review, Lean build, or axiom audit is represented as having occurred.
- Kuehn, [Minor-universal planar graphs without accumulation points](https://web.mat.bham.ac.uk/D.Kuehn/unvacc3.pdf),
  Theorem2, especially condition(iii) for nested cycles and disconnected
  graphs. This is the CITED primary input to Session140's inherited
  honeycomb subdivision proof; no effective-layout claim is added.

All results newly proved here are independent of the external six-colour
lower bound. No downloaded papers or full extracted source texts are
redistributed. Internal provenance and frozen input hashes are separate
from literature-priority claims.
