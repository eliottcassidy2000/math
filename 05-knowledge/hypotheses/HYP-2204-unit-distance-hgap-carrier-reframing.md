# HYP-2204: Unit-Distance H-Gap Carrier Reframing

**Status:** OPEN, supported by S627 carrier computation and synthesis.

**Claim.** The meaningful relation between the unit-distance tournament model
and the permanent tournament gaps `H=7` and `H=21` is not a literal
unit-distance tournament whose Hamiltonian-path count is `7` or `21`.  It is a
carrier-forgetting test: raw unit-edge counts can visibly equal `7` or `21`
only after spine, tile/bulk, frontier, and exact-vs-lattice side channels have
been collapsed into one scalar.

In other words, `H=7`/`H=21` should be used as forbidden scalar-collapse
certificates.  If a quotient asks a unit-distance carrier to be "the same as"
a tournament `H=7` or `H=21` object, that quotient has forgotten which pieces
make up the visible count.

## Evidence

S627 adds `04-computation/unit_distance_hgap_carrier_s627.py` and stores
`05-knowledge/results/unit_distance_hgap_carrier_s627.out`.

Incoming HYP-2203 supplies the complementary Moser-lane traceability signal:
the unit spine/tie-order is itself a side channel that must be retained.  S627
uses that as one of the carrier packets and asks a different question: what
goes wrong when the spine, tile, frontier, and H-gap channels are collapsed
into a raw scalar?

The edge echo table compares:

- exact small planar values `u(n)` through `n=14`;
- triangular/Harborth lattice lower-bound values;
- the unit-spine decomposition `u = (n-1) + tiles`;
- the S599z literal unit-tournament Hamiltonian-path counts where available.

The two load-bearing rows are:

| Row | Scalar echo | Carrier split | Literal tournament reading |
|---|---:|---|---|
| exact `n=5` | `u(5)=7` | `7 = 4` spine `+ 3` tile/bulk | S599z unit-tournament `H=15`, not `7` |
| Harborth `n=11` | lattice edges `=21` | `21 = 10` spine `+ 11` tile/bulk | exact planar `u(11)=23`; `21` is lattice edge echo, not literal `H=21` |

Thus the unit-distance side does not refute the H-gap theorem.  It explains why
the same visible integer can be legal in an additive edge-count carrier while
remaining impossible as a tournament Hamiltonian-path evaluation.

## Relation To H=7 Impossibility

The H=7 impossibility says there is no tournament whose Hamiltonian-path count
is `7`.  Unit distances have `7` unit edges at `n=5`, but that count is already
decomposed as a Hamiltonian unit spine plus three extra unit tiles.

So S627 treats `u(5)=7` as an equidecomposability signal, not an
equinumerosity theorem.  The scalar is the same, but the pieces and allowed
operations are different:

- tournament `H` is constrained by OCF / strong-component completion;
- unit-distance `u` is an additive edge count with geometry and boundary
  side channels;
- the raw equality `7=7` destroys exactly the information that matters.

This extends HYP-2178/HYP-2181: H-gaps are transferable as coimage
side-channel warnings, not as raw scalar numerology.

## Tournament Analysis

Vertices are carrier packets/proof obligations rather than points:

- unit-spine Hamiltonian path;
- tile/bulk unit-edge packet;
- boundary frontier shell;
- Moser/non-lattice exact carrier;
- OCF H-gap guardrail;
- round-LRC worry-set channel;
- literal unit-tournament `H`;
- raw edge-count quotient;
- equidecomposability ledger.

Pairwise observable: which carrier better preserves traceability, the H-gap
guardrail, and exact-vs-lattice scaling data.

The three single gauges are transitive, but they disagree substantially:

- geometry vs H-gap flips `23/36` edges;
- geometry vs scaling flips `17/36`;
- H-gap vs scaling flips `22/36`.

The majority carrier tournament has score histogram
`{0:1, 1:1, 3:2, 4:1, 5:1, 6:1, 7:2}`, `5` directed 3-cycles, a large SCC
`[unit spine, tile/bulk, frontier, Moser carrier, OCF guardrail, literal H,
equidecomposability ledger]`, and `25` Hamiltonian paths.  This is the useful
LRC-style signal: no single scalar order decides the proof route once geometry,
H-gaps, and scaling are all retained.

## Assumption Challenge

Alternate vertices considered: points, unit pairs, nonunit pairs, distance
classes, Hamiltonian paths, frontier additions, centered hex shell events, OCF
cycle packets, LRC residues, and proof obligations.

Chosen vertices: carrier packets/proof obligations.  This preserves the
predicate "unit-spine traceability plus no forbidden-H scalar collapse" and
destroys continuous embedding families, all optimal graph representatives,
fine nonunit-distance inequalities, and literal `H` values beyond computed
rows.

Challenged assumption: a scalar equality such as `u(5)=7` means the same thing
as `H(T)=7`.

## Next Tests

1. Compute literal unit-tournament `H` for exact Moser/planar carrier rows
   `n=9..14`.
2. Test whether known non-lattice `n>=22` candidates keep a unit-spine
   Hamiltonian path.
3. Define a shared carrier-entropy variable for Sawin's `n^1.014` side and
   the feasible OCF/proof-obligation side, avoiding raw scalar matching.
