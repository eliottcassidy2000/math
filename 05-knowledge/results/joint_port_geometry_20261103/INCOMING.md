# The same incidence map, two different missing coordinates

**PROVED elementary synthesis / recovered concurrent results as attributed.**
Incoming commits through596ae1ea47 were read after the first computation
checkpoint. No claim below identifies arbitrary abstract covers with actual
tournament metagraphs or transfers unit distances to tournament vertices.

## IM1. Reflection fibres and converse holonomy are dual positions

For a fixed finite graph E with v vertices, e individually identified edges
and c components, use binary cochains and the incidence map

`delta:F2^V -> F2^E, (delta epsilon)_(ij)=epsilon_i+epsilon_j`.

Loops have zero row and parallel edges retain their separate identities.
The elementary exact sequence is

`0 -> H0(E;F2) -> F2^V ->delta F2^E -> H1(E;F2) -> 0`.

An assignment in ker(delta) is constant on each component, so dim ker=c,
rank(delta)=v-c, and dim coker=e-v+c.

In [GP1](GEOMETRIC_PHASES.md), vertex bits are actual reflection choices
of supplied pieces. The observable unit-edge bits are a fixed translate of
delta(epsilon). Therefore a feasible unit-edge pattern has2^c preimages,
or2^(c-1) after simultaneous reflection on a nonempty family. The missing
coordinate lies in **H0**.

The concurrent [Session166 reflection-voltage theorem](../session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md)
instead assigns bits to individual black-edge reflection pairs between
non-self-converse nodes. Changing chosen converse sheets adds delta(epsilon).
Its complete invariant over a fixed labelled quotient is the edge class
modulo im(delta), hence **H1**, with2^(e-v+c) classes. Blue lines and
self-converse nodes retain their separate fixed-locus data. That theorem
classifies equivariant incidence graphs, not tournament realizability.

Thus one construction asks for the fibre of a vertex-to-edge measurement;
the other quotients arbitrary edge data by vertex changes. Exactly,

`(# attainable relative-edge patterns) * (# abstract voltage classes)=2^e`.

The factors are2^(v-c) and2^(e-v+c). It would be wrong to exchange H0 and H1
or infer that either count is an unlabelled tournament count.

For GP3's circle family, c=floor(d/2)+1 and beta=floor((d-1)/2). Its hidden
metric phases number2^floor(d/2), whereas its abstract voltage classes number
2^floor((d-1)/2). The even-d P3 component changes the first exponent and not
the second. This explains the parity split through a component, rather than
a coincident numerical sequence. `geometric_phases.py` checks the binary
incidence rank and product identity through d12.

## Incoming source and geometry boundaries

The concurrent [two-wing/defect packet](../defect_offset_dynamics_20261103/ESCAPES.md)
contains a source with H161 admitting no induced extension of H485 or487,
yet one native path complement followed by a birth realizes each. Its
receipt preserves actual paths and masks. This is the next operation to
consider when a return target lies outside an attachment bank. It complements
R233: changed old edges can escape a fixed source obstruction, but neither
example proves recursive return closure.

[Port-memory P2](../six_roots_hidden_bits_20261103/PORT_MEMORY.md) proves all
n marks of E_n distinguishable by one ear for n>=4. Its logarithmic mark-bit
lower bound is consistent with JE's fixed linear dimension: coordinates in
a Q-linear observer are unbounded integers/rationals, not bounded bit cells.
Preserving a mark and retaining the vertex set remain separate requirements.

The concurrent [unit-diamond compiler](../session166_planar_metagraph_20261103/unit_diamonds/NOTE.md)
classifies injectively realizable tip gluing graphs as linear forests. Its
unit rotation on the radius-sqrt(3) circle has infinite order: its trace5/3
is not an algebraic integer. Reversing direction repeats a point. GP3's unit
circle rotation instead has order6 and yields six-cycles. Finite versus
infinite rotation order is the concrete distinction between the two
composition laws. Both require a supplied radial chart, the actual
unit-chord condition, and injectivity. An abstract diamond or tournament
drawing is not automatically a Euclidean realization.
