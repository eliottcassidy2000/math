# Count corridors, reframing, and a non-universal route to coverage

**Status: PROVED elementary criterion; FINITE-EXACT experiment through order7.
All-order continuation is OPEN.**

A framed tournament has the fixed directed path0→1→...→n−1. Its free
coordinates are all pairs(i,j) with j≥i+2. Consider the following moves:
(a) flip one or two free arcs simultaneously;
(b) switch a cut of the free-edge graph K_n minus that path;
(c) complement all free arcs;
(d) reframe/relabel the SAME unframed tournament by another Hamilton path.
Only retain an edge when its two counts differ by0 or2. Reframing has zero
count change. No claim is made about intermediate states in simultaneous
moves. Moves(a–c) keep the named Hamilton path, so every endpoint is physical.
These are Session274's chart/gauge objects with a new scalar transition filter.

## Integer intermediate value, with an actual initial state

**Lemma.** On any connected graph of odd integer-valued objects whose edge
changes have absolute value≤2, the attained values are an odd integer interval.
A path joining values a≤b must hit every odd value in[a,b]: at the first
vertex whose value is≥c, its preceding value is≤c−2, hence the new value is c.
This proves the assertion, including edges of change0 and nonmonotone paths.

Let C be the component containing the six-vertex tournament with outgoing
rows(14,60,56,48,33,1), H=23. Enlarge the graph by neutral source or sink
padding and their inverses, and all the moves above at every finite order.
**Sufficient coverage criterion:** if H is unbounded on C, then every odd
integer≥23 is realized. An arbitrary finite path to an object of H≥target
is an exact selector by the first-crossing argument. The initial value23
and the continuation obligation are separate. Padding preserves H and is
physically explicit; the criterion does not assume that a gauge commutes
with padding. No arbitrary subtraction, formal negative module, or limit
object is needed.

The criterion is not asserted necessary for Chan–Pak. It asks for one
unbounded count corridor, not a grammar reaching every tournament, and not a
preselected formula for each integer. It could be proved by a growing sequence
of overlapping corridors or a recurrent component escape construction.
Unbounded SIZE in C is immediate by padding; unbounded H is the OPEN part.

## Exact structural experiment

`experiment.py` enumerates ALL fixed-path charts for n=3..7, respectively
2,8,64,1024,32768 objects. Vertex-subset dynamic programming computes H.
Canonicalization ranges over all permutations within equal-outdegree blocks,
with the blocks ordered by degree. This is a complete isomorphism test:
any isomorphism preserves degree and hence lies among those candidates.
It recovers2,4,12,56,456 ordinary classes. Canonical-identical framed objects
are joined neutrally only for the `reframed_*` tests.

The full move set has exactly as many components as maximal consecutive-odd
runs in the actual spectrum at each tested order:1,1,2,5,10. Thus each such
component contains every framed representative in its run, not just one witness
per value. No global count-spectrum claim is extracted from this census.

The effect of restoring the frame choice is substantial. At order7, fixed-chart
full moves have111 components; allowing reframing leaves10. Reframed single
arc moves leave159 components; allowing pairs leaves17; cut gauges and full
complement together with those pairs leave10. Therefore ordinary single-flip
Morse reasoning loses genuine bridges even before the all-order question.

At order7 the component containing23 spans exactly23..61. Another spans65..105.
The value63 is absent at this order; a bridge between them necessarily needs
new vertices. This identifies a concrete first stabilization problem rather
than another scalar range forecast. At order6 the corresponding component
spans23..33; vertex growth can enlarge corridors but the all-order rule remains
unproved. Small source-padding orbits alone should not be mistaken for such a rule.

## Provenance and scope

Historical outputs `morse_minima_proof_89.out`, `morse_reverse_proof_89.out`,
`homological_morse_23.out` inspired studying connected count levels rather than
only individual witnesses. Their global proposed minima statements are not
used. Session274 supplies the exact framed/gauge semantics. The elementary
intermediate-value lemma is not claimed as a new general mathematical result;
the contribution is the exact physically specified corridor target, its initial
state, and the operation-sensitive experiment exposing which frame data matter.

Reproduce from repo root:
`python3 05-knowledge/results/session275_exact_selection_20261111/landscape/experiment.py`
The output must equal RESULTS.json. The use of an assertion to check interval
connectivity is only a consistency check, not part of computing components;
an independent path-based review is requested before final publication.
