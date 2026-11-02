# The next integer exposes an erased edge: 32 to36 to38 states

Status: **PROVED by complete finite computation**, for the explicitly active
source universe and operation language below. This extends the inherited
[fixed-G25 complement carrier](../lower_frontier_complement_fibonacci_20261102/COMPLEMENT_CARRIER.md).
It is not an unrestricted tournament carrier or an all-size recurrence.

## The actual source and the new consumer

Let G25 be the square-sum graph on1..25. Retain all its4,864 orientations
that admit at least one of its20 directed Hamilton paths. The inherited
minimal observer has32 states on this active universe; including the dead
fibre gives33 states, but that enormous dead fibre is NOT enumerated here.
Commands K_P are the actual tournament complement along a currently legal
named path P: retain P's arcs and reverse every other arc. On guard bits
they are x xor(1_G xor E(P)). Active sources remain active after every command.

Adding the literal integer26 creates exactly two new square-sum edges:

    10+26=36,     23+26=49.

A terminal command supplies their two orientation bits and reports the
guarded Hamilton-path count H_G26. All other new tournament arcs can be
fixed identically in both sources; they do not enter this guarded count.
The complete directed G26 bank has24 paths: four with26 at an endpoint,
and twenty with26 internal. An internal26 leaves a two-path spanning forest
when deleted, not an old Hamilton path:10+23=33 is not a square.

**Finite theorem.** Preserving the old32-state observation and the four
possible immediate terminal counts requires36 active-source states.
Preserving every old legality observation and every future count after
an arbitrary legal word of old K commands requires exactly38 active-source
states. One refinement round beyond36 suffices: the counts are36,38,38.

The old active observer fits in five fixed bits; this enriched active
language needs six. This comparison does not claim that the old full
33-state observer ever fit in five bits.

## A one-edge hostile with unchanged sources and ports

Use the lexicographic32-edge order recorded in extension26.json, with bit1
meaning the smaller endpoint points toward the larger. Take guard masks

    x=1620718146,     y=1620750914=x xor2^15.

The sole changed edge is{6,10}. Both sources have exactly the same legal
old path

    (18,7,9,16,20,5,11,25,24,1,8,17,19,6,3,22,14,
     2,23,13,12,4,21,15,10),

and every legal old complement word keeps the same old observer output.
This is an isolated-path state in the inherited compatibility matching.
Use the SAME two new arcs10->26->23 in both. The first source then has
H_G26=1; the second has H_G26=2. Its extra path is

    (18,7,9,16,20,5,11,25,24,12,4,21,15,1,8,17,19,
     6,10,26,23,13,3,22,14,2).

The common path is the old path extended by26. The added path actually uses
6->10, exposing the old ignored bit. Every edge and orientation is checked
both by masks and by literal directed-edge membership. This is a native
source-preserving boundary, not an abstract context with arbitrary new labels.

## Why the computation proves both bounds

The forced leaf18 makes exhaustive depth-first enumeration from18 complete;
reversal supplies all directed paths at both orders. Each old path fixes24
of32 bits. Enumerating its2^8 completions and taking their union gives
exactly4,864 active sources. All four assignments to the two new edges are
examined for every source,19,456 actual source/attachment tests.

Begin with the inherited old state and the vector of four terminal counts.
There are36 distinct outputs. Repeatedly refine by the labelled target
class of every legal K_P. The complete transition table stabilizes at38.
A second initialization uses only current legal path names and the terminal
vector, not the inherited old partition; it reaches the same equivalence.

Stable classes prove sufficiency because outputs and every legal transition
are constant on a class. Distinct classes prove necessity: a first output
difference is already observable, and each refinement split prepends its
distinguishing K command to a distinguishing future observation. This
induction gives a finite distinguishing word, so merging any two classes
would lose a named permitted consumer. Hence38 is minimal for this language.

There are128 old-state/fixed-attachment groups. Exactly four split under the
immediate query, involving four old states; closure accounts for the further
two states. All counts and the explicit witness are in extension26.json.
Normal and optimized Python agree, with explicit exceptions enforcing checks.

## The surviving extension question

The20-object K4 map remains an exact path-level reversal/compatibility map.
Every query on one of those20 paths can of course be translated using its
full lookup table. The failure here concerns a quotient of ORIENTATIONS,
which discarded edges unused by its old path bank. A new path bank can
observe those edges. These are different source types.

Next include old dead sources that become active at26, or allow complements
along newly created26-paths. Neither is included in the38-state theorem.
A useful all-size grammar must retain enough spanning-forest port states
to handle internal new vertices, not only append old Hamilton paths.
