# Complete descendant fibres at Hamilton counts 31 95 and 97

**FINITE-EXACT complete bounded-value certificates, with all word lengths
included.** The operation language is the incoming triangle-seed two-live-port
language, not arbitrary ears, arbitrary source graphs, or free re-marking.

## L1 Operations and universe

Start with the directed triangle and two named slots. A command chooses
one slot as root, adds a vertex x with root->x and x->every other old
vertex, then forgets one of the three eligible names. Forgetting removes
eligibility for future root choices, not the vertex. Encode a command by
3*root_slot+discard; discard=2 forgets the newborn. There are six commands.

[LIVE_PORTS L1-L4](../adaptive_state_induction_20261103/LIVE_PORTS.md)
proves a positive 24-coordinate path-cover carrier, sharp rational observer
rank 17, and strict H growth. Its previous finite saturation realizes every
allowed odd value through 999. Those are inherited results.

The current experiment saturates the entire 24-coordinate orbit at H<=295.
It has 30,087 states and 39,300 transitions that stay within this bound.
Every state has a decoded literal source, of order at most 11. A separate
C++ algorithm recounts H on every one of these 30,087 actual graphs by
subset path DP. It also independently checks every parent certificate,
every matrix transition, closure of the proposed set, and the descendant
support tables below.

This is all-length completeness at a value bound: every actual command
strictly raises H, so a path that leaves the bound cannot return. No
vertex-size limit or heuristic profile filter is used in generation.
Equal carrier states have the same future counts by the inherited exact
transfer theorem; merging them for this observation does not license
replacing their graphs for arbitrary operations.

## L2 Full target supports

For a state with H=h let bit j indicate existence of a descendant with
H=3h+2j, j=0,1,2, reached by any positive number of these commands. The
following histograms count all reachable carrier states in each fibre.

| Descendant mask | H31 fibre | H95 fibre | H97 fibre |
|---:|---:|---:|---:|
|0: none|6|41|122|
|1: 3H only|0|37|84|
|2: 3H+2 only|3|11|8|
|3: 3H and 3H+2|0|10|0|
|4: 3H+4 only|10|8|8|
|5: 3H and 3H+4|18|2|2|
|6: both positive offsets|0|1|0|
|7: all three|0|0|0|
|Total|37|110|224|

The computation propagates target bits backward through the entire finite
DAG in descending H order. It is not limited to one or two births.

**Consequences within this language:**

1. No H31 source state has both positive-offset descendants, even though
   all three scalar successors 93,95,97 occur elsewhere in the language.
2. No H97 state has both positive offsets. Exactly one H95 state has both
   287 and 289 as descendants, and it lacks 285.
3. No source in any of these three fibres supports all three branches.

The receipt includes an actual source program, graph, retained pair, and
descendant programs for every nonempty histogram category. These positive
witnesses are decoded and independently counted, not merely inferred from
the support bits.

## L3 The quantifier boundary and the next construction

Scalar coverage through 999 says that each target has *some* source.
A closed induction needs one current source with all required continuations,
whose chosen children again satisfy the induction interface. The table
separates these statements at H31, long before a large-value test.

The result does not rule out a different seed collection that avoids H31,
state-dependent radix rules, additional retained ports, or actual native
resets between supplied sources. In particular, the earlier four-level
ternary construction uses a larger source alphabet and actual births,
complements and deletions; it is not contradicted by this append-only result.
The [incoming native excursions](../adaptive_state_induction_20261103/METAGRAPH_RESETS.md)
remain candidate operations, with their source and endpoint guards intact.

This is a stopping reason for the specific proposal to treat the two-port
coverage experiment as an already ternary-closed source bank. No global
Hamilton-path-count omission besides the inherited 7 and 21 is inferred.
