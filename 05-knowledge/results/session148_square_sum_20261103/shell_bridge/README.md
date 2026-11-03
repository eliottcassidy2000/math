# A cardinality bridge and its exact graph obstruction

## Objects and map

Let G18 have vertices 1,...,18, with distinct x,y adjacent iff x+y is a perfect square. Delete vertices 13,14,15 and their incident edges. Its three components are:

- 6 vertices: {1,3,6,8,10,17}; edges form the path 10–6–3–1–8–17.
- 4 vertices: {4,5,11,12}; edges form the path 12–4–5–11.
- 5 vertices: {2,7,9,16,18}; edges 2–7, 7–9, 7–18, 9–16.

On the other side, let prime shells be π(Fi22) = {2,3,5,7,11,13} (size 6), π(Baby Monster)\π(Fi22) = {17,19,23,31,47} (size 5), and π(Monster)\π(Baby Monster) = {29,41,59,71} (size 4). These follow from [Fi22](https://brauer.maths.qmul.ac.uk/Atlas/v3/group/Fi22), [Baby Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B), and [Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M) order factorizations. The actual groups need not themselves form an embedding chain for these set differences to be defined.

Match the size-6, size-4 and size-5 source components respectively with these 6, 4 and 5 prime shells. Pair each source component's increasing sequence to its prime shell's increasing sequence. This imposes the bijection:

- size 6: 1→2, 3→3, 6→5, 8→7, 10→11, 17→13;
- size 4: 4→29, 5→41, 11→59, 12→71;
- size 5: 2→17, 7→19, 9→23, 16→31, 18→47.

The size match is exact. Increasing-rank matching is unique after that convention is imposed, but it is not a canonical mathematical identification from cardinality alone. There are 6!·4!·5!=2,073,600 shell-respecting bijections before selecting one.

## Exact square-sum comparison

On the prime side, use the simple induced square-sum graph on all 15 shell primes, with adjacency also defined by distinct values whose sum is a square. Under the rank map, among the 105 unordered source pairs:

- 5 are edges on both sides;
- 7 are source edges but not prime-side edges;
- 9 are prime-side edges but not source edges;
- 84 are nonedges on both sides.

For example, 1–3 is a source edge because 1+3=4, but its image 2–3 is a nonedge because 2+3=5. Conversely 1 and 9 are not adjacent (sum 10), while their images 2 and 23 are adjacent (sum 25). All preserved/lost/created pairs are in `results.txt`.

The source graph has 12 edges and three connected components. The prime-side graph has 14 edges and is connected. Thus no graph isomorphism exists under any bijection whatsoever, not merely the chosen rank map. Also no isomorphism exists within any matched component/shell: edge counts at sizes 6, 4, 5 are respectively 5,3,4 on the source and 3,2,2 on the prime side.

The 14 prime-side edges split into seven within-shell and seven cross-shell edges:

- Within Fi22 shell: 2–7, 3–13, 5–11.
- Within Baby-minus-Fi22 shell: 17–19, 17–47.
- Within Monster-minus-Baby shell: 29–71, 41–59.
- Across shells: 2–23, 2–47, 5–31, 5–59, 7–29, 13–23, 23–41.

These are exact finite checks. The equal cardinalities provide a useful comparison and do not by themselves imply a moonshine, group-action, or graph-theoretic correspondence. A different target graph or added structure would require its own definition and proof; the present obstruction concerns precisely the two square-sum graphs above.

## Replay

Run `python check.py` and `python -O check.py`; both outputs equal `results.txt`. The check uses only the standard library, constructs every unordered pair directly, verifies every stated count, and uses explicit exceptions rather than `assert`, so optimization does not bypass its checks. It is read-only.

## Source pins

Prior repository snapshot commit `74e81ed37c9175d94e8fad57fc5809cf78cee572`: Session132 `prime_graphs/PROOFS.md` Git blob `79ca61cacb8e1c9ad6a4b07077fd76a94c76d2cc`; Session131 `moonshine/README.md` Git blob `dd14df76a81c029e95d41715fa82f2ce2d264944`. These source pins identify context only; this packet's finite graph check does not import them or the earlier prime packet. The group-order data comes from the linked ATLAS pages.
