# An exact parent-only refinement prediction of A000568(8)

24 October 2026. **PROVED:** the conditional unpointing/error mechanism and its
order-eight specialization using the cited holonomy theorem. **FINITE-EXACT:**
the supplied parent atlas, computed refinement sum, and resulting integer.
**OPEN:** an unrestricted all-order predictive correction or size-only
recurrence. No novelty, literature-priority, or performance claim is made.

## 1. Inheritance and the retained information

The frozen Session50 growth packet supplies the parent attachment-subset
quotient and exact unpointing trace. Its executable also knows how to recover
root orbits through a deletion/re-rooting kernel. We retain its complete
parent representatives and automorphism permutation actions, but deliberately
remove that kernel from the prediction stage.

The separate Session50 cell-holonomy packet supplies the structural theorem
that every strong tournament on eight vertices has stable refinement cells
equal to its automorphism orbits, as well as the SCC error-series identity.
The signed triangle-cell mechanism and sixteen-class matching family have
Session39/40 antecedents, as documented in that packet. The contribution here
is a separately executable, restricted-input application, not a claim to have
invented those antecedents.

`source-pins.json` points to literal archived copies of the growth source,
its theorem/validation/manifest, and the frozen holonomy theorem/manifest.
The old executable is present for input generation and input validation only;
`predictor.py` neither imports nor reads it.

## 2. Input contract and pointed states

For each n from zero through seven, the input lists one representative T of
every ordinary n-vertex tournament class, with its complete orientation-
preserving permutation group G=Aut(T). Rows are outneighbor bitsets, and an
action g is an old-vertex to new-vertex map. Completeness of the atlas and of
each action list is part of the contract, not inferred from group closure.

For an order-seven parent T, let S run through one representative of every
G-orbit on subsets of its vertices. Form U(T,S) by adjoining a new vertex *
and orienting v to * exactly when v belongs to S. These pairs are in bijection
with ordinary isomorphism classes of pointed order-eight tournaments:
removing the distinguished vertex recovers a parent, and an isomorphism fixing
the distinguished vertex restricts precisely to a parent isomorphism carrying
the attachment subset. Thus completeness and nonduplication require only the
supplied parent atlas and its actions, not any target-order comparison.

Let C_p be the stable directed color cell of the new vertex in the extension
for pointed state p. Stable refinement starts from one color, equivalently
from outdegrees after its first step, and repeatedly records the old color
and outgoing counts into each current cell. Incoming counts are determined
in a tournament. Define

    U_N = sum_(pointed N-vertex classes p) |C_p|.

The blind program computes U_8 from the parent subset action and these local
refinement calculations alone.

## 3. Exact error accounting

For an unpointed tournament U, its pointed classes correspond to Aut(U)-vertex
orbits O. Let C(O) be the stable cell containing O, and define

    e(U) = sum_O (|C(O)| - |O|).

Since stable cells are invariant, this is a nonnegative integer. Summing
|O| over the vertex orbits gives |U|, so, with a_N the number of ordinary
unpointed classes and E_N the total e(U) over those classes,

    U_N = N a_N + E_N.

This identity holds at every order. By itself it is not a predictive scalar
recurrence: E_N retains missing target-order symmetry information.

Write A(z)=sum a_n z^n, a_0=1, and E_s(z) for the error series over strong
classes, including singletons as strong. The unique ordered strong-component
decomposition preserves the stable cell partition within each component and
makes the full automorphism group their direct product. Hence e is additive,
and marking the component contributing the error proves

    E(z) = A(z)^2 E_s(z).

This is ordinary class enumeration of ordered strong-component sequences,
not a Burnside fixed-object sum and not a labelled generating function.

## 4. Why the order-eight correction is known without order-eight tests

The archived holonomy theorem proves:

1. Every tournament through order six has e=0.
2. Every nonregular tournament on seven vertices has e=0.
3. Every strong tournament on eight vertices has e=0.

The proof uses regular odd stable cells, uniform crossings between coprime
cell sizes, transitivity of regular five-vertex tournaments, and the fact
that at most two triangle cells have no cycle obstruction. A seven-cell at
order eight leaves a uniformly adjacent singleton, making that tournament
nonstrong. This is a structural argument independent of the order-eight
census or the count sought here.

The blind program computes E_s,7 from the known parent atlas, not from a
hard-coded exception count. It refines each known parent and forms its vertex
orbits from the supplied actions. Exactly one seven-vertex parent has positive
error: rows (74,84,97,22,37,11,56), with Aut order 3, orbit sizes 1,3,3, and
one stable cell of size 7. It is strong and contributes

    (7-1) + (7-3) + (7-3) = 14.

All other input-parent errors vanish. With E_s,8=0 supplied structurally,
the series identity yields

    E_8 = [z^8] A(z)^2 E_s(z)
        = (a_0 a_1 + a_1 a_0) E_s,7 = 28.

No order-eight canonical representative, order-eight automorphism action,
child deletion, root-changing kernel, or target count is fed into this step.

## 5. Exact blind computation

The 456 order-seven parents yield 54,256 attachment-subset orbits. The stable
cell of the new root has size 1 in 53,882 states, size 3 in 352, size 5 in 12,
and size 7 in 10. Therefore

    U_8 = 53882 + 3*352 + 5*12 + 7*10 = 55068,
    a_8 = (U_8 - E_8)/8 = (55068 - 28)/8 = 6880.

These numbers are emitted results and separate controls, not values used to
construct U_8 or the correction. The numeric answer/trace literals appearing
in the source boundary's blacklist are checked against injection and do not
enter the prediction arithmetic. The runtime reads only its own source and
the pinned input atlas; deleting all archived sources and result receipts
still leaves the same blind result.

## 6. Scope and failure boundary

This is an exact finite prediction one order above a complete input atlas.
It eliminates target-order canonicalization and deletion/re-rooting
isomorphism tests, but it does not eliminate the parent isomorphism structure,
its true group actions, or the labelled adjacency needed for color refinement.
It does not compute from n alone or from an unmarked blue/black line graph.

At order nine, signed interactions among three triangle cells can be
unbalanced and force rigidity while refinement retains three cells of size
three. The archived theorem proves E_s,9 at least 72 from an explicit family.
Thus its order-eight zero is not a persistent vanishing rule. The executable
supports only target order eight and rejects nine rather than applying an
unsupported correction. Any independent higher-order extension must compute
or prove its missing strong error separately.
