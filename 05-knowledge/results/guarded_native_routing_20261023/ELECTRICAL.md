# Actual electrical guards: window-local routes and retained target owners

23 October 2026. **FINITE-EXACT** on the four declared original order-eight
source fibres. The all-size theorem in PROOF.md concerns count-residue
guards, and does not authenticate an electrical factorization. Here the
original potential and native owner map are used directly.

## 1. The move menu and its relation to incoming work

The incoming [native guarded-move packet](../native_guard_moves_farey_20261023/native_moves/NATIVE_GUARDED_MOVES.md)
already enumerates these four full Hamiltonian fibres and supplies a local
strong-four transposition, selected-port rectangles, and the important
distinct-port capacity obstruction. We read its proof and correction on
arrival. We do not claim our overlapping path census or existence of a
non-automorphism move as a separate discovery.

Its cost is the NUMBER of changed owner clocks, allowing distant changes.
Our different menu replaces a CONTIGUOUS block of at most b internal clock
positions. Every intermediate state must be a Hamiltonian path in the same
actual source, with the same endpoint owners and original electrical phase.
We also test the stricter version that fixes the exact merged target owner.
For a supplied block, use uniform conditional resampling over accepted
paths agreeing at every outside clock. This is a reversible PSD heat-bath
chain on the actual finite fibre.

Both menus are atomic observer reorderings. Neither is asserted to be an
old elementary metagraph edge or a globally balanced wiring switch.
Window width, changed-owner count, cell-phase support, event length and
electrical work are different costs. In particular a two-owner strong-four
swap spanning four clocks has window width four, not two.

## 2. Complete finite universe and sharp window requirements

We retain all 1,866 paths, with no source symmetry quotient. The four source
automorphism groups are trivial in the incoming complete census. Partition
first by (first owner,last owner,original phase); separately refine by the
exact merged target owner. Every nonempty fibre is included. For each one,
test every window width from two through six and every accepted block
transition. The least connected width is exact: connectivity supplies an
upper certificate, and the smaller-width components supply the lower one.
Singletons have cost zero.

| Source | Endpoint/phase fibres | Nontrivial | Require width six | Nontrivial after fixing target | Target-fixed fibres requiring six |
|---:|---:|---:|---:|---:|---:|
|1238|221|98|44|46|19|
|2242|219|99|41|42|17|
|2246|283|142|67|29|9|
|2315|291|146|64|38|11|

The 216 width-six endpoint/phase fibres are not a contradiction to the
incoming 105 fibres needing SIX CHANGED OWNERS: contiguity is a stronger
restriction. A short sparse reordering can need a longer containing window.
Both retain the actual clocks. Allowing width six always gives a complete
refresh inside a two-endpoint fibre, so that trivial upper bound alone is
not evidence of locality; the lower-budget disconnections are essential.

For every nontrivial fibre we also construct shortest legal routes and
calculate their exact rational all-pair, length-weighted congestion. At the
minimum connected width, the worst values by source are respectively

    78/5, 162/7, 18, 28.

Thus this finite collection has a native window-chain gap at least 1/28.
With the target owner also fixed, the respective worst bounds are

    6, 4, 4, 4.

These are guaranteed gap LOWER bounds from the declared routes, not exact
eigenvalues or evidence for the same constant at larger orders. All costs
use the actual fibre measure, not uniform target classes. Singleton fibres
have no nonconstant spectrum. The maximum shortest-route event length in
the first collection is six and in the target-fixed collection is two.

The subsequently received [universal electrical guard theorem](../session44_guarded_native_moves_20261023/electrical-duality/UNIVERSAL_ELECTRICAL_GUARDS.md)
gives this second test a sharper interpretation. On the authenticated full
order-eight atlas, at a fixed source, agreement of ALL original modulo-eight
harmonic characters is equivalent to equality of the merged target owner.
Thus our target-fixed window chains retain that entire bank. This equivalence
uses the incoming finite native target-recovery result; it is not inferred
from one selected potential or extended to every order.

## 3. A small move that retains the source AND target

In source 1238, the actual paths

    P=(0,4,1,5,3,2,6,7),
    Q=(0,4,3,1,5,2,6,7)

differ only at clocks 2,3,4. Both retain endpoint owners 0 and 7, merged
target 3442, and the original electrical phase 4. Every actual Hamiltonian
arc is checked from the source matrix. The unique permutation carrying P
to Q by their clock positions is not a source automorphism. This is a
genuine three-clock local reconfiguration beyond kernel-lift transport,
while keeping the target owner as a port mark.

Each of the other three sources has a saved analogous width-three
non-automorphism move with its target retained. If only the original phase
and endpoint owners are marked, there are also local changes of target:
in source 2246,

    (0,7,6,5,2,4,3,1) -> (0,7,6,5,4,3,2,1)

keeps phase 6 and the outside clocks but changes target 2654 to 1296.
That move is excluded by the target-fixed menu. No per-port target
preservation is inferred merely from a matching phase.

## 4. Same phase is not complete wiring conservation

The incoming rectangles preserve selected source/target counts and the
original mod-eight Kirchhoff equations, but their higher dyadic debts
remain. Their apparent doubled aggregate capacities fail because an
actual target has only one distinct path. This is a stronger constraint
than the one-path reconfiguration problem here.

Our uniform random intermediate in PROOF.md routes one specified return
demand; it may reuse an observer on different trials. It is not an injective
assignment of many simultaneous selected ports. Applying it to a complete
wiring would require capacities for distinct observers and compatible
source/target incidence, as well as the electrical ledger. This connection
identifies a new obligation rather than silently replacing it by flow mass.

## 5. Reproduction and exact input handling

`electrical.py` rebuilds each source from the committed canonical bit word,
enumerates every permutation, checks all Hamiltonian junctions, obtains
both owners using the ORIGINAL forward-bit atlas, and evaluates the frozen
potential. It checks literal reversed-converse paths too. Transition groups
retain the full outside-owner tuple; BFS routes and exact rational edge
capacities give the congestion receipts.

The first provenance check exposed a text-serialization difference: the
inherited potential hash was exactly the CRLF view of the committed LF
file. Both hashes are recorded and explicitly checked, with no mathematical
change. Concurrent Session44 independently repaired the original metadata
in [its receipt review](../session44_guarded_native_moves_20261023/INCOMING_RECEIPT_REPAIR.md).
Our binary owner-map hash remains an exact raw-byte check. The checker
accepts only the two documented text views and pins the LF content itself;
it does not ignore arbitrary hash changes.

The checker consumes inherited exact atlas data; it does not re-prove the
entire order-eight class enumeration or harmonic-potential construction.
