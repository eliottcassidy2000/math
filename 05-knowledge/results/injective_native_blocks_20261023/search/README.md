# Native current circuits with injective paths and shared owner coordinates

23 October 2026. **FINITE-EXACT:** two four-source circuits, their six-source
elimination, and a complete eight-slot binary-menu calculation. **PROVED:**
the current-transfer formulas and the stated menu-relative minimality from
the exact certificates. These are selected native inventories, not a new
realization of the complete all-path tournament metagraph.

## 1. The inherited obstruction tells us what to change

Session44's `rectangle-bank/THEOREM.md` excludes the old two rectangles:
the full mod8 bank requires even amplitude, but the actual scarce target
ports allow only amplitude one. Merely combining more copies recycles a
path. The new `../electrical/` lane additionally proves that the four old
sources have no nonzero balanced unit-entry character relation. This is a
source-side obstruction to every such trade on those four sources, not only
to the two previously selected target pairs.

We retain the source charge, actual path coordinates, and target charge,
and change the source menu. The resulting mechanism is a relation among
several **different** sources: individual rectangles fail, while a larger
block has zero total electrical response.

## 2. A literal four-source block

Let U=1965 and V=2234 be original order-eight converse-merged target owners.
The ordinary source representatives, literal paths, forward masks, endpoint
pins, an indexed arc, and **every** unchanged owner clock are saved in
`selected_blocks.json`. Owner labels are local to each canonical source.

|Source|Before path|After path|Target change|Changed clocks|
|---:|---|---|---|---|
|1477|7 6 3 1 4 5 2 0|7 6 5 1 4 3 2 0|U to V|2,5|
|2837|1 0 6 2 7 5 4 3|1 4 6 2 7 5 0 3|U to V|1,6|
|1425|5 2 3 4 7 6 1 0|5 6 3 4 7 2 1 0|V to U|1,5|
|1921|7 6 4 1 5 3 2 0|7 2 4 1 5 3 6 0|V to U|1,6|

Clocks are zero-based. All four source automorphism groups are trivial.
Each transition is a genuine two-owner Hamiltonian-path exchange. The
unchanged endpoint and indexed-arc coordinates are authenticated in the
original source adjacency matrices. The previous local strong-four-vertex
exchange criterion supplies a useful discovery primitive; it does **not**
by itself certify any electrical phase.

All eight involved source-target edges are pure black. Their black
conductances, ordered U,V, are

    1477: (2,2); 2837: (4,2); 1425: (4,2); 1921: (2,2).

Each selected ordinary path and its literal reversed path in the converse
tournament contribute conductance2. The eight ordinary paths across the
two endpoint states are distinct by source/path; each endpoint inventory
has four distinct ordinary paths and four distinct converse partners.
The checker enumerates all source paths, rather than inferring capacity
from a black weight.

Put

    wA = e1477 + e2837 - e1425 - e1921,
    v  = e1965 - e2234.

The selected before-to-after current matrix difference is

    DeltaL = 2(wA v^t + v wA^t).                         (1)

The factor2 is the actual converse-paired black multiplicity. Sources and
targets are disjoint; all row/column degree margins are preserved.

For the original saved integer representative X of the mod256 harmonic
potential,

    X(U)=X(V)=146,
    X(1477)+X(2837)=152+223=158+217=X(1425)+X(1921).

Consequently **DeltaL X=0 over the integers**, not only modulo256. This
does not make X an integer harmonic function of the entire graph: the
original statement is LX=0 mod256. It says the selected block changes none
of this particular integer vector's currents.

The electrical lane supplies explicit full image witnesses for 2wA and 2v
modulo256. Direct multiplication by all438216 original edge rows verifies
them. Image/kernel duality therefore proves that (1) annihilates **every
original mod256 harmonic character**. Together with the independently
authenticated exponent256 of the original 2-primary critical group, the
same holds at every dyadic modulus. It asserts inclusion of the old kernel
in the kernel of an algebraically perturbed Laplacian, not equality of
kernels, and not preservation at odd moduli.

## 3. Circuit elimination produces a six-source block

A second actual four-source relation uses the same target pair:

    wB = e1477 + e1841 - e1394 - e2185,
    X(1477)+X(1841)=152+99=113+138=X(1394)+X(2185).

Its complete paths and pins are in the same selected packet. Reverse this
block and combine it with A. Their source1477 contributions cancel:

    wC = wA-wB
       = e2837 + e1394 + e2185 - e1425 - e1921 - e1841.

This is a **six-source** actual current circuit, with

    223+113+138 = 158+217+99 = 474.

The image witness for 2wC is the difference of the witnesses for 2wA and
2wB. No new harmonic-bank assumption is needed. The six native before
paths and their six after alternatives are distinct, retain all declared
endpoint/arc/unchanged-clock pins, and preserve the chosen per-slot phase.
All relevant actual target capacities and literal path lists are recorded
in `receipt.json`; there is no repeated source in this six-slot inventory.

There are exactly two admissible states among its2^6 supplied binary
choices, once target margins, original X, and the authenticated paired
character tests are imposed: all before, and all after. Hence its **sharp
simultaneous slot budget is six within this declared menu**. Its owner-clock
occurrence cost is13. This is not a global lower bound allowing other paths,
extra vacant slots, changed pins, or other source tournaments.

## 4. The shared-owner eight-slot packet and its honest simplification

Keep both occurrences of source1477 in the two four-slot blocks, making
eight slots over seven sources. At this source block A uses P to Q while
the reversed block B uses Q to P. Both endpoint inventories contain P and
Q exactly once. We impose two actual shared owner equalities:

    slot0[2] = slot4[5],
    slot0[5] = slot4[2].

They change from values(3,5) to(5,3). Thus these are dynamic shared owner
coordinates, not merely two slots carrying the same constant label. Moving
either original four-slot block alone violates these equalities and also
repeats one path at source1477.

The complete2^8 check nevertheless finds **four**, not two, admissible states:

    00000000, 01110111, 10001000, 11111111.

The two source1477 slots can exchange separately, and the other six can
perform circuit C. The accepted state graph therefore has sharp connected
move budget6, with independent two-slot and six-slot moves. It would be
incorrect to claim that the eight-slot packet is irreducible. This is a
positive reduction: retain the shared coordinates while eliminating an
internal source to obtain a larger genuine six-source current circuit.

## 5. Preserved data, lost data, and inventory scope

**Source:** actual canonical tournaments, paths, fixed endpoints, literal
indexed arcs, unchanged owner clocks, and the two explicit crossed shared
coordinates. **Target:** an injective selected converse-closed port inventory
and its source/target incidence matrix. **Map:** the displayed path exchanges,
with actual K_P target owners authenticated by the original owner map.

**Preserved:** every declared shared coordinate; source and target counts;
pure-black pairing; chosen per-slot mod8 phase; all original dyadic aggregate
electrical currents; and the original X current vector over the integers.

**Not preserved:** universal electrical phase per individual slot. The genuine
original harmonic character4e1965 modulo8 changes by4 between U and V, so it
distinguishes every target-changing slot. Aggregate cancellation must not be
reported as per-slot universal preservation. Odd-modulus currents, kernels
of newly constructed graphs, and the full all-path inventory are not proved
unchanged.

The checker also authenticates injectivity at the **target** end. In each
state the16 ordinary/converse normalized forward words are distinct, their
16 complementary inverse target words are distinct, and these two sets are
disjoint. The literal K_P involution restores the original matrix on the
same owners. No arbitrary isomorphism of a target class is needed to name
or distinguish its inverse ports. The union over before and after contains
56 actual endpoint words; every one is saved in the receipt.

Let R be any unchanged selected background inventory disjoint from **all
56 source and inverse target endpoint ports** in that reservoir. Then R plus the before choice
and R plus the after choice are both injective; no path is cloned. The
complete native all-path inventory already contains both sets, so replacing
a subset inside that fully occupied inventory would require additional
global reassignments. We do not claim such a full native rewiring.

## 6. Bounded discovery and reproduction

Discovery used all3352 non-self-converse merged sources at order8, taking
the canonical `class_a` representative, and all1040912 of their Hamiltonian
paths. It grouped paths by endpoints, chosen original phase, and one fixed
internal arc at clock1 through5. It retained target changes with equal paired
mod8 character code. There were11138 source/target-pair menu entries across
1398 target pairs. Pair-sum lookup in the full retained source code found252
four-source candidates with a target X difference divisible128 and no legal
balanced two-slot subblock. These are discovery-universe counts, not a
classification of all native moves or all pin schemas.

The large discovery caches are intentionally untracked. To reproduce them:

    python 05-knowledge/results/injective_native_blocks_20261023/search/explore.py all
    python 05-knowledge/results/injective_native_blocks_20261023/search/find_blocks.py
    python 05-knowledge/results/injective_native_blocks_20261023/search/select.py

The decisive small checker needs none of those caches. It independently
enumerates2739 paths of the seven selected source tournaments, authenticates
converse paths and actual K_P matrices, all target capacities, all binary
states and the stated exact response identities:

    python 05-knowledge/results/injective_native_blocks_20261023/search/check.py
    python -O 05-knowledge/results/injective_native_blocks_20261023/search/check.py

Run the separate electrical checker for complete image certificates and
the all-dyadic exponent/completeness statement. Text hashes normalize CRLF
to LF; compressed atlas hashes retain their original bytes.
