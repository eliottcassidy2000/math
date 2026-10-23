# A genuinely coupled four-source native circuit with no recycled ports

23 October 2026. **FINITE-EXACT:** an injectively realized selected native
packet changing eight merged incidences, retaining every declared physical
owner/arc pin, and preserving the full original harmonic banks modulo eight
and modulo 256. Each of its two constituent rectangles fails both the full
mod-eight bank and the chosen modulus-256 potential separately; their debts
cancel. **PROVED:** the conditional all-size source-circuit compilation rule
in section 6. No complete canonical-atlas rewiring is claimed.

## 1. Inheritance and the changed carrier

The original `native_guard_moves_farey_20261023/native_moves` packet kept the
actual owners and converse factor two but its fixed source pair
1238/2246 had scarce one-port targets. Session44's `rectangle-bank` proved
that its two rectangle charge orders are four modulo eight: the next
amplitude is electrically required but has no second distinct actual port.
Session44's `coupled-moves` retained shared variables and exact matching;
its strongest fully guarded examples did not change merged incidences.
The Session44 `electrical-duality` row-image criterion supplies the correct
universal-bank test. The new operation changes the source pair, not the
meaning of a port or of the guard.

Concept board: full actual path / explicit common physical-owner frame /
K-converse resource tile / full-bank image charge / higher potential debt /
primitive coupled move. The key representation change is from one source
difference to a signed four-source charge. Neither the graph's total black
capacity nor a repeated local owner number is used as a resource certificate.

## 2. The explicit four-source packet

Let the targets be U=1965 and V=2234 and put

    q = e1425 + e1921 - e1477 - e2837,
    v = e1965 - e2234.

The positive slots move U to V; the negative slots move V to U. The original
canonical-owner paths, correctly oriented, are:

| Source | Before | After | Target change |
|---:|---|---|---|
|1425|5,6,3,4,7,2,1,0|5,2,3,4,7,6,1,0|1965 to 2234|
|1921|7,2,4,1,5,3,6,0|7,6,4,1,5,3,2,0|1965 to 2234|
|1477|7,6,5,1,4,3,2,0|7,6,3,1,4,5,2,0|2234 to 1965|
|2837|1,4,6,2,7,5,0,3|1,0,6,2,7,5,4,3|2234 to 1965|

Every pair changes exactly two owner clocks. Every unchanged owner clock
is literally pinned, so each slot has six fixed indexed owners, including
both endpoints and an internal actual arc. The complete corresponding
six-pin Hamiltonian family is exactly the displayed pair: only two owners
remain for two clocks, and both permutations are actually Hamiltonian.

The checker recovers the complete source path inventories independently
by DFS and by all 40,320 owner permutations. Their sizes are respectively
313,323,259,559. Each source has a trivial automorphism group, independently
tested on all permutations, and ordinary and converse source classes differ.
These are not invisible source-automorphism changes.

### Actual common physical frame, not a numerical coincidence

Use physical owners 0 through 7. The explicit bijections from canonical
local owners 0,...,7 to these physical owners are:

| Source | Local-to-physical owner map |
|---:|---|
|1425|7,1,2,3,4,0,5,6|
|1921|7,4,1,2,3,5,6,0|
|1477|7,3,1,2,4,5,6,0|
|2837|1,0,4,7,2,5,3,6|

The source matrices are transported by these same bijections, not left in
their old gauge. `packet.json` retains all four original and physical 8 by 8
matrices, both complete owner paths per slot, and every indexed pin.
Every slot has physical endpoint owners 0 and 7. Every slot also pins the
same physical actual arc 3 to 4: at clock two in slots 1425,1921,2837 and
clock three in slot1477. Thus the shared owner/arc coordinates have explicit
occurrence maps. The other two fixed owner clocks per slot are retained too.
No unspecified equality between unpinned clocks is invented or discarded.

For example, the physical after-path in source1425 equals the physical
before-path in source2837. This does NOT reuse a port: the fixed physical
source matrices differ, and their normalized forward words and resource
tiles differ. A resource is not an untyped permutation of owner labels.

### Full literal and resource injection

For an actual path P in a physical tournament T, retain all four literal
records

    (T,P), (K_P(T),P), (T^transpose,reverse(P)),
    (K_P(T)^transpose,reverse(P)).

The checker verifies every matrix arc, the K involution, and the full
converse square. In the original forward-word gauge, F complements the 21
nonconsecutive bits and R reflects their clock indices. The four words
{w,Fw,Rw,FRw} are a pure-black resource tile of size four. Resource identity
is the whole orbit, represented by its minimum word, not a path count.

The before/after tile identifiers for the four slots are:

| Source | Before resource | After resource |
|---:|---:|---:|
|1425|284821|317909|
|1921|660606|923390|
|1477|423741|645314|
|2837|284377|21593|

All eight tiles are distinct. Each endpoint therefore has four selected
ordinary-source paths, eight source-facing ordinary/converse presentations,
and sixteen ports in its fully K-and-converse-closed packet. Across both
endpoints there are 32 distinct literal matrix-and-path ports. An independent
resource audit separately reconstructs these facts from the original atlas.

The original black weights before/after are (4,2),(2,2),(2,2),(2,4), with
zero blue weights throughout. The complete ordinary target capacities are
(2,1),(1,1),(1,1),(1,2). Thus every selected use is real, and none is justified
by pretending that aggregate weight creates an additional port. Including
the actual converse partner contributes black coefficient TWO per slot.

The nonzero signed edge change is therefore

    -2: (1425,1965), (1921,1965), (1477,2234), (2837,2234);
    +2: (1425,2234), (1921,2234), (1477,1965), (2837,1965).

Source degrees and target multiplicities are preserved. This changes eight
actual selected incidences, rather than permuting an unchanged edge multiset.

## 3. Universal banks and the exact higher-layer cancellation

Let L be the original symmetric integral Laplacian. The signed change is

    Delta L = 2(q v^t + v q^t).                                 (1)

`image_witnesses.json` supplies two full grounded vectors z,w with entries
in 0,...,255, plus complete integral correction ledgers a,b, satisfying

    Lz = 2q + 256a,       Lw = 2v + 256b.                        (2)

The checker multiplies these vectors by every original exported edge. It
does not trust the witness generator, its elimination pivots, or a purported
complete harmonic bank. Symmetry of L gives, for every x with Lx=0 mod256,

    2 q.x = 2 v.x = 0 mod256.

Equation (1) then gives Delta L x=0 mod256. Reduction of (2) modulo eight
proves the same inclusion for EVERY original harmonic character modulo
eight, not just those reducing from a chosen saved potential. This is
preservation of the old harmonic bank, not equality of old and new kernels.
A full-bank all-dyadic extension may be combined with the independently
authenticated exponent-256 theorem; it is not inferred from one character.

For the original chosen primitive modulus-256 lift X, the six values are

    X1425=158, X1921=217, X1477=152, X2837=223;
    X1965=X2234=146.

Hence q.X=158+217-152-223=0 and v.X=0 over the INTEGERS. The complete frozen
current debt is exactly zero, a stronger statement than vanishing mod256.
Every individual source slot also keeps its original phase modulo eight.

The two constituent rectangles are (1425,1477;U,V) and
(1921,2837;U,V). Their frozen debts, in target coordinates (U,V), are

    (12,-12),        (-12,12).                                 (3)

Every other coordinate is zero. Each rectangle separately fails the chosen
mod256 guard and already has a nonzero frozen mod-eight aggregate debt.
Their sum cancels it exactly. This is not a bundle of independently safe
rectangles, nor an amplitude-two reuse of one scarce source pair.

There is also an independent universal-bank separator. The original actual
harmonic character 2 times `radical_lifts_mod4[2]`, authenticated against all
original edges, gives each rectangle target defect (4,4) modulo eight.
These nonzero defects cancel when the rectangles are coupled. Lower-order
characters and the selected higher potential are checked separately.

**The per-slot guard distinction is essential.** This packet preserves the
original per-slot X phase and every original aggregate Kirchhoff harmonic
bank, NOT every per-slot phase for every harmonic character. The actual
character y=4e1965 modulo eight is harmonic: every incident nonloop original
weight at1965 is even. Its target difference y1965-y2234 is four, so each
individual U/V port phase changes by four. The paired coefficient two kills
that difference at the aggregate level. The checker authenticates this
separate hostile explicitly; no nonunit is cancelled.

## 4. Exact minimum simultaneous budget in the complete pinned menu

There are exactly sixteen Cartesian choices from the four complete two-path
six-pin families. Target multiplicities must remain two U and two V.
Every nonempty proper balanced change replaces one positive and one negative
slot, hence is one of four two-source rectangles. The actual separator 2Y2
rejects every one of them. Odd numbers of changed slots fail the target
margin. Consequently the universally mod-eight-safe balanced fibre consists
of exactly the all-before and all-after packets.

Thus four slots must change simultaneously within THIS complete pinned
menu; there is no legal intermediate pinned packet. Each changed owner
permutation costs at least two clocks, and all four attain that lower bound.
The total cost eight counts the four ordinary-source slots and is minimum
for a move changing these four permutation slots. Counting their actual
converse copies would give sixteen changed clock occurrences; retaining the
full K/converse closure gives thirty-two. These are different cost conventions. This is not a global minimum over all native menus or all orders.

This also identifies why the fixed-pair tariff is evaded. The target rows
see the combined source charge q, not four independently tested pair charges.
The scarce per-slot capacities remain one, but the missing parity and higher
current are supplied by distinct physical ports at new source owners.

## 5. The expanded finite source menu

`search_native_menu.py` independently scans all original nonselfconverse
canonical source owners. It retains the declared trivial-automorphism
carrier (equivalently, twice the complete path count equals the original
merged word fibre), excludes the two target owners, and checks every actual
U/V path pair for common endpoints and an indexed internal arc. The exact
source census and a best full path pair per eligible source are retained in
`target_pair_inventory.json`; it contains 28 eligible sources and 28 distinct
mixed electrical labels. Each record retains complete target capacities.
This is NOT a census of selfconverse or nontrivial-automorphism carriers.

The 28-source menu supports substantially more source relations than the
old fixed pair. Its relation count is not inferred from black capacities.
The four-source packet above is explicitly constructed and checked rather
than justified only by a pigeonhole estimate. An abstract source-relation
bound requires actual two-target options at each source; a count of all
3528 network vertices is not such a native option inventory.

## 6. Conditional all-size source-circuit compilation

**Theorem.** Fix any finite native network with symmetric Laplacian L,
modulus m, distinct target vertices U,V, and distinct source vertices S_i
disjoint from those targets. For each slot i supply two actual observers
P_i^U,P_i^V of source S_i with those actual target owners. Retain their
matrices, all shared-owner occurrence maps, literal pins, and K/converse
closure. Assume both endpoints described below are compatible, every
selected resource tile is distinct, and the actual paired convention
contributes coefficient c. These are hypotheses about real objects, not
consequences of an electrical equation.

Choose signs sigma_i in {+1,-1} with sum sigma_i=0. Start the positive slots
at U and negative slots at V, and switch them all. Set

    q=sum sigma_i e_(S_i),       v=e_U-e_V.

Then this is an injective balanced selected native switch. It preserves
EVERY original harmonic character modulo m if and only if

    cq and cv belong to im(L mod m).                            (4)

If a fixed potential X also has q.X=v.X=0 over the integers, its entire
current is unchanged over the integers. If every new selected tile avoids
an unchanged global inventory, the switch extends injectively over that
inventory. The full canonical all-path inventory need not satisfy that
avoidance hypothesis.

**Proof.** The actual supplied paths and compatibility hypotheses authenticate
the two endpoints; tile distinctness authenticates injection and the stated
avoidance gives the extension. The sign sum gives the target margin and
each source still contributes one slot. Adding actual edge currents gives
Delta L=c(qv^t+vq^t). At any source row the defect is +/-c(v.x); at target U
it is c(q.x), and at V its negative. Thus vanishing on the whole harmonic
kernel is equivalent to cq and cv annihilating that kernel. For a symmetric
matrix over Z/m, finite perfect character duality identifies the annihilator
of its kernel with its image. This proves both directions of (4). The exact
integer-potential assertion follows from the same formula. QED.

The theorem applies in every finite order but does not assert that every
source has suitable observers, that independent shared-frame kernels glue,
or that a uniformly small circuit exists. Repeated source coefficients need
additional distinct actual ports and full shared-lift compatibility; they
cannot be implemented by repeatedly naming the same slot. All-size native
menu existence, bounded complexity and complete canonical rewiring remain
open.

## 7. Reproduction, provenance and deliberate non-consequences

`SOURCE_PINS.json` fixes the original atlas, owner map, edge export and
electrical receipt by byte hashes. The certificate checker uses explicit
exceptions, never optimization-sensitive assertions or floating arithmetic.

    python3 -B check.py --repo /path/to/math --out receipt.json
    python3 -O -B check.py --repo /path/to/math --out receipt.optimized.json
    cmp receipt.json receipt.optimized.json

Optional full image-witness regeneration and expanded-menu census:

    python3 -B solve_images.py --repo /path/to/math --modulus 256 --out image_witnesses.json
    python3 -B search_native_menu.py --repo /path/to/math --out target_pair_inventory.json

The image generator uses exact bit-plane elimination and then independently
checks its vectors. Its ranks are not premises of the certificate checker.
The path checker uses two independent enumeration methods, actual color
counts, the complete 16-state hostile-control table, explicit shared-owner
maps, literal K/converse matrices and the complete original edge action.
No claim here settles LRC(14), Collatz, unrestricted all-size wiring, or
realizability of a modified full canonical all-path atlas. The guarded
selected-packet construction is the precise positive result.
