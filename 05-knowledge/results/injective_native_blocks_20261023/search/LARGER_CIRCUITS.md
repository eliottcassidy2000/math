# Sixteen sources can pay together when every proper supplied subblock fails

23 October 2026. **FINITE-EXACT:** injective8-,10-,12-,16-source circuits in
the original **order-eight** tournament metagraph. **PROVED:** menu-relative
minimality by complete finite state checks and the exact current-code
criterion. Sixteen refers to selected source slots, not tournament order.

## 1. The strongest retained example

All sixteen sources have the same two eligible target owners,
U=1965 and V=2234. Each actual source offers its two saved Hamiltonian paths,
with actual endpoints, an indexed shared arc, every unchanged clock, its
chosen original mod8 phase, and its literal converse partner retained.

Eight sources move U to V:

    176,1477,1691,1974,2071,2087,2185,2905.

Eight sources move V to U:

    1394,1425,1457,1638,1841,1921,2553,2581.

Their original integer potential values satisfy

    197+152+114+81+139+107+138+88
      =113+158+77+94+99+217+104+154
      =1016.

Since X(U)=X(V)=146, the entire selected current difference annihilates X
**over the integers**. The exact path table is
`larger16_selected.json`; `larger16_receipt.json` independently authenticates
all6734 Hamiltonian paths of these sixteen canonical source tournaments.
Every source automorphism group is trivial. Every candidate source/target
edge is pure black, and its complete actual port capacity is saved.

The source and target inventories are both injective. In either endpoint
state there are32 distinct normalized ordinary/converse source words and32
distinct inverse target words, with these sets disjoint. Their union over
the two states contains128 distinct endpoint words, saved as the exact
reservoir to exclude from an unchanged background. Thus no path or inverse
target port is recycled. The inverse ports are named by literal free-word
complementation and actual K_P matrices, not an arbitrary target isomorphism.

## 2. Why the whole block passes all original dyadic currents

For each of these nonzero unit-entry balanced source vectors w and
v=eU-eV the paired-black difference is

    DeltaL = 2(wv^t+vw^t).

The disjoint source/target supports make preservation of every original
modulus-m harmonic character equivalent to

    2w in im(L mod m),  2v in im(L mod m).

The complete effective paired source code is authenticated in
`../electrical/ELECTRICAL_CODES.md`:

    Sigma(j)=(X(j) mod128,Y1(j) mod2,...,Y4(j) mod2)
      in Z/128 x (Z/2)^4.

For zero-sum integer charges d, `Sigma(d)=0` is equivalent to preservation
of every original dyadic harmonic character by the paired charge2d. This
uses the actual native2-primary critical group, not a guessed collection of
vertex signatures. Both v and the displayed16-source w have zero code.
The search checker tests the source relation exactly; the code theorem
supplies complete-bank sufficiency. For the12- and16-source cases there is
also an independent, more direct upper certificate:
`../electrical/larger_image.json` saves full vectors z with
`Lz=2w mod256`, including their complete integer correction ledgers. The
checker multiplies them by all438216 actual original edge rows, and also
replays the separate image witness for2v. Thus these two positive moves do
not depend on inferred completeness of a proposed bank. Extending from256
to every dyadic modulus still uses the authenticated exponent256 theorem.

This does **not** preserve every character's phase at each individual slot.
The authentic original mod8 harmonic4e1965 separates U and V. Aggregate
current cancellation and per-slot universal phase preservation are different
requirements. This construction establishes the former, while retaining
the specifically selected original phase at every slot.

## 3. Exact minimum under the supplied binary menus

Choose before or after independently for each of the sixteen slots. All
65536 possibilities retain their local native legality and pins; sources
are distinct, so the individual path resources do not collide. Impose the
exact target degree counts and the original universal dyadic current bank.
Exactly two possibilities survive:

    0000000000000000, 1111111111111111.

Thus every proper nonempty subblock fails a degree or current obligation.
The sharp simultaneous slot budget in this supplied menu is16. The total
changed owner-clock occurrence count is36. This is a finite statement about
these actual choices: it does not exclude detours through extra paths,
different pins, vacant slots, or other source classes.

The obstruction is genuinely sensitive to higher dyadic information.
The exact accepted state counts, always retaining target margins, are:

|Universal original modulus|4|8|16|32|64|128|256|
|---|---:|---:|---:|---:|---:|---:|---:|
|Accepted states|378|194|86|38|18|4|2|

For2<=k<=8 the complete paired code at modulus2^k is
`(X mod2^(k-1),Y1 mod2,...,Y4 mod2)`. Above256 the effective dyadic code
saturates at the original critical group's exponent256; one must not invent
extra bits by reducing this particular saved X modulo ever larger powers.
The layered table makes precise how low-modulus repairs disappear when
their missing information is restored.

## 4. Smaller independently retained controls

`larger_selected.json` also contains primitive8-,10-,12-source circuits,
each with disjoint actual sources and the same target pair. Their two
integer source sums are respectively484,693,786. Complete checks of their
256,1024,4096 binary states leave only the two endpoint states, giving sharp
menu-relative slot budgets8,10,12. All actual path, converse, inverse-port,
pin, capacity and current checks are independent of the16-source case.

These do not supersede the earlier shared-owner eight-slot packet in
`README.md`. That packet has four states and factors into a2-slot exchange
and a6-slot circuit. Keeping the reducible example is a useful hostile
control against mistaking a large written block for an irreducible move.

## 5. Bounded search, recurrence, and scope

The source menu was the28 actual eligible sources for U,V found by the
earlier bounded order-eight search. The larger search samples equal-sized
source subsets with fixed seed20261023 and buckets them by their **exact**
X sum and four binary code sums. On a disjoint collision it tests every
proper equal-cardinality signed subset against the full effective code.
This is a finite meet-in-the-middle construction; it makes no random-model
claim about the actual metagraph.

The8-,10-,12-source cases required139,425,210 sampled subsets. The16-source
case required1014 sampled subsets and five disjoint balanced candidates.
The search stopped at16 as a deliberate session bound; neither16 nor any
other value here is asserted to be a maximum. No new global wiring rule is
claimed.

The reusable mechanism is **combine relations, retain actual ports, then
test proper subrelations**. Ordinary graph degree-preserving moves can be
split into rectangles, but those rectangles need not preserve the electrical
character code. The larger zero-sum source relation is the extra information
that lets several source pairs pay their current debts together.

This changes selected native inventories with an explicit unused reservoir.
It does not replace the fully occupied all-path atlas by a new realizable
atlas, assert equality of old/new harmonic banks, or preserve odd-modulus
currents.

Reproduce the bounded deterministic discoveries and small exact checks:

    python 05-knowledge/results/injective_native_blocks_20261023/search/larger_search.py
    python 05-knowledge/results/injective_native_blocks_20261023/search/larger_check.py
    python -O 05-knowledge/results/injective_native_blocks_20261023/search/larger_check.py
    python 05-knowledge/results/injective_native_blocks_20261023/search/larger_search.py 16
    python 05-knowledge/results/injective_native_blocks_20261023/search/larger_check.py --sixteen
    python -O 05-knowledge/results/injective_native_blocks_20261023/search/larger_check.py --sixteen

The checker does not depend on discovery caches. Its all-bank inference
depends only on the pinned original native atlas and the separately audited
electrical code theorem; individual character checks alone are not presented
as a complete bank.
