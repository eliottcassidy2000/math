# Independent erasure-section audit

Status: **VERIFIED** proof mechanism and exact finite computation,2026-10-24.
Reviewed `check.py`, `receipt.json`, `audit.py` and the completed
`ERASURE.md`. No producer changes were needed. The finite table and the
subsequent universal finite routing consequence are distinguished below.

## Exact recurrence and independent representation

Removing an actual source slot of charge a from a complete subset-count
table satisfies E_k(q)=D_k(q)-E_(k-1)(q-a). The producer uses this recurrence
with the already-computed lower row and omits the old top degree, exactly as
required. Sequential removal applies it once per distinct erased source;
no record or source multiplicity is silently reused. Nonnegative counts,
total binary assignment count and complement-symmetric coverage profiles
are checked after every deletion.

The independent checker imports no producer counting routine and takes its
charge coordinates from the original electrical potentials, not shortened
producer labels. It recomputes each subset image in the Boolean group
algebra. Its parity-bit permutation followed by cyclic shift implements
translation in Z/128 x F2^4;4096 singleton controls check the bit order and
wrap. All1351 pin-set profiles agree with the integer-count recurrence.

Normal and optimized Python reproductions both pass:

* `check.py`:5,207 checks;
* `audit.py`:5,242 checks, all1351 deletion sets and all1140 repairs.

## What the unary pins mean in the actual host

The20 cheap core sources have two fixed full records differing at exactly
two owner clocks. At either changed clock, fixing the literal U owner or V
owner selects exactly that record from the complete inherited binary menu.
Such a new pin does not overwrite any inherited common-clock pin. The core
cells have private owner sets in the retained XOR host, so distinct such
unary pins do not create an untested overlap equation. An additional direct
audit reconstructed these40 possible unary pins using the actual global
owner frames and verified singleton acceptance in every case.

Thus freezing e of these slots is faithfully represented by deleting them
from the *variable* bank. Their physical slots remain occupied. If p of the
frozen records select V, subtract p from the required variable occupancy and
subtract their actual Q-charge from the required charge. No source, edge,
converse partner or fixed contribution has disappeared. More general
relational pins cannot be replaced by whole-slot deletion without proof.

## Quantifiers and repair boundary

For all pin sets of sizes0,1,2, the universally complete remaining grade
bands are respectively5..15,6..13,7..11. At size3 there is no one complete
grade shared by every pin set. Four triple deletions have no complete grade;
two other triples have complete grades7 and10 but incomplete grades8 and9.
These are distinct phenomena, so neither symmetry nor neighboring complete
grades implies central completeness.

For every one of the1140 triple deletions, adding the already existing
unfrozen core source1841 to the receiver bank yields complete grades8,9,10
on18 free slots. Its switch costs3 actual owner-clock occurrences, it is
outside the20 cheap sources, and no new native resource is created.
The count proof is the union of the17-slot grade-k image and the shifted
grade-(k-1) image; the Boolean audit independently rebuilds the18-slot bank.

The section table alone applies at receiver grades8..10. The completed
`ERASURE.md` additionally proves the global routing result in the next
section, using receiver migration in the interior and direct feasible
endpoint replacement outside it. Neither statement covers erasing1841 as
an additional pin, arbitrary source deletions, or four or more frozen cheap
slots. The fixed external helper owner4 pin remains unchanged throughout.

## Global routing with the three extra unary pins

After accounting for the three frozen core records, there are25 free core
slots: the18-slot receiver and seven outside sources. Its total possible
clock work is37=17*2+3. For free core V population t in8..17, every required
canonical receiver grade max(8,t-7) lies in8..10. The receiver migration
conditions are therefore satisfied with the genuine native records.

Entry below grade8 or above grade10 is a simultaneous cardinality-balanced
replacement. The usual subset bound, using at most2d+1 clock work for d
changed receiver slots and4 per changed outside source, gives the stated
49 upper bound (some endpoint constraints sharpen it, but no sharpening is
needed). Consecutive receiver grades differ at at most17 slots, costing at
most35; one outside change gives39. An outside Johnson exchange costs at
most8 plus37 for the receiver, giving45. A two-helper exchange costs at
most15 plus37, giving52. All accepted events change at most20 slots. These
are uniform bounds for arbitrary choices of valid receiver sections, so
the prescribed endpoint normalizations can be reversed exactly.

For t<8 or t>17, the selected or unselected free-core endpoint subsets have
size at most seven. Their symmetric difference has at most14 slots, and
the two helper pairs differ at at most four. Direct replacement of the
two already feasible endpoints is therefore an accepted event with at
most18 slots and85 clocks. No auxiliary grade or unguarded complement event
is introduced. Consequently the completed proof gives20-slot/52-clock
connectivity in the interior, and20-slot/85-clock connectivity for EVERY
nonempty fixed-margin/full-Q fibre under the three extra unary pins.

The frozen records' selected count and charge are constants in this proof;
they are not silently discarded. All-Q coverage is guaranteed only in the
interior, whereas boundary connectivity is conditional on feasible
endpoints. This is a valid additional routing corollary of the actual
18-slot repair, beyond the finite grade table itself.
