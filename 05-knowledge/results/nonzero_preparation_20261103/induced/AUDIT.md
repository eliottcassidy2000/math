# Independent audit of the complete induced-extension obstruction

**PASS — 2026-11-03.** The coordinator audited the mathematical reduction,
the complete search source, its exact arithmetic, and an independent replay.
This accepts the statements of [PROOF.md](PROOF.md), including weak final
extensions. It does not accept an obstruction to old-arc surgery or to all
tournaments with the same Hamilton count.

## Analytic audit

The three terms in the arbitrary-birth law distinguish a first, last, and
internal new vertex. The internal term retains both ordered endpoints;
diagonal endpoint entries vanish. Expanding cut bits gives the stated
quadratic recurrence with the correct minus sign. Fixed insertion before
the first beaten old path vertex (or at the end) is injective after deleting
the newcomer. Consequently every prefix of a requested endpoint survives
the count bound, even when that prefix is weak. There is no class quotient,
port restriction, or early strongness filter.

The primary Busch theorem was read directly: the strictly increasing strong
minimum is 5625 at order17. Thus the search through five new vertices is
complete for strong cores of counts4289 and4291. The weak reduction is also
valid: the strong old source lies wholly in one SCC, whose count is at least
1429. An outside odd factor greater than1 would be at least5, because the
requested totals are not divisible by3. It would make the containing core's
count smaller than1429. The outside factor is therefore1, and every outside
SCC is a singleton. This pays all larger weak orders without enumeration.

The labelled weak +4 endpoint formula is independently recovered. A
one-birth core has d choices for its label, (d−1)! orders of outside labels,
and d splits between prefix and suffix, giving d!*d. A two-birth core has
binomial(d,2) label choices, six placements in the old y clone block,
(d−2)! outside orders and d−1 splits, giving 3*d!*(d−1). Their sum is
d!*(4d−3). The root response scales by 2 per prefix vertex and is unchanged
by suffix vertices. Neither nontripling one-birth F4677 nor a positive
prefix can meet F4689; clone cores with empty prefix do meet it.

The inherited clone response was checked against its ordered path slots.
The equation for another enlargement has only the exponent pair1,3, so one
physical +4 does not close the branch. The one-birth all-clock formulas have
the correct differences 9z−18 and 9z−30; their clock-one coincidence is
explicitly outside an all-clock conclusion.

## Source and independent computation

The audited source is
[`induced_complete.cpp`](../../../../04-computation/nonzero_preparation_20261103/induced_complete.cpp),
SHA256 `3e2647525cc79a33788031e5b2e11dad425eb6744fa7cefb93215bc0ba3d83f1`.
An independent coordinator build with `c++ -std=c++17 -O2` completed all
**150,939,321 always-active gates**. Its output is byte-identical to the
producer's -O0/-O3 frozen receipt
[`induced_complete.json`](../../../../04-computation/nonzero_preparation_20261103/induced_complete.json),
SHA256 `9247ce871dde82aaca8b441e23fb2a07bd69c0a018f13e619f7a33af6c551cfa`.
The source recomputes joint matrices at every retained parent and directly
recounts all reported strong endpoints. Counts, ordered covers, and signed
cut intermediates fit the documented 64-bit bounds through order16.

The separately written, no-producer-import
[`root_audit.py`](../../../../04-computation/nonzero_preparation_20261103/root_audit.py)
uses arbitrary-precision endpoint deletion and the direct cut sum. It
independently enumerates all2048 first cuts and all151552 second cuts,
recovering retained-parent counts37 and164 and all target multiplicities,
including the first-birth +4 core. Together with literal small one-exit
controls it passes **154,382 gates** in normal and optimized Python modes.
Both modes match
[`root_audit.json`](../../../../04-computation/nonzero_preparation_20261103/root_audit.json).
The source/receipt SHA256 pins are respectively
`68c809eeca40d312dfe38b5b116b31f4cf4c519b586b9ada47bf1d559e06e726`
and `9338c31ffa48971a8758d0eb17fa3d910e362a4672284801fc55b04740146b43`.

```text
c++ -std=c++17 -O2 04-computation/nonzero_preparation_20261103/induced_complete.cpp -o /tmp/nonzero_induced_audit
/tmp/nonzero_induced_audit > /tmp/nonzero_induced_audit.json
python3 -B 04-computation/nonzero_preparation_20261103/root_audit.py --output /tmp/nonzero_root_audit.json
python3 -B -O 04-computation/nonzero_preparation_20261103/root_audit.py --output /tmp/nonzero_root_audit_optimized.json
```

This is an all-order theorem built from a **FINITE-EXACT** certificate and a
**CITED** order bound. Status/link-only promotion of the proof does not change
the frozen source or receipt. No Lean formalization is claimed.
