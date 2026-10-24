# Independent audit: fixed-occupancy native bank resilience

23 October2026. **PASS** for the mathematics and implementation in
`RESILIENCE.md`, `resilience.py`, and the decompressed full receipt.
Normal and optimized replays each passed527,396 checks. The complete source
sidecar and original core sidecar are pinned by normalized text hashes.

## Bounded hitting search

Any hitting set intersects the first witness. After choosing i there,
every remaining witness avoids i, and any additional hitter must meet the
first remaining witness. The common-intersection test is necessary and
sufficient for one final element to finish the hitting set. Accordingly:

- the initial intersection detects exactly the size-one case;
- branching on i then taking the residual intersection detects size two;
- branching also on j in the first residual witness then intersecting the
  still-unhit witnesses detects size three.

The size-three branch cannot wrongly reuse i or j: every remaining witness
excludes both. An empty remaining family would already have supplied a
size-two solution in the previous stage. Thus the early minimum searches
are complete and the returned masks have the stated sizes. The independent
literal enumeration on2,223 nonempty uniform hypergraphs checks this logic.
When no size<=3 solution exists, the output correctly records only a lower
bound4, with the greedy hitter as an upper-bound witness. It does not report
the greedy size as the minimum.

## Coverage and signed restrictions

Forced-U loss is the intersection-of-witnesses gate; forced-V loss is the
complement of their union. Both are independently replayed by scanning all
actual fixed-cardinality witnesses. The reduced48/k3 and28/k4 banks have
complete per-charge compact sections whose empty intersection proves every
single-U gate. Their inclusion-minimality is checked after greedy deletion;
no global cardinality minimum is inferred.

For the48/k4 bank, every one of the48 helpers can individually be forced to
either U or V while retaining all2,048 charges. Every charge also has hitting
number greater than3, proving the universal up-to-three-forced-U guarantee.
The calculation does not establish four-forced-U robustness. The53/k4 bank
has the same stated lower bound. I separately checked that the original
32/k4 bank is a subset of the53 pool, validating the optional monotonic
argument for that pool's forced-U half.

The forced-V half follows analytically from deletion-robust triple coverage:
fix h at V, solve q-label(h) with three distinct slots avoiding h, then adjoin
h. This is the forced-padding equivalence in `../controller/PADDING_PROBE.md`.
It does not rely on false unqualified monotonicity of subset-sum coverage.

For two forced-V helpers at occupancy4 there are only C(46,2)=1,035 or
C(51,2)=1,275 possible remaining pairs. Each is below2,048, so arbitrary
two-bit coverage is impossible. The stated mixed U/V obstruction is also
sound: with j required and i forbidden, the remaining triple must avoid
both i and j. The full receipt directly checks the translated missing
charge from the two-erasure triple hostile. Thus single arbitrary-bit
resilience is sharp even though one-sided U tolerance is larger.

## Moving restrictions and clocks

The dynamic result quantifies over core arrivals whose actual interactions
permit every fixed-occupancy helper choice except the listed restrictions.
It retains this substantive compatibility hypothesis. Choosing a different
certified section at each arrival remains valid; both helper subsets have
size4, so at most8 helper slots change, plus the two core slots. The stated
support10 and owner-clock-occurrence bound64 follow. This does not imply a
route of legal one-slot microsteps or a distinct-global-variable clock bound.
All helper slots remain occupied, and final prescribed helper records must
themselves be compatible with the final core state.

One wording precision was requested from the author: the phrase “recycling
a consumed resource” should mean simultaneous reuse or double counting at
an arrival, not prohibit later return to a released resource. The controller
explicitly uses the latter reversible-return contract. This does not affect
the mathematics or receipts; no producer file was changed by this audit.
