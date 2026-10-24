# Independent review of the neutral-fibre lifting theorem

23 October2026. **PASS** after proof and implementation review of
`../compatibility/NEUTRAL_FIBRE_LIFTING.md` and `../compatibility/check.py`.
Normal and optimized replays each passed661,673 checks after the final
14 direct final-C3-permutation completeness controls were added. The
theorem and previously audited interval/quotient logic were unchanged.

The interval compiler correctly uses up(F) and P\down(Z). Since the order
means that selecting x forces its successors, both boundaries are upsets.
The hypothesis that the old all-A/all-B assignments pass the pairwise
old-old constraints is essential: it restricts forbidden binary pairs to
01/10 and makes implication compilation exact. The note keeps that
hypothesis and does not encode the separate charge equation as2-SAT.

The induced neutral fibre is exact even across visible/neutral implications:
visible predecessors force their neutral successors through the lower
closure; neutral predecessors of a forbidden visible coordinate are removed
through the upper closure. Thus no cross-boundary implication is lost by
restricting to N. Removing minimal differing elements and then adding maximal
differing elements gives a lawful one-component route inside the interval.
Every neutral component must have zero charge AND cost at most the budget;
the text correctly keeps heavy neutral components visible.

For two neutral intervals the proposed witnesses

    Wa=La union (Lb intersection Ra),
    Wb=Lb union (La intersection Rb)

are legal upsets. Direct set subtraction gives Wa\Wb=La\Rb and
Wb\Wa=Lb\Ra, proving the exact weighted distance formula. This uses the
same neutral carrier and nonnegative coordinate-additive weights. Zero
weights may give zero distance between disjoint intervals; the proof
explicitly preserves that distinction.

The quotient theorem then follows in both directions. Actual edges pay at
least the interval distance. Conversely, connected light-neutral fibres
reach the attaining witnesses, and the simultaneous helper/visible/neutral
replacement has exactly the claimed cost. It is not necessary for the
intervals to intersect. The event bound (k+1)|N|+k follows by adjusting each
neutral coordinate at most once per intermediate or final fibre route.
This is budget connectivity, not equality of shortest-event distances.

Costs here are per-slot changes or owner-clock OCCURRENCES, so summing old
and helper costs is legitimate even when owners are shared. A metric that
counts distinct global owner variables only once would require another
cost theorem; the present statement does not claim it.

The native C3 control retains actual source owners, prefix pins, normalized
resource multiplicity, and the reflected converse clock2. Equality of the
first final-C3 owner forces all supplied slots into one common phase. The
third phase cannot avoid the full simultaneous-slot lower bound. The
all-size argument depends on the inherited endpoint-fixed block maps,
while the checker directly reconstructs r=1,2,3; this is correctly scoped.

No correction was required. The theorem supplies an exact compiler for an
authenticated helper supply; it does not assert that such a supply exists
for every charge fibre or every tournament order.
