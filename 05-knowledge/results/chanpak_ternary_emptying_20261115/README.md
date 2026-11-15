# Chan–Pak: endpoint reflection and the emptying invariant

**Universal coverage remains open.** The positive result here is a
shorter structural route to the known exclusions, with explicit
repairs and exact phase laws explaining why arithmetic alone
doesn't connect the two triadic chains.

1. **[Exact endpoint reflection](ENDPOINT_REFLECTION.md).**
   The induced-subtournament recursion gives
   T_v(a)=(-1)^(m+1) S_v(-a−1), with a counting complementary
   path slots. At one slot this gives mod 3; the optimal
   general unsigned modulus at fixed order is 2a+1.
   The next 3-adic digits count marked odd blocks emptied
   from the rooted complement. The
   [squarefree path proof](ROOTED_SLOT_RECIPROCITY.md) gives
   a second proof by canceling induced-set terms.
2. **[Sharp matched step](MATCHED_TERNARY_REPAIR.md).**
   Holding the relevant other root channels fixed genuinely
   forces three odd half-steps. Strong counts 327 and 333
   have the same continuation state and attain that bound.
   The [different two-tape mechanism](TRANSPORT_TERNARY_READOUT.md)
   has primitive repairs in all three half-step residues,
   with exact width-by-width counts. The step of three
   is an invariant of a specified readout, not the spectrum.
3. **[Empty to a single outlet](EMPTYING_EXCLUSIONS.md).**
   The exact formula H_new=H_old+2E accounts for the
   reopened odd cycles. Its small rooted table, together
   with a [Busch-free order-seven cutoff](SEVEN_CUTOFF.md)
   and SCC multiplication, re-proves the all-order
   impossibility of 7 and 21. For the known 63 witness the
   same emptying gives 25+2·19 from two coupled wings.
4. **[One-apex run law](RUN_EMPTYING.md).**
   The rotating mod-3 gap at fixed run count is exact: k
   one/zero run pairs omit k mod 3 and occupy the other two
   classes at every ternary precision. A fully symbolic
   argument excludes 7 and 21 in this model and gives precisely
   two converse words for 63. Strong 5, 15 and 45 have
   direct words; changing an outside guard moves the
   factorization from H to H+2 and realizes 63 through
   65=5·13. The family 2^(3r)−2^r+3 containing 63 has an
   exact 3-adic lifting law.
5. **[Triangle insertion](TRIANGLE_LIFT.md).** A literal
   replacement by three cyclic copies multiplies by 3
   only when the old vertex is a singleton SCC. Otherwise
   two and three visits give a positive excess. This is why
   the arithmetic chain is not an inverse rule in the
   strong part of the spectrum.

The useful constructive target exposed by the exclusions is
now explicit: for the desired odd N, choose an actual
rooted, preferably strongly connected remainder with
H+2E=N and the outlet in the first SCC, so it reaches the whole remainder.
The small cases show how to reject pairs (H,E), and the
wing identity shows how changing a guard creates an
otherwise unavailable pair. A general positive selector
for those rooted pairs is still missing.

## Reproduction

Run the five Python replay files with Python 3. The text
receipts give exact checks for the rooted deletions, induced
identities, physical matched repair, transport residue counts,
ternary lifts, and cycle-slot construction. Replays do not
replace the all-order proofs.
