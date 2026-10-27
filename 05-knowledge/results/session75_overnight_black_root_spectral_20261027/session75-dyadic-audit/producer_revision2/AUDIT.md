# Self-audit and independent-audit handoff

## Status boundary

- The universal matching, insertion, local-bound and selected-tree claims have complete elementary proofs in README §§2–4. Finite checks are hostile controls, not their proof.
- All energies are raw squared bin sums. No normalization or conditional-expectation norm replaces the native R3 ledger.
- Multi-terminal energies are computed in an orthogonal direct sum, never by coherently summing across terminals.
- The selected tree is charge-dependent and changes the source partition. Its right-hand side is terminal energy, not source norm. It does not prove native `D^+<=S/2` or universal native entry.
- Generic sharpness and unbounded-fixed-tree families are not asserted to be prescribed native `A_N` configurations. The separate native local witnesses and 433/79 instance are constructed from the literal prescribed charges and native map.

## Exact computations

Normal and `python -O` runs both PASS and are byte-identical. The script uses explicit `check` exceptions, not `assert`. Runtime is about 1.3 seconds per replay in the current environment. It requires only Python's standard library and the included pinned frontier file.

Declared universes:

1. All 625 ordered four-tuples from `(-1,0),(0,0),(1,0),(0,1),(2,-1)`.
2. All perfect matchings for one explicitly generated integer-vector sample at each arity 2,4,6,8.
3. All 40320 permutations of eight labelled charges, four `+1` and four `−1`; these are actual nested block trees. Mean energies are `0,32/7,48/7,8`, their debt is 0, and mean realized debt is `8/7`.
4. All 256 length-eight sign words for the compatible minimum-matching construction, plus equal-charge sharpness controls.
5. Original native source sets at N=512,1024,4096. Complete first-arrival depths are 42,43,67. At N=4096 additionally depth cutoffs 1,3,10. No assertion is made about untested scales or all rectangular depth/clock cutoffs.
6. Exact native local witnesses `(N,q,u,r)=(512,1,435,0)` and `(4096,4,443,0)`, including their empty fourth/third bin, all original charge/clock/port/depth records, and actual low-bit pairing.
7. A same-terminal native pair at N512, sources799 and1065, with norm cross term0 and marked cross term6.

Independent arithmetic paths:

- Literal B-chain Eisenstein charge versus the three independent real source rows.
- Raw bin energy differences versus direct enumeration of every same-terminal labelled source pair and its valuation-defined LCA.
- Native middle energy versus explicit `(0,2)|(1,3)` grouping of four actual residue bins.
- Source-separated energy versus the sum of original source norms.
- Empty bins are implicit arithmetic zeros but explicit conceptual bins; the receipt records occupied-label and unoccupied counts. Zero-charge source labels remain in the generated source ledger and in pair enumeration.

## Recovered 433/79 instance

At N512 complete arrival: 768 sources, source norm1244, terminal energy1313, level energies

`1313,1274,1268,1261,1255,1242,1252,1247,1244,1244,1244`.

The deltas are `39,6,7,6,13,−10,5,3,0,0`, so positive whole-level debt is79. Positive parent-bin interactions total293; positive original labelled-pair contributions total433. This is a reproducible native instance matching the user history; identity with the unlocated historical computation is not claimed.

## Suggested independent audit pressure points

1. Check the exact insertion identity (3.1), especially P<C.
2. Check that the sharp local expected bound maximizes a convex function on the correct simplex, and that C>=P/4 is inherent.
3. In the global selected-tree proof, check the direction of the coarsening gain and the use `E>=P/m` for direct-sum terminal vectors.
4. Verify that the native witness parent bins contain exactly the claimed source labels, not a selected subset omitting other contributors.
5. Confirm the source guard `s !=5 mod8`, full clock sign, and full denominator with pending sources at partial cutoffs.
6. Check that the correlated three-pattern whole-level hostile is not misrepresented as excluding arbitrary independent parentwise choices.

## Reproduction

```bash
python3 verify.py > independently_replayed.json
python3 -O verify.py > independently_replayed.optimized.json
cmp independently_replayed.json independently_replayed.optimized.json
cmp independently_replayed.json receipt.json
```

`MANIFEST.json` freezes the milestone hashes. `SOURCE_PINS.json` authenticates the extracted canonical sources and frontier at commit `96e35642110104cc547516a63f37eebcf930f707`.

## Receipt serialization correction (superseding freeze v1)

Independent audit found that the original final receipt reused variable `r`
for a parent residue after loading source799, serializing the marked hostile
as `[0, source1065]`. The original source ledgers and the preceding numerical
check were correct; no mathematical claim changed. Distinct source/residue
variable names now prevent the overwrite. A regression validates the final
JSON-round-tripped pair records against complete expected records for799 and1065,
including source IDs, terminals, clocks, depths, entry ports, addresses and charges.
Normal and optimized receipts were regenerated and agree byte-for-byte.
All three full source-ledger SHA-256 hashes are unchanged.
`MANIFEST.superseded-v1.json` preserves the prior freeze; the new manifest binds
the corrected checker and receipts and explicitly supersedes their old hashes.
