# Independent audit of the externally pinned native capacity obstruction

Status: **VERIFIED proof and complete finite checker audit.** Reviewed
`README.md` and `check.py`; normal and optimized replay both pass **72**
checks with exact frozen receipt equality. The count is not the census size:
the computation explicitly traverses all 2^21 normalized words, all 3,528
paired codes, and every stated pair or triple universe.

## Logical scope

At an independently owned source with one always-occupied paired unit,
preservation of every original dyadic source current under target replacement
is exactly equality of the paired target code. The factor two is essential.
The cited saturated code theorem identifies its complete target fibre as
{570,1965,2234}. This necessity would fail without the stated prohibition on
other variable incidences compensating at that same source row. The report
retains that condition, separates source and receiver roles, and does not
exclude more general multi-unit compensation.

The upper envelope is correctly one-sided. A source whose V choice can vary
must also reach U or W. All 120 such source classes are obtained from actual
normalized words; bitwise complement determines the actual partner target.
Therefore self/class-converse cases and possible loop options are included,
not silently removed by a simple graph representation. The independent
100-class U/V envelope matches the earlier 28+72 inventory.

Granting already occupied sources, repeated source labels, all presentations,
and no pins or resource exclusions can only enlarge the possible pair image.
Dropping the separate W-row equations also enlarges it. The remaining 60
missing charges are thus genuine necessary-condition failures. Residual 37
and fixed-leaf charge 1290 combine in Q to total 1327; this is not ordinary
integer addition of charge encodings.

The proof separately handles permanently V-only slots: they contribute fixed
constants, leaving at most one variable label, rather than expanding the
variable alphabet to all 381 V-incident classes. This distinction is required:
the deliberately excessive 381-class pair image is complete and would give a
false positive if incident edges were mistaken for switchable unit menus.

## Positive comparison and boundary

The 51 free helpers at occupancy three give all 2,048 residual charges by
complete enumeration of 20,825 triples. Translating by the fixed leaf charge
is correct. This produces helper occupancy five and is explicitly a separate
comparison family, not a legal occupancy-four repair. Its six-slot/43-clock
bound uses the six greatest actual switch costs, not a count of labels alone.

The grade-four core comparison retains both exact margins. Its 48-clock bound
is a valid upper bound using actual endpoint and hub selected weights; the
new native packet improves this by routing helper pairs through Johnson
exchanges. At high core grade, only core bits are complemented; helper charge
and the external owner pin remain unchanged. The proof states this sign
condition and makes no false full-system disconnectedness inference.

No correction requested. The finite obstruction, occupancy-five comparison,
and core-funded escape are compatible statements at different interfaces.

Replayed from the repository root:

    python -B 05-knowledge/results/external_pin_sections_20261024/coverage/check.py
    python -O -B 05-knowledge/results/external_pin_sections_20261024/coverage/check.py
