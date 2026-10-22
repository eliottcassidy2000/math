# Independent audit of the S51 packet (THM-4748, THM-4749, HYP-9359, note)

These are the independent auditor's own scripts and outputs, kept verbatim as audit evidence. None imports the session's code or the S50 audit helpers. `alib51.py` holds the shared helpers, and `c0_selftest.py` checks them first. Each script prints ALL CHECKS PASSED. The report is `REPORT.md`; the session's corrections are recorded in MISTAKE-657.

Run `c2` before `c8`, and `c4` before `c5`, `c5b`, `c6` and `c9`. They share caches in the session scratchpad, which is not committed.

| script | content |
|---|---|
| `c0_selftest.py` | helper self-test |
| `c1_two_arc.py` | THM-4748 (1)–(3): all 18,564 two-arc configurations with `m ≤ 16`; every tournament with `fa ≤ 2` on `m ≤ 16`; the mod-8 lines; no span-1 back arc in a minimum ordering |
| `c2_census.py` | the class census `n ≤ 9` (191,536 classes); least feedback arc numbers |
| `c3_three_arc.py` | THM-4748 (6): the three-arc census `m ≤ 13` |
| `c4_inventory.py` | an independent ear census and canonical form below 2100 (723,544 classes), compared with the session's inventory |
| `c5_support_law.py`, `c5b_generic_variants.py` | HYP-9359: the regression, the decomposition (exact covers, `3 | |Aut|`), R1–R7 and G1–G4, generic-part variants |
| `c6_degree_law.py` | THM-4749: the Session24 formula for every odd `h < 2100`, `n ≤ 30`; degrees, onsets, leading coefficients; the tables |
| `c7_chronology.py` | the pre-registration chronology: commit times against file modification times |
| `c8_refuted_leads.py` | the refuted leads, including the mod-8 transition mixing |
| `c9_fan_uniqueness.py` | THM-4748 (4) and the fan uniqueness at `k = 6, 9, 11` |
