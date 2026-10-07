# Reader scripts, logs and witnesses: openai/math second reading (mac-mini-2026-10-06-oaimath2)

The six parallel readers of the session wrote and ran these scripts while reading openai/math manuscripts. The results are summarized and typed in [the results note](../../../05-knowledge/results/oai2_openai_math_second_reading_20261006.md).
* They are kept as run, so absolute scratch paths may need adjusting.
* No OpenAI code or LaTeX is included; the manuscripts are public at github.com/openai/math.
* The session's own independent re-checks are the four `../oai2_20261006_*.py` scripts.

| folder | manuscripts | main files |
|---|---|---|
| `plane_ramsey/` | #158, #172 | `s1_moser_certificate.py` (#158's placement certificate), `s2*_kappa*.py` (κ(q)), `s3*_residue_colorings.py` and `s4_multiquad.py` (Lemma R patch checks), `s5_hn_tower.py`, `s6_ramsey_idempotent.py` (THM-4559 certificates), `s7_fivecycle_words.py` |
| `hindman_ramsey/` | #164, #189 | `finv_cegar*.py`, `finv_direct.py` (the fully gap-determined Erdős 592 game), `finv_check.py` + `finv_witness_n3_t{4,5,6}.txt`, `finv_3_7*.log` (undecided runs), `davenport_level.py`, `smooth3_*.py`, `book_checks.py` |
| `barnette_cancellation/` | #180, #047 | `kp_enum.c`, `kp_certify.py`, `kp_minsets.py`, `kp_extend.py` (HYP-9215 data), `board_*paths.py` (8 × 8 board), `barnette_*.py`, `canc_verify*.py`, `canc_lnd.py` (#047 identities), `jc_fibres.py` (THM-4562 counts) |
| `littlewood_tiling/` | #076 (+2), #155 | `tower_merit.py`, `switched_tower.py`, `walsh_transfer.py`, `s_invariant.py`, `paley_rows.py`, `bent_check.py`, `bent_search.py`, `cyclo_partition.py` |
| `lattice_ds/` | #090, #022 | `c1..c11` (node sets, torus energies, the `N = 7` torus, deep-well overlaps, `u_A(N)` searches), `witnesses_D7.json` |
