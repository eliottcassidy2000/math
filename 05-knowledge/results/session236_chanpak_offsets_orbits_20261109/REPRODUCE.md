# Reproduce Session236

From the repository root, using standard-library Python:

    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/last_letter/verify.py --rebuild
    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/orbit/verify_orbit.py
    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/replay.py
    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/orbit_audit/independent.py
    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets_audit/independent_check.py

The first checker also imports the inherited THM-5070 engine. The independent
matrix audit reconstructs the alphabet from physical tournament operations
without importing that engine. Expected counts: 57 matrices from 1,105 marked
blocks; four tower powers, 52 separated last letters and one excluded preimage;
615 cone inequalities; sharp growth292/193; 53,732 states at h<=3000, missing
10,31,61,106,217,274,337,346,607,865.

The two proof checkers use explicit verification gates and support --mutate;
mutation must fail both normally and under -O. FINAL_REPLAY.json records the
recovered final replay, and SHA256SUMS pins its packet. Finite enumeration is a
control; the all-order claims follow by mathematical induction.
