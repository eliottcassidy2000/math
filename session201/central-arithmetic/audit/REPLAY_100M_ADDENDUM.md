# Additional independently replayed finite inventory

**FINITE-EXACT: PASS through 100,000,001. Universal occupancy remains OPEN.**

After the main audit was frozen, the auditor ran the same reviewed `inventory.cpp` and independently reconstructed `central_forms.h` at the program's documented maximum cap. This is a new replay of the reviewed source, not adoption of the separately archived exploratory program or its output.

Inputs retain these SHA-256 hashes:

- `inventory.cpp`: `a2a5b627d22864f37de28ec32fde825d5b7169d6900b4b2ad4eacad45711798c`.
- `central_forms.h`: `105eb8f7f929560097d9c96b267615f004ef744f668b68e5afbce5f8ca703844`.

The completed output is `CPP_REPLAY_100000001.json`. It enumerates 19,856,206,889 index pairs. At maximum central lengths six, seven and eight it finds respectively 35, 31 and 29 gaps. Their lists are identical to the already replayed 10,000,001 lists. The largest length-eight gap remains 2647; the largest length-six/seven gap remains 17383.

The all-index product bounds, transpose-closed full menu, independently checked generated matrices, complete loop limits and explicit 64-bit safety bound were reviewed in `REPORT.md`. No claim about untested targets above 100,000,001 follows. No universal occupancy statement is promoted by this addendum, and no main proof or source code was changed.

From the packet root:

    g++ -O3 -std=c++17 inventory.cpp -o /tmp/central-inventory
    /tmp/central-inventory 100000001 > check-100m.json
    cmp check-100m.json audit/CPP_REPLAY_100000001.json

The run completed successfully in about four minutes in the audit environment. Elapsed time is descriptive, not a complexity bound.
