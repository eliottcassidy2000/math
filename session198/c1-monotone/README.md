# Monotone C1 source inverse and exact coverage through50,001

`CERTIFICATE_PROOF.md` additionally proves NP membership for each strong C1/C2/C3 scalar language using O(log²N)-bit physical construction certificates. This is a verification bound, not a polynomial-time inverse.

`PROOF.md` proves that the physical C1 readout F=D−E+U strictly increases under every source move. It gives a complete target-sensitive inverse, retaining F<=N instead of the larger U<=N bank.

Run `python3 verify.py` and `python3 -O verify.py` from this directory. Both standard-library runs reproduce the same validation report and repair graphs. The50,001 run enumerates2,074,518 physical source profiles and independently replays24,971 positive-target source recipes. It does not need a network, external graph bank, or repository checkout.

C1 has exactly29 missing odd targets between3 and50,001, largest2,647. Explicit C2 witnesses fill26 of them, and a C3 witness fills67. All27 actual repair graphs are independently counted, with maximum order14. The resulting C1/C2/C3 union misses exactly7 and21 in this finite interval.

`REPAIR_WITNESSES.json` stores the source profiles and recipes for every positive exceptional repair. `REPAIR_GRAPHS.json` stores their actual adjacency rows. The full large positive-target recipe table is reconstructed by the verifier rather than uploaded; its canonical digest is in `VALIDATION.json`.

Universal and cofinite odd coverage remain open. The all-order theorem is the monotone inverse, not a claim that the finite missing list has stabilized forever.
