# Session52 MRB / supersingular / mod-15 research bundle

Read `RESEARCH.md` for proofs, exact scope, citations, and the connection/loss ledger.

Main result: A(8) determines a genuine CM elliptic curve of j-invariant 8000. Its supersingular-reduction filter is modulo 8, not an identification with the fifteen Ogg primes or an explanation of the actual MRB defect. Modulo 15, every unit residue has both an Ogg and a non-Ogg prime. The exact mod-15 MRB projection retains parity (period 30); without it the unit-filtered series diverges as (4/15) log²N. A certified inclusion-exclusion decomposition refutes weighting C by the unit density 8/15.

Run `python3 check.py` and `python3 -O check.py`. Both produce identical `results.json` and logs with 157,050 active checks. Run `python3 test_rejections.py` and its `-O` counterpart for eight explicit negative controls. No network is needed; Python standard library plus mpmath suffice. All code checks remain active under optimization.

The final packet also recovers the exact Bolza quotient from THM-4012, preserving the -1 twist and rejecting any JC transfer. A full-profile fifteen-point Fourier projection identifies the fourteen modes lost by the density approximation. `INHERITED_BOLZA_EXCERPT.md` pins the bounded prior proof.

`results.json` contains the finite arithmetic witnesses and directed intervals. `results.normal.json` is the normal-mode replay snapshot. `SOURCES.json` distinguishes external literature, current checkout, and archived Session31 proofs. `MANIFEST.json` hashes the deliverable files; `verify_manifest.py` checks them. Claims are typed separately as PROVED, CITED, FINITE-EXACT, VERIFIED-INTERVAL, or OPEN.

No git operation or publication was performed by this worker. The parent owns integration and independent review.
