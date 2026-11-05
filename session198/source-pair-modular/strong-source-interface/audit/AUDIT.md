# Independent audit: fixed-unit source interface

Verdict: PASS on frozen addendum bytes.

STRONG_SOURCE_INTERFACE.md SHA256: 10d7839eb16aab86397b05d19ffa36ad05040d9555ae01b56c5303b39fd76e95
check_interface.py SHA256: 0096f9f4f49e2ea9a65e3635063f833db8bb990cbebb9708cc3e44304a722e61

The literal source construction and all counting steps are valid. Ending-root
covers of T are subsets of its transitive deletion: either the empty subset,
or a subset whose last base vertex beats the root. This gives1+2(L-1)=2L-1.
Converse therefore has exact profile(g(L)+L-parity(L),1,2L-1,1), with both
root endpoint counts one, and it is strong C1 forL>=2. The singleton boundary
is correctly excluded.

The selection u=(U0+1)/2 is division of an odd representative plusone, not
an inverse of2 moduloM. With M even it fixes the required parity, and the
inherited joint(g,L) compiler yields the stated independent odd D,U residues.
The stronger cover congruence modulo2M follows directly. UnboundedL implies
unboundedD, so the infinite family assertion is substantive.

This strictly extends the old dyadic-only fixed-unit-interface conclusion;
the stronger old dyadic theorem about independent g,L,Q remains separate.
There is no claim of exact-value independence or unrestricted source profiles.

Isolated normal and optimized verifier replays agree byte-for-byte:127 physical
sources forL2..128,2,870 joint residue pairs over even moduli2..40,55 full
integer replays, and the excludedL1 boundary.
