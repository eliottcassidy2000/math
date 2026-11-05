# Independent audit: fixed-left C2/C3 all-modulus families

Verdict: PASS, no correction required.

PROOF.md SHA256: c23f009564b6bb347062a182d24b5bc110cca482fe6f8754f7c2023761efddc7
verify.py SHA256: e6b304eea77d49287e09b79493d7d070fd6dfbc99b2fb58d011577ee47ff0dd3

The actual right source is the converse of the reviewed endpoint-one gadget,
not a formal profile. The fixed-left response rows give the stated formulas.
The Bezout selections preserve even L and odd g residues simultaneously for
M=2m, without assuming any coefficient is invertible. The inherited joint
compiler realizes those residues by actual binary words. Arbitrarily largeL
produce unbounded counts by the explicit linear lower bounds, so infinitely
many distinct counts in each desired residue class follow. The dyadic witness
order and general quadratic-in-modulus order bounds are consistent with the
inherited compiler. Neither is an exact-target bound.

Overlap check: session197/odd-modular-c2/PROOF.md asserts only odd moduli;
session197/dyadic-lowcycle/PROOF.md explicitly leaves mixed even moduli open.
The old joint-modular gadget theorem already gives all-modulus C1. This packet
correctly claims the direct fixed-left transfer to C2 and C3, not a new C1
statement or literature priority. Targeted searches found no older duplicate
of these two fixed-left all-modulus transfers.

Isolated normal and optimized replays pass identically, and all release
manifest hashes match afterward:63 physical source checks,126 actual graft
checks,2,348 modular checks,44 full integer replays. A separately authored
permutation counter additionally verified both families for L=2,3,4,5,
without importing the packet counter; see INDEPENDENT_LITERAL_AUDIT.json.

The assertion is congruence coverage only. In particular both families meet
19 modulo every modulus but cannot equal19, by their lower bounds29 and31.
There is no implication of cofinite or universal exact odd-target coverage.
