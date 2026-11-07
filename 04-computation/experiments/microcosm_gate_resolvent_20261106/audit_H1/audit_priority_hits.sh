#!/bin/sh
# audit_H1: prints the earlier in-repository appearances (file:line) of items that THM-4900 / HYP-9490 / the note present as new or leave uncredited.
cd "$(dirname "$0")/../../../.." || exit 1
show() { echo "--- $1"; shift; grep -n "$@" | cut -c1-260; }
echo "== phi(T) = F_T(1) is THM-1950's total inverse-response s(T) = 1^T(I+K)^(-1)1, K = A - A^T (u^T X u = u^T X^T u, so (I+M)^(-1) and (I-M)^(-1) give the same scalar)"
show "THM-1950 (klein-2026-07-21-S400)" -E "inverse-response|s = x|Velocity-addition|s\(T\)=n|regularity coordinate|A singleton has" 01-canon/theorems/THM-1950-h-ge-disc-reduced-to-strongly-connected.md
show "klein S400 reflection: strong minima of s for n = 3..8" -E "min s|0.667" 07-reflections/the-skew-determinant-composes-by-velocity-addition-klein-S400.md
show "THM-4051: order-7 strong range of the inverse response" -E "inverse response ranges" 01-canon/theorems/THM-4051-tournament-order-seven-strong-base-exact-frontier.md
show "THM-3729: det(I+K) u^T(I+K)^(-1)u = odd rooted-Pfaffian square sum" -E "det\(I\+K\)u\^T" 01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md
echo "== gate operator, all-n size-weighted cycle expansion (sw_1 = 2 c3, sw_2 = 2 c5 + 4 pairs), one-vertex tower values"
show "Session150 transfer (2026-11-02)" -E "Explicit two-vertex wrapper|T -> x -> y -> T" 05-knowledge/results/session150_hamiltonian_spectrum_20261103/transfer/TRANSFER.md
show "Session159 cycle_profiles (2026-11-03)" -E "F_T\(r\) = sum_K|r\^\(n−2\) coefficient|F_\{G\(T\)\}\(r\)|repeated wrappers" 05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md
show "Session159 profile_operator (2026-11-03)" -E "F_r\(G\(X\)\) =|P_G = \(z\+z|h\(X_m\)  |This transfer is entrywise nonnegative" 05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md
echo "== HC = 7 / HC = 21 class counts at the least orders"
show "Session212 cycle-levels n7_hist.tsv (cycles, classes)" -E "^7	" session212/cycle-levels/n7_hist.tsv
show "Session212 cycle-levels n8_hist.tsv (cycles, classes)" -E "^(7|21)	" session212/cycle-levels/n8_hist.tsv
echo "== F_21 = doubling plus translation mod 7; Paley primes; doubling tower keeps F_21"
show "THM-4556 remark" -E "F_21 = " 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
show "mersenne_switch_parity_f21 note (2026-10-06)" -E "generate .F_21|Modulo 7, .A_1" 05-knowledge/results/mersenne_switch_parity_f21_compression_20261006.md
show "seven_twentyone note (2026-10-06)" -E "F_21 = \{x" 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md
show "THM-4644 part 3 (Paley primes)" -E "Paley primes" 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md
show "THM-4557 / HYP-9162" -E "^title" 01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md
echo "== Session35 all-depth padding law (E-tower)"
show "Session35 section 6" -E "all-depth law" 05-knowledge/results/session35_dimensions_dualities_20261023/recursions/RESULT.md
