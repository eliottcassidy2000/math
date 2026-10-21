# Independent exception-tower transfer audit

Status: APPROVED mathematically and computationally. No required mathematical correction.

Scope: read-only inspection of math repository instructions, current tournament frontier, active guardrails, targeted mistakes, THM-3121-path-cover-walk-content-substitution-kernel.md and THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md. No repository or Git mutation. Reviewed PATHCOVERS.md, all four scripts and eight JSON receipts.

## Proof gates

- Ordinary transfer is correctly attributed to THM-3121 and translated from unordered pc_r to c_r=r! pc_r. Encounter-order run decoration is invertible, including cuts between two paths in the same block, and requires no SCC assumption.
- Fixed transfer retains exactly the missing coordinate. Reversal of the ordered path list makes even cover size impossible in odd order. The fixed block contributes d_(2s+1); each swapped pair contributes one ordered local c_r, not its square or a factorial. In the central quotient walk the unique central letter is the outer fixed vertex. The last-halfword arc into that vertex is necessary and correctly enforced by A e_0. Matrix multiplication orientation is correct.
- For C3, W and the specialized V=(1+Y+Y²)/(1-XY²) are correct. General anti-involution transfer follows from the same explicit decoration inverse, not from extrapolation of numbers.
- Full wreath normalization is properly restricted to canonical maximal strong modules. C3, R5 and P7 are prime outer templates here. Centralizer equation allows e(B) choices on the fixed block and a(B) per swapped block pair, and e(A) quotient choices. The exponents in the tower solutions agree with these recurrences.
- Oddness of tournament automorphism groups and Sylow conjugacy establish a/e anti-involutions. A fixed path orbit has one anti-involution stabilizer per path, yielding b=d1/e by double counting. Freeness gives h=H/a. The single-vertex seed may equivalently be handled directly or by a formal signed extension.
- THM-4731 Statement 3 explicitly identifies the fundamental black MULTIGRAPH degree at a self-converse class as (tc-bm)/2. Here tc=h and bm=b. Loops retain their degree multiplicity; this is not distinct-neighbor degree or adjacency reconstruction. The packet's caveat is necessary and correct.
- Triangle full c and d transfers use O(m^4) and O(m^3) arithmetic respectively. Coefficient dependency ordering is well founded; O(m^3) and O(m²) storage bounds hold. Fixed-template coefficient counts and cover counts have O(m log m) bits, so the stated polynomial bit-complexity conclusion is valid. No arbitrary-tournament polynomial counting claim is made.

## Exact reproduction

Replayed all four scripts normally and under python3 -O. All eight newly generated JSON outputs match the supplied outputs byte for byte. All entries of producer SHA256SUMS validate. Explicit exception guards remain active under -O.

The full tower vectors at 3 and 9 are independently controlled by vertex subset-cover DP and symmetric vertex-permutation/cut enumeration. Order-15 H is controlled by full vertex Held–Karp. Fixed H at 15, 21 and 27 is controlled by vertex pair-mask DP. At 21 and 45 the ordinary counts remain quotient-transfer results, not independent vertex-census results, exactly as scoped.

All displayed H, d1, h, b and paired degrees at 3,9,15,21,27,45 agree with replays; 81 and 243 receipt normalizations also agree exactly. The complete replay files are in this review directory.

Additional audit-owned hostile.py implements quotient words recursively and truncated polynomial convolution separately from producer's transfers. It compares complete ordinary cover vectors against vertex subset DP on C3[K1,T2,C3] (unequal block sizes), T3[C3,C3,C3] (off-family transitive outer), and T2[C3,T2]. It also compares the entire fixed vector of T3[C3] against symmetric vertex permutations/cuts. All pass normally and under -O; hostile receipts match.

## Nonblocking editorial clarifications

1. Section 5 says at 243 only H,d1 and normalizations are computed. The script actually computes full d243 transiently via fixed_triangle, retains only d1, and does not compute c243. Say “only H,d1 and normalizations are retained; the full c243 vector is not computed.” This does not affect complexity or values.
2. Cite THM-4731 plus its full slug explicitly for the inherited fundamental-black degree convention, rather than leaving the source implicit. Its statement 3 exactly supports the claimed convention.
3. PATHCOVERS.md is not included in the producer SHA256SUMS. The audit allowlist separately binds it, so this is not a publication blocker.

APPROVED_SHA256SUMS binds exactly the mathematical packet: document, four scripts, eight receipts and producer checksum manifest. Excludes caches and incidental stdout. If document wording is edited, obtain a new document hash; unchanged scripts and receipts remain approved.

## Prose repair verification

Producer implemented all three editorial clarifications. Verified the new explicit THM-4731 Statement 3 citation and loops-count-twice convention, the accurate transient full-d243 wording, and the refreshed producer manifest including PATHCOVERS.md. No code or numeric receipt changed. Approval is refreshed to the hashes in APPROVED_SHA256SUMS; document hash e18b58b2aa8b0da2714a770f106aa95f8a6cbae1081a9ffdcce08be3fb6ba97a.
