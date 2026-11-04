# Exact second-generation limits of the two-port ternary fan

**FINITE-EXACT profiles with an all-order cutoff.** The explicit fan in `PROOF.md` has parent N_k=5+7·2^k and children 3N_k+2 and 3N_k+4. At the next level, three of its four scalar branches have fixed-core realizations with the same clock. Only one has a unit-start wrapper realization. The fourth has no fixed-core same-clock realization at any order, including nonstrong cores that could yield a strong wrapper.

This is a statement about the exact wrapper architecture and coefficient matching. It is not a missing value of the global Hamilton-path spectrum.

## 1. Coefficient comparison

For a fixed core D and root r, use the same W_k construction and k for every member of the family. By the proven wrapper identity,

    H(W_k(D,r)) = H(D) + 2^k F_r(D).

The four twice-iterated ternary target families require exactly these profiles:

| Target | Required (H(D),F_r(D)) |
|---|---|
| 3(3N_k+2)+2 = 9N_k+8 | (53,63) |
| 3(3N_k+2)+4 = 9N_k+10 | (55,63) |
| 3(3N_k+4)+2 = 9N_k+14 | (59,63) |
| 3(3N_k+4)+4 = 9N_k+16 | (61,63) |

Comparing any two distinct k establishes these coefficient requirements; neither the graph order nor endpoint data can replace them.

## 2. Complete result, including the unit-start requirement

The exact inventory gives:

| Profile | Strong core exists? | Some W_k has a unit-start vertex? |
|---|---|---|
| (53,63) | Yes | Yes, with an eight-name core |
| (55,63) | Yes | No, for any core or k>=1 |
| (59,63) | Yes | No, for any core or k>=1 |
| (61,63) | No | No realization |

Explicit witnesses for the three existing profiles are:

- H=53, r=5, rows (2,16,3,7,45,79,159,63). Its complete J-vector is (1,1,2,4,1,9,27,18), so v=0 is a unit-start port of every wrapper.
- H=55, r=6, rows (2,24,11,33,13,87,31). Its J-vector is (3,3,7,5,5,9,31).
- H=59, r=5, rows (6,12,32,37,15,83,31). Its J-vector is (4,4,3,4,15,3,30).

The exact J-vector for each witness is generated in `second_receipt.json` and verified against literal permutation/cut enumeration. By the full start-vector formula from the fan proof, the old starts are J(v,r), the outer start is H(D), and the clock starts are powers of two times F_r(D). Since H(D),F_r(D)>1, a unit-start vertex exists if and only if some J(v,r)=1. The nonexistence statements test every rooted realization of the specified profile, not just these displayed examples.

### The surviving lower continuation can preserve the physical ports

The H=53 witness admits a relabeling to

    D=(250,244,249,226,8,16,48,112), r=1, v=3.

Its induced first six names are exactly the lower child core A from the fan. The two new names 6,7 are dominated by names 0,1,2,3, dominate names 4,5, and satisfy 7->6. Thus the first lower child really has a two-vertex extension realizing the next lower child while retaining the original physical root and unit-start port. This is one verified continuation, not a closed rule.

## 3. Why this finite computation is complete at all orders

For a strong tournament of order n, Busch's cited bound is

    H >= 5^((n-1)/3).

If n>=9, then H^3>=5^8=390625>61^3=226981. Thus every strong tournament with any of the four target H-values has order at most eight.

The inherited, independently audited bank contains all 6,880 order-eight isomorphism classes. Its exact SHA-256 is

    55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903.

Every smaller tournament occurs, up to isomorphism, as an induced subtournament of one of these representatives: extend the smaller tournament arbitrarily to eight vertices and then relabel to its bank representative. The verifier enumerates all nonempty induced subsets of size at least three and removes only literally identical row lists. This gives 33,278 input matrices; these are not claimed to be distinct isomorphism classes. For every strong input it computes H, every root's F, and, for the target profiles, every J-coordinate.

The retained profile-record counts at H=53,55,59,61 with F=63 are respectively 39,22,65,0. These count records in the induced-matrix input, not isomorphism classes. Exactly two records for H=53 have any J-coordinate equal to one; neither the H=55 nor H=59 records do. These results prove the strong-core assertions using the cutoff, with no large-order extrapolation.

### Nonstrong cores cannot evade the odd-F obstruction

W_k(D,r) is strong only if r belongs to the final strong component C of D. Necessity follows because that final component otherwise cannot reach r or w; sufficiency follows because every vertex of D then reaches r. Write D=X->C with r in C. The ordered-bin concatenation bijection gives

    F_r(D)=P_X(2) F_r(C),

where P_X(2) counts all ordered spanning two-bin covers, allowing empty bins. If X is nonempty, swapping the two bins is a fixed-point-free involution, so P_X(2) is even. Every target here has odd F=63. Therefore X must be empty and D must already be strong. This extends all exclusions above to every core that could produce a strong wrapper.

## 4. Reproduction and source scope

Run from the packet directory:

    python3 verify_second.py
    python3 -O verify_second.py

The script compiles `generate_bank8.cpp`, regenerates the complete bank with no data-file input, verifies its exact pinned hash, constructs the complete induced-subtournament input, compiles `verify_second.cpp` with standard C++17, and retains all target root records and full J-vectors. The bank generator starts from the empty tournament, tries every one-vertex extension, and canonicalizes by exhaustively permuting vertices within sorted outdegree cells. Isomorphisms preserve those cells, so this is exact canonicalization; induction proves complete coverage. It reproduces the inherited bank byte for byte in under one second in the tested environment. No bundled bank dump is required. A separate Python recurrence rechecks every retained target record and all retained H=61 matrices. Literal permutation/cut enumeration additionally verifies the three explicit witnesses. Both runs produce byte-identical output and the same receipt hash. The checks use explicit exceptions, not optimizable assertions.

The inherited receipt and review document the bank's independent canonicalization and labelled-orbit-mass completeness check. They do not replace independent review of this new profile scan. Source: [Session189 order-eight bank audit](https://github.com/eliottcassidy2000/math/blob/dcbe43e7fc9caf175baba4fef7b402fd9e3c51d9/05-knowledge/results/session189_response_frontier_20261104/prime8_audit/REVIEW.md). The lower bound is A. H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, https://doi.org/10.37236/1141. The source use is inherited from Session190's checked published abstract; no new independent verification of its full proof is claimed here.

The same receipt also computes the full rooted-F and total-two-bin value sets at H=61. Their possible use for arbitrary clock shifts is a separate extension and is not needed for the same-clock theorem in this note.
