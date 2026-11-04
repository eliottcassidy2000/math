# Independent audit of the named three-state gate packet

**PASS**, with the precise finite verification scope below. This review does not certify universal reconstruction or exclusion of the prime twelve-core residual. No full-response collision is claimed.

## Mathematical review

1. The wildcard-outgoing interpretation of E_v+sum_w C_vw is exact: each ignored successor identifies the unique ordered two-path cut, while the terminal case contributes E_v. The two insertions of x,y plus the x-first case exhaust the calibrated singleton query. Thus F(e_v,0)=H+2S_v is sound.
2. The twin-swap cancellation accounts for separated twins and both possible adjacent orders. The surviving block contracts to exactly the wildcard score of the quotient, giving a strictly positive increment. Consecutive ranks in any transitive module are ambient twins; their positive successive increments recover the entire actual physical order. This does not recognize modules or recover unrestricted tournaments.
3. The occupied-gap boundary certificates count every skeleton and every nonempty run assignment once. The maximum-three-run explanation applies to this boundary marked-permutation problem, not arbitrary two-port probes. The binomial recurrence solution has the correct index and terminal correction. In particular the new 4^(k-1) channel follows from the exact recurrence, not numerical curve fitting.
4. The six exterior-zero profiles imply the displayed fixed global-converse gaps after q_1=1, q_2=2^r-2 and q_3=3^r-3*2^r+3 are substituted. All gaps are positive even at r=1. Scalar H identifies state j, the exterior0 singleton identifies the external orientation, and module singleton rankings identify its arbitrary named transitive order. The described 6r! family is correctly bounded and all queries retain names.
5. The inherited deletion-cover coefficient applies to the physical first pair5->6, exterior i=0,j=3. Contracting this pair leaves the root first in the module. Freezing the cut boundary and forbidding its gap gives exactly the two stated cover profiles. Their difference(2,4,1) produces 2^m+2 for m>0 and2 at m=0, hence precisely the stated formula for m=r-2. The nonzero coefficient proves degree4; it does not imply a low-degree separator is impossible.
6. The whole-converse cubic and internal-reversal quartic are consistently separated. There is no claim that the complete full response equals only its H,S,E,C sidecar.

## Independent computations

`audit.py` imports no author modules. It independently constructs physical gate rows from pair orientations and uses a memoized predecessor Hamilton DP, distinct from the compiled/physical-query forward DP.

- All15 displayed finite skeleton profiles are recomputed by explicitly enumerating every permutation and every subset of occupied gaps, then testing the expanded word's actual arcs. This does not use the author's mandatory/permissible-set binomial shortcut.
- The exact unique packet universe is checked:21 whole-converse packets and18 pair-reversal packets, covering all3 states and the prescribed r ranges. It retains221208 paired-battery count slots; repeated base counts in different comparisons are not advertised as independent computations.
- Both values of five diverse physical battery assignments and all four dense query differences are independently recomputed for every packet.
- Every named module ordering at r=3, all3 states, and both exterior orientations is checked, giving36 full reconstruction-family controls. Every module singleton is evaluated literally; sorting recovers the actual order, and the scalar and exterior0 probes recover the stated choices.
- Every coefficient of all9 packets at r=1,2 is independently computed by literal Hamilton-word Boolean expansion. This retains signed factors for every old-port adjacency, rejects backward old-old and y->x steps, and compares the entire resulting named polynomial. It is not another interpolation of the supplied table.

The independent checker passes normally and under `python -O`, with byte-identical receipts. There are no removable assertion gates.

## Author replay and implementation checks

The author verifier passes normally and under optimization, and both replay outputs match the supplied receipt. It inverts every retained battery, checks the complete emitted coefficients, and evaluates its additional exact controls. The C++ generator and both Python DP implementations were source-reviewed. Full compiled regeneration of all large batteries is not repeated in this audit; fresh checks above independently cover selected physical queries in every large packet and every coefficient in the small packets.

Signed64-bit C++ arithmetic is safe here: the largest literal tournament has14 vertices, so each path-count entry is at most14!, and each four-variable finite difference sums at most16 such values, far below2^63. Named mask sizes are also within32-bit shifts.

A nonblocking hardening suggestion was sent to the author: explicitly check the complete unique (state,r,kind) packet universe and all displayed skeleton-profile constants in verify.py, rather than merely checking39 packets and emitting some profile values. The supplied data already pass these stronger independent checks; this is not a discovered mathematical or data failure. The author subsequently implemented both checks without changing the mathematical statements or data; its normal/optimized rerun is reported separately.

## Reproduction

    python audit.py > receipt.json
    python -O audit.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

Author replays are retained as author.normal.json and author.optimized.json. The independent receipt pins the full_response.json SHA256. All reads and writes were targeted; no forbidden ledger was read and no packages were installed.
