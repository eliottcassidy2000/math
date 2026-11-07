# Eleven adjacent odd counts from common module coefficients

**PROVED by explicit finite certificates and an all-order bijection; independently audited.** Chan–Pak completeness remains open. This gives a uniform eleven-value interval family, not an enlarged unrestricted coverage range and not a single-edge-flip staircase. No novelty claim for individual scalar counts is made.

## Construction theorem

For each odd d from −7 through 13 there is an explicit strong rooted tournament (G_d,v_d), of order at most eight, such that for every nonempty tournament U,

    H(G_d[U at v_d]) = 36 Q_U(2) + d H(U).

Here H is the number of directed Hamilton paths, and

    Q_U(2) = sum over S subset V(U) of H(U[S]) H(U[V(U) minus S]),
    H(empty) = 1.

The substitution replaces the marked vertex by U, retaining uniform incidences toward every outside vertex. In particular for every integer k>=1 and every odd d in [−7,13] there is a **strong tournament of order at most k+7** with

    H = 36*2^k + d.

Thus all eleven consecutive odd integers from 36*2^k−7 through 36*2^k+13 occur. The construction selects one of eleven contexts; it does not assert that consecutive outputs differ in one edge.

## Literal cores

Rows are outgoing adjacency bitmasks, vertices numbered from zero. The first column a equals d+72. The last column is a directed Hamiltonian cycle, including its closing edge to the first vertex.

| a | root | outgoing rows | cycle |
|---|---|---|---|
| 65 | 0 | 4,5,88,19,35,15,59 | 0,2,6,3,4,5,1 |
| 67 | 5 | 2,12,25,97,11,23,55 | 0,1,2,3,6,5,4 |
| 69 | 0 | 2,96,11,19,7,29,61 | 0,1,6,5,2,3,4 |
| 71 | 4 | 2,4,9,163,15,87,31,119 | 0,1,2,3,7,5,6,4 |
| 73 | 3 | 2,4,161,7,15,91,31,123 | 0,1,2,7,5,6,4,3 |
| 75 | 0 | 2,12,9,113,39,71,23 | 0,1,2,3,4,5,6 |
| 77 | 0 | 2,12,25,81,35,15,55 | 0,1,2,3,6,4,5 |
| 79 | 3 | 2,4,97,7,15,91,155,63 | 0,1,2,5,6,7,4,3 |
| 81 | 4 | 2,4,9,99,15,87,151,63 | 0,1,2,3,5,6,7,4 |
| 83 | 3 | 2,4,49,7,43,203,31,95 | 0,1,2,4,5,7,6,3 |
| 85 | 0 | 2,8,3,53,71,151,47,95 | 0,1,3,4,6,5,7,2 |

The machine-readable copy is interval_witnesses.json. Each row list is a tournament, and each displayed cycle is valid, so each core is strong. Replacing a vertex of a strong tournament by any nonempty tournament remains strong: paths between module vertices can leave the module, follow a directed cycle through the hole in the quotient, and re-enter; all other reachability follows from quotient paths. This proves strongness even if U itself is not strong.

## The exact coefficient certificate

Let F_U(j) count ordered spanning covers of U by j nonempty directed paths. Thus F_U(1)=H(U) and F_U(2)=Q_U(2)−2H(U). For a fixed core with hole v, remove v and choose any permutation p of the other vertices. Its linear gaps include the initial and final gaps. A gap may be empty if it is at an end or its adjacent core vertices form a forward edge. A gap may be filled by a nonempty U-path if its preceding core vertex beats v and v beats its following core vertex; omit a condition when the corresponding endpoint is absent.

For each permutation multiply e_g+f_g z over all gaps, where e_g,f_g are these zero-one permissions, then sum over all permutations. Let c_j be the coefficients. Every Hamiltonian path of the substitution has a unique core order, occupied-gap set, and ordered U-path cover. Conversely these data produce one path. Therefore

    H(G[U at v]) = sum_j c_j F_U(j).

Every listed core has the exact coefficient certificate

    c_1=a, c_2=36, c_j=0 for every j>=3.

The value c_0=H(G-v) is recorded in VERIFY.json but does not contribute for nonempty U. The finite permutations exhaust each fixed core, so the certificate establishes an identity for arbitrary U, not interpolation from finitely many input sizes. Substituting F_U(2)=Q_U(2)−2H(U) proves the first formula.

If U is the transitive k-tournament, each ordered bipartition S,S-complement has exactly one path in each part, so Q_U(2)=2^k and H(U)=1. This gives the interval family.

## Verification and controls

Run from this directory:

    python3 verify.py > VERIFY.json
    python3 -O verify.py > VERIFY_OPT.json
    cmp VERIFY.json VERIFY_OPT.json

The verifier recomputes all core certificates, validates the literal cycles, and independently uses ordinary vertex-subset path dynamic programming on 891 literal composites: for each of eleven cores, all 75 labeled nonempty tournament seeds through order four and transitive seeds of sizes one through six. Thus the literal check does not use the claimed module formula to count an output. It also checks the exact eleven-value interval of core counts. It uses explicit exceptions rather than removable assertions.

Hostile scope controls: empty U is excluded from the theorem and generally does not follow the displayed formula. Negative spectral coefficients d are allowed but all underlying cover coefficients a and 36 are positive. The result does not supply all integers between the separated dyadic windows, nor an induction on arbitrary odd counts.

## Inheritance and provenance

The general path-cover substitution kernel is inherited from THM-3121-path-cover-walk-content-substitution-kernel.md. The positive core-gap expansion is already used in session203/transitive-module-rebuilding/PROOF.md and session203/arbitrary-module-rewrite/PROOF.md; Session215's multi-exit and fixed-caps packets explain the Stirling/evaluation connection. These methods and arbitrary-module lifting are not claimed new here.

The contribution is the literal common-coefficient bank giving eleven consecutive offsets, with strong outputs and full arbitrary-input scope. Candidate cores were found by scanning the existing inherited_tournaments_through8.txt in session194/weak-parent-certificate. The existence theorem does not need that bank to be complete: every necessary witness and every coefficient is included here and independently replayable. Exploratory scans and their restricted nonexistence findings are separate from this theorem.
