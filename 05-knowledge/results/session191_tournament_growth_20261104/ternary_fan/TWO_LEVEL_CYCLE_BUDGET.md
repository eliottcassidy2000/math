# A complete physical two-level ternary fan, with sharp cycle budget seven

**PROVED construction; FINITE-EXACT cycle inventory with an all-order cutoff.** Start from N_k=9+11·2^k, k>=1. The cycle-crossing fan and the continuations below construct both ternary children and all four second-generation children through actual strong tournaments. Every intermediate has at most seven Hamiltonian cycles and lies below its own route's final path-count target. No third-generation or universal closure statement is made.

Seven is necessary within the specified fixed-core, same-clock, single-root wrapper architecture: the target family9N_k+14 requires exactly seven cycles in every such realization. This is not a lower bound for arbitrary tournament realizations of those scalar counts.

## 1. The missing lower-parent continuations

The first lower child from `CYCLE_CROSSING.md` has core

    A=(14,52,40,2,45,9), root1,
    (H,F_1,e_1)=(29,33,3).

Write L_k=29+33·2^k. Its lower child is obtained by two births: first add6 with sole old exit6->3; then add7 with sole old exit7->6. The cores and profiles are

    (78,116,104,2,109,73,8),
        (H,F_1,e_1)=(51,57,3),
    (206,244,232,130,237,201,8,64),
        (H,F_1,e_1)=(89,99,3).

For its upper child, the order of operations matters:

1. Add6 with sole old exit6->5, obtaining

       (78,116,104,66,109,9,32), profile(65,73,4).

2. Reverse the old arc3->1 to1->3, obtaining

       (78,124,104,64,109,9,32), profile(27,27,1).

3. Add7 beating exactly old vertices {2,3,4,6}, obtaining

       (206,252,104,64,109,137,32,92), profile(91,99,2).

The root remains the same physical vertex1 throughout. Every old exterior arc stays fixed. Each new core name is dominated by Z and w. All displayed cores are strong. Reversing3->1 before the first birth would make vertex3 a sink in the six-name core, so that unlisted ordering is not a strong route.

With x=2^k the lower route has counts

    29+33x ->51+57x ->89+99x=3L_k+2,

and cycle counts3,3,3. The upper route has

    29+33x ->65+73x ->27+27x ->91+99x=3L_k+4,

and cycle counts3,4,1,2. The strict differences between the upper final count and its two intermediate counts are26+26x and64+72x. The lower final count exceeds its intermediate by38+42x. Thus neither route overshoots its final target. The final outputs each have k+9 vertices.

These statements use the already proved identities H(W_k)=H(core)+2^kF and hc(W_k)=e_root(core). The finite core coordinates are certified independently by literal permutation/cut enumeration and subset recurrence in `verify_left_continuation.py`.

## 2. The complete two-level tree

Together with `REROOTED_CONTINUATION.md`, which supplies both children of the first upper parent, the exact tree is:

| Node | Entire path-count family | One chosen final cycle count |
|---|---|---:|
| Root | N_k=9+11·2^k | 1 |
| L | 3N_k+2=29+33·2^k | 3 |
| R | 3N_k+4=31+33·2^k | 3 |
| LL | 9N_k+8=89+99·2^k | 3 |
| LR | 9N_k+10=91+99·2^k | 2 |
| RL | 9N_k+14=95+99·2^k | 7 |
| RR | 9N_k+16=97+99·2^k | 2 |

Every edge of this abstract target tree now has a concrete physical route. Only the R-to-RL route changes the attachment root; it uses the two-in-neighbor bridge already proved in the re-rooting note. That bridge has six cycles, so all intermediates throughout the whole tree still have at most seven. The other routes have maxima at most four.

The assertion is stronger than seven unrelated scalar witnesses: each route retains all prior physical names, specifies its added vertices and changed old arcs, and verifies strongness at every stage. It is weaker than an iterating grammar: no outgoing routes from the four leaves have been supplied.

## 3. Exact sharpness of the cycle budget

Require each entire leaf family to be realized as W_k(D,r) for one fixed rooted core and the same clock k. Comparing constant and 2^k coefficients forces

    (H(D),F_r(D))=(89,99),(91,99),(95,99),(97,99),

respectively. The exact possible e_r(D), hence the exact possible leaf cycle counts, are

    H=89: {2,3,6},
    H=91: {2,6},
    H=95: {7},
    H=97: {2},

always with F_r=99. These are all-order results, not just an order-eight observation.

**Complete finite reduction.** Busch's bound excludes every strong core of order at least ten, because H>=5^3=125>97. All smaller cores are covered as follows.

- `generate_bank8.cpp` regenerates the complete 6,880-class order-eight bank by all vertex extensions and exhaustive degree-cell canonicalization. Its bytes match the independently audited inherited bank's SHA-256 55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903.
- All induced subsets give 33,278 distinct literal row lists of orders3 through8. These are not asserted to be isomorphism classes. Every smaller isomorphism class occurs.
- Exactly 701 order-eight representatives have H<=97, including both strong and nonstrong representatives. Extend each in all 256 ways, testing 179,456 physical order-nine matrices and every root.

This last reduction needs only Hamilton-path insertion monotonicity, not a strong-deletion theorem. Given any old Hamilton word and a new vertex, insert the vertex at the first legal position: first if it beats the first vertex; last if all old vertices beat it; otherwise between the first consecutive change from old->new to new->old. Deleting the new vertex recovers the old word, so H(parent)<=H(child). Thus every order-nine target has an order-eight parent among the 701 retained representatives. Relabeling that parent and selecting the corresponding extension realizes the target in the scan.

The order-nine scan contains 18,18,18,75 strong records with H=89,91,95,97 respectively; none contributes a new F=99 profile. These counts are extension records, not isomorphism-class counts. They are retained to distinguish a genuine nonempty scan from an omitted order.

**Nonstrong cores.** A strong wrapper forces its root into the final strong component C. Writing D=X->C gives F_r(D)=P_X(2)F_r(C). If X is nonempty, swapping the two ordered bins shows that P_X(2) is even. Since 99 is odd, D must already be strong. Thus nonstrong padding cannot evade this inventory.

The RL leaf consequently forces seven cycles in this architecture. The displayed physical tree achieves at most seven at every stage, so seven is sharp for this stated task. Other exteriors, clocks depending on the target, cores depending on k, or arbitrary graph realizations are outside the lower bound.

## 4. Reproduction and dependencies

    python3 verify_left_continuation.py
    python3 -O verify_left_continuation.py
    python3 verify_cycle_budget.py
    python3 -O verify_cycle_budget.py

The first script is standalone Python: six independent core profiles, 36 actual wrappers, 30 named transitions, and six direct lifted-cycle enumerations. The inventory script needs only the adjacent `generate_bank8.cpp` and `verify_cycle_budget.cpp`, a standard C++17 compiler, and Python's standard library. It generates every input in memory; no bank dump is needed. Both normal/optimized pairs produce identical stdout and receipts, with explicit checks active under `-O`.

The preceding two fan notes and the re-rooting note provide the other physical routes and their independent verifiers. Busch's lower bound is cited with the same source scope as those notes: A. H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC13(2006), N3, https://doi.org/10.37236/1141. No global spectrum theorem, full unbounded tournament classification, or priority claim against unexamined prior work is made.
