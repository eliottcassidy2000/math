# Independent audit: named defect transport and two-birth census

2026-11-04. PASS. No mathematical correction required.

## Algebra

Reviewed the partition of marked two-slot covers into joinable and bad-join covers. Joining/splitting immediately after the physical mark is bijective, including the empty second slot. This proves F_w−H=sum_(v->w) J(w,v) exactly. The birth recurrence follows by expanding the named endpoint updates and taking the coefficient of the single new variable; the product of the two correction terms vanishes because z_x^2=0. In particular the c_w constant in K_w correctly records the newborn as a new in-neighbor. All right-hand data are from the old state; the next birth must use the updated state.

The proof has the correct scope: an exact unbounded data state, not a realization theorem for either affine branch, and not a legal native edge-complement compiler.

## Independent complete extension census

`independent_audit.cpp` is independently written and imports none of the author's code. It uses a permutation-tree DFS, not subset DP. A permutation with no backward adjacency contributes one to H and to every marked F. A permutation with one backward adjacency contributes one only to F at that bad adjacency's tail. Permutations with two backward adjacencies cannot contribute. This computes all eight marks together. Strongness is independently tested by Floyd-Warshall transitive closure.

It exhausts all 64*128 two-birth incidence choices from each literal six-vertex source, retaining every old row and physical name. Results:

- C2: six strong fixed-mark +2 hits, twelve strong all-mark +2 hits, zero +4 hits in either menu.
- C4: no strong hit for either digit at any output mark.
- Without strongness, exactly two fixed-mark C4 +2 hits: (57,127), global sink7; (63,57), global sink6.

Every hit and the complete totals are in `independent_audit.out`. The six fixed-mark cuts and added marks 2,6,7 agree exactly with the proof. This independently checks the negative census as well as the displayed positive witnesses.

Compile with the existing compiler:

g++ -O2 -std=c++17 -Wall -Wextra -Werror independent_audit.cpp -o /tmp/s181-grammar-independent
/tmp/s181-grammar-independent

## Normal and optimized controls

Replayed the author's permutation/polynomial verifier normally and with Python -O. Outputs byte-match. It checks 4,306 exhaustive named birth cases (all labelled orders 1..4, all existing marks, every cut), the defect identity, birth recurrence, six fixed-mark witnesses and twelve all-mark witnesses. Exceptions, not removable assertions, enforce checks.

## Provenance

Fetched the pinned S180 `phases/RESULT.md`: its exact C2/C4 rows, physical mark0 and responses (17,21)/(19,21) match this packet. Also inspected the named later `zero_ready_strong_switch_20261103/DEPTH.md` and README at commit 9f5bdfaea3a2debeaac8d076c0367d5006aca71a. The later actual-source results have different scope and are appropriately acknowledged. The bounded prototype failure does not contradict them.

The upper bound F<=8*8! used for uint64_t safety is conservative but valid. No extrapolation to large-order fixed-width arithmetic is made. No installations or forbidden-file access were used.
