# Two grounded six-vertex pieces close the remaining finite grammar gaps

Session 195, 2026-11-05. **PROVED physical returning transfers and exact finite certificates.** The global odd Hamiltonian-path spectrum conjecture with exceptions `{7,21}` remains **OPEN**. This packet enlarges a physically composable grammar, not its seed list.

## 1. Starting point and scope

The named Session 195 `rooted-interfaces/PROOF.md` establishes a returning rooted one-seam interface and the complete bank of 65 matrices coming from all rooted tournaments of order at most five. The named `prime-growth/PROOF.md` integrates stationary-root near-sink and sink moves. Their combined grammar misses precisely `{7,21,67,173}` through 20,000. Its absence claims at a fixed count bound are exact over arbitrary composition depth, not just over graphs of a bounded order.

The present packet makes two changes: add two **actual six-vertex rooted tournaments** as reusable pieces, and retain exactly the same four physical coordinates and elementary operations. It independently verifies every constant, the resulting literal adjacency rows, the exact bounded spectrum, and a restricted minimality result. No priority claim or all-integer extrapolation is made.

## 2. Physical interface and operations

For a finite tournament `A` with marked vertex `a`, put

    x(A,a) = (h,e,R,d)
    h = H(A), with H(empty)=1;
    e = the number of Hamiltonian paths ending at a;
    R = the number of ordered spanning path covers (P,Q),
        where P may be empty and Q is nonempty and starts at a;
    d = H(A-a).

For a fixed rooted tournament `(G,b)`, orient all cross arcs from A to G except `b -> a`, and mark the resulting tournament at b. Write

    (g,s,L,E,Q,B,delta,dG)
    = (H(G), starts_b(G), ending_b_two_covers(G),
       ends_b(G), starting_b_two_covers(G), P_(G-b)(2),
       [|G|=1], H(G-b)).

Here `L` counts ordered spanning covers `(P,Q)` where P is nonempty and ends at b, while Q may be empty. `P_T(2)` permits both path slots to be empty. The exact returning transfer is

    h' = g h - s e + L R,
    e' = E h - delta e,
    R' = Q h + B R,
    d' = dG h.                                           (1)

For completeness, a Hamiltonian path avoiding `b->a` has one A run and one G run. There are `gh` candidates, with exactly `se` forbidden boundaries. A path using the seam has the uniquely determined sequence of runs `A_before,G_before,A_after,G_after`, with A_after starting at a and G_before ending at b. The other two runs may be empty, giving `LR`. A path ending at b cannot use the seam, giving `Eh-delta e`. In a rooted two-cover, an unused seam gives `Qh`; a used seam must be the first edge of its second path and gives `BR`. Deleting b removes the seam and leaves the ordered join `A -> (G-b)`, proving the fourth row. Thus these are genuine graph identities for **every rooted input A**, not numerical transitions fitted to the two target values.

The retained-root operations from the starting grammar are

    K_old(h,e,R,d) = (h+R-e, d, 2R-d, d),
    D(h,e,R,d)     = (h, 0, 2R, d).                       (2)

`K_old` adds one vertex beating only a and retains the old mark a; D adds a universal sink and retains a. The enlarged class starts only from a rooted singleton, and closes under (1) for the original 65-piece bank plus the two pieces below, and under (2). This is a recursive, grounded graph classification. A graph belongs exactly when it has a construction of this kind; no unique decomposition or classification of all tournaments is asserted.

## 3. The two new pieces

Rows are outgoing-neighbor bitmasks, with vertices numbered from zero.

### Piece A

    rows = [52,17,2,7,12,30],  root = 3,
    (g,s,L,E,Q,B,delta,dG) = (29,5,43,4,50,84,0,9).

It supplies the reusable operation

    A(h,e,R,d) = (29h-5e+43R, 4h, 50h+84R, 9h).          (3)

Applied directly to the singleton it gives

    (1,1,1,1) -> (67,4,134,9).

The literal seven-vertex output is

    rows = [110,104,34,4,15,24,60],  root = 4.

### Piece B

    rows = [56,1,3,6,14,30],  root = 1,
    (g,s,L,E,Q,B,delta,dG) = (23,1,75,8,24,84,0,9).

It supplies the reusable operation

    B(h,e,R,d) = (23h-e+75R, 8h, 24h+84R, 9h).           (4)

After one sink move from the singleton it gives

    (1,1,1,1) --D--> (1,0,2,1) --B--> (173,8,192,9).

The literal eight-vertex output is

    rows = [246,252,224,5,12,24,56,120],  root = 3.

These pieces are small structured objects rather than replacements of the current tournament. Deleting vertices `{0,1}` from A leaves a transitive tournament. Deleting vertex `{0}` from B leaves a transitive tournament. The larger state language still contains arbitrary iteration and unbounded graph order.

`gadgets.json` stores the pieces and constants. `witnesses.json` stores the two full output graphs. Independent subset DP and literal permutation/cut enumeration agree for all rooted profiles of both pieces and both outputs. The separate rooted-interface worker independently checked the same constants and both displayed output states.

## 4. Exact count-bound decision procedure

Every physical state satisfies

    1 <= d <= h <= R,    0 <= e <= h.

The inequality `d<=h` follows by inserting the deleted vertex into each Hamiltonian path: if it beats the first vertex insert it first; otherwise take the last vertex that beats it and insert immediately afterward. If a following vertex exists, it is beaten by the inserted vertex. Deletion recovers the original path, so this is injective. Splitting each Hamiltonian path just before its marked vertex injects into the R-covers, giving `R>=h`.

For each nonempty actual gadget,

    L>=g>=1, Q>=g>=1, 0<=s<=g, B>=1.

Consequently every root-moving step has `h'>=max(h,R)` and `R'>R`. The retained-root K step has `h'=h+R-e>=max(h,R)`. D preserves h and doubles R. All operations preserve nondecreasing h.

For a bound X, breadth-first search starts from `(1,1,1,1)`. For each operation it records its output count whenever `h'<=X`, and queues the resulting state only if `R'<=X` and it has not been visited already. There are finitely many states with `1<=d<=h<=R<=X` and `0<=e<=h`, so this terminates even when graph growth produces a state loop.

**Completeness.** In any successful construction of a count at most X, discard all terminal D steps: they do not change h. Before the final non-D step, every R must be at most X because the next non-D step has output h at least that R, and intervening D steps only increase R. Thus every required preterminal state is visited, and the final count is recorded even if its final R exceeds X. The procedure therefore decides exact arbitrary-depth membership up to X.

This proof is why absent small values in the certificates are true exclusions from the specified grammar, while the absence of new gaps in a long interval does not prove infinite coverage.

## 5. Exact finite spectrum and restricted minimality

With the original 65 matrices, the two new matrices, K_old, and D:

- Through 1,000: **21,353** retained four-coordinate states, with missing odds exactly `{7,21}`.
- Through 20,000: **2,584,754** retained four-coordinate states, with missing odds exactly `{7,21}`.

Thus every odd integer from 23 through 19,999 belongs to this single fixed physical grammar. These statements are not a proof for any integer above the tested bound.

The required maximum gadget order **to cover both 67 and 173, and hence the stated certified interval, is exactly six** within this operation architecture. Sufficiency is given by the explicit pieces. Necessity follows because the exhaustive bank containing every rooted gadget of order at most five still omits 67 and 173 even at arbitrary depth.

There is a sharper finite result: **at least two additional rooted gadget matrices are needed if all added gadgets have order at most six**. Exhaustive enumeration of all `2^15=32,768` labelled six-vertex tournaments and all six roots gives `196,608` rooted gadgets and **280 distinct full eight-constant signatures**. Deduplication includes dG, not merely the first seven constants. For each of these 280 possible additions separately, run the complete decision procedure through 173 with the original 65 matrices and K_old,D. No single addition attains both 67 and 173. The exhibited pair therefore attains this restricted minimum.

The full enumeration is independently reproduced in Python with `--full-inventory`; the C++ enumerator is an additional independent implementation. `single_addition_audit.json` records every candidate's bounded state count and which of the two targets it reaches. This result does not rule out a different operation family, a larger single gadget, or another non-matrix architecture.

## 6. Verification and remaining task

Run from this directory:

    python3 verify.py --full-inventory --bound 1000
    python3 -O verify.py --bound 1000
    g++ -O3 -std=c++17 -Wall -Wextra -pedantic census.cpp -o census
    ./census bank67.txt 20000

To regenerate the independent C++ six-vertex inventory:

    g++ -O3 -std=c++17 enumerate.cpp -o enumerate
    ./enumerate 6 > inventory6-regenerated.json

The Python checker uses explicit exceptions, so `python -O` does not remove its checks. It checks the representative constants by independent DP and literal enumeration, both full output graphs by those same independent methods, the new transfer laws on all rooted input graphs of order at most three, the exhaustive six-vertex inventory, every single-candidate minimality case, and the expanded grammar closure. The C++ census uses 64-bit arithmetic for transitions and 16-bit packed retained coordinates; it explicitly rejects bounds outside `1..65535` and validates physical state inequalities before packing. All retained coordinates are bounded by X, so there is no silent truncation.

The substantive open problem is now arithmetic rather than a missing finite gadget at 67 or 173: prove that this fixed recursively composable state language hits every odd target other than 7 and 21, or identify an additional obstruction. This work supplies no such theorem. Every update is still integer-linear, so a state whose four coordinates have a common divisor cannot escape that divisor. In particular, a universal claim that every marked parent has both literal `3H+2` and `3H+4` continuations remains false. An eventual global proof must choose suitable presentations or exploit another structural argument.
