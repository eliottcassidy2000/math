# Five coordinates are exactly necessary to track paths and cycles

Session 197. **All-order physical interface and exact continuation-rank obstruction.** This is a separate supplement to the frozen H=67 classification and rooted-start grammar. The full scalar spectrum remains open.

## 1. The returning state

Use the one-seam operation of `ROOTED.md`: all cross arcs A->G except b->a, and move the root to b. Begin with a rooted directed triangle, so every subsequent graph has order at least three. Keep

    x(A,a)=(h,e,R,s,C).

For the gadget constants defined there, the exact update is

    h' = g h - t e + L R,
    e' = E h - delta e,
    R' = Q h + B R,
    s' = d s,
    C' = E s.                                        (1)

The first three rows are the inherited S195 path interface. The final two are the separately proved physical gluing bijections in `ROOTED.md`. Thus every fixed actual gadget supplies a five-dimensional integer linear map. The outputs h and C are linear state coordinates. No order-two exception occurs, because the seed has order three and every step adds vertices.

## 2. Minimality for simultaneous path and cycle output

**Theorem.** On the one-seam grammar seeded at a directed triangle, permitting rooted gadgets of order at most four, no rational linear state of dimension at most four with fixed linear gadget maps and linear outputs for both H and C can represent all continuations. The five-state interface (1) is therefore minimal under these assumptions.

Use the triangle seed `([2,4,1],0)`, then take the following five actual inputs:

1. The seed itself: state `(3,1,3,1,1)`.
2. Append singleton `([0],0)`: state `(5,2,6,1,1)`.
3. Append transitive edge `([0,1],0)`, marked at its sink: `(9,3,9,1,1)`.
4. Append transitive triple `([0,1,3],0)`, marked at its sink: `(15,3,15,1,1)`.
5. Append rooted strong four-tournament `([8,1,3,6],1)`: `(41,6,54,3,2)`.

Each non-seed input is obtained by a single step from the seed. All five are strong. Use the following five continuation/output queries:

- H immediately;
- H after one singleton step;
- H after one sink-marked transitive-edge step;
- C immediately;
- C after one sink-marked transitive-edge step.

Their exact evaluation matrix, with rows indexing inputs, is

    [ 3   5   9  1  1 ]
    [ 5   9  17  1  1 ]
    [ 9  15  27  1  1 ]
    [15  27  45  1  1 ]
    [41  89 149  2  3 ]

Its determinant is -24. A d-dimensional rational linear realization factors this matrix into an input-state matrix with d columns and a continuation/output matrix with d rows. Its rank is at most d, so the nonzero determinant forces d>=5. The initial-state encoder need not itself be linear. Equation (1) supplies the matching upper bound. QED.

This is a minimality result for a specified linear realization problem, not a lower bound on every possible nonlinear algorithm or every graph invariant. It strengthens the inherited dimension-three result for H alone by tracking C as well.

## 3. The old cycle count does not substitute for the start response

There is also a direct physical collision, stronger than merely observing that a proposed formula has not been found. Consider

    A1=[32,17,11,19,37,14], root 3;
    A2=[32,17,11,35,13,22], root 2.

Both have

    (H,e,R,H(A-root),C)=(37,5,71,13,4).

Nevertheless their root-start counts are respectively 9 and 8. Appending the identical singleton seam gives

    H'=103, e'=32, R'=108

in both cases, but their new cycle counts are respectively 9 and 8.

Therefore the inherited four-coordinate path/deletion state `(H,e,R,H(A-root))`, even augmented with the current C, does not determine the next cycle count. This is a literal pair of graphs with equal proposed data and unequal outcomes. The new start coordinate s is not recoverable from those data on all rooted tournaments.

The collision does not contradict the dimension-five interface: its coordinates include s rather than deletion H. The section does not assert a minimal dimension when deletion H must also be an output.

## 4. Verification

Run

    python3 verify_cycle_interface.py > cycle.normal.stdout
    python3 -O verify_cycle_interface.py > cycle.optimized.stdout
    cmp cycle.normal.stdout cycle.optimized.stdout

The checker independently counts the collision inputs, collision outputs and five rank-witness inputs by literal vertex permutations and cut enumeration for R, and compares against the graph subset DP in `verify_rooted.py`. It counts the continuation graphs by subset DP and evaluates the determinant by its full signed permutation expansion, without numerical linear algebra or external packages. `cycle_interface_certificate.json` records every adjacency matrix and root.

The exact five-state transfer rests on the bijective proofs in `ROOTED.md`; the nonzero continuation minor proves the lower bound. The finite counterexamples are independently checkable and do not rely on a guessed global low-cycle classification.
