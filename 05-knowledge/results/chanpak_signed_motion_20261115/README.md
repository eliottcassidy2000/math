# Signed count motion and interacting repairs

Research session, 2026-11-15. Chan–Pak remains open. Every displayed family here comes from explicit strong tournaments. We use selected physical witnesses and symbolic identities, not a larger numerical coverage bound.

## The main result: two blocks with the same common factor still give +2

Use two adjacent copies of the same five-vertex core, each with 13 Hamiltonian paths. An independent local repair has an extra factor of 13, but the coordinated move changes the two cross-contact terms as well. The three low changes are **+104, +15, −117**, leaving **+2**. The contact counts of both blocks change; the two *weighted combinations* read by the guards remain equal. At the smallest guards the counts are **39329 and 39331**, and the difference stays two for arbitrary independent positive left/right guard lengths. Only the control contacts change. Separating the two cores by even one transitive vertex destroys this particular unit.

[Mirror norm and repeatability test](MIRROR.md) explains this as the energy **2hJ + L² + T²** of a block and its converse. The same signed drift twice would require an isotropic direction modulo h; the exhibited drift (1,−1) cannot repeat at h=13. Adding a third copy of the same factor bank also fails with a single control, by unavoidable 13-divisibility. [Literal two-block fixtures](signed_move.md#two-active-nontransitive-blocks-cancel-the-amplification-or-bypass-it) include a second mechanism in which two amplified local moves are **−1422 and +1424**, again leaving +2.

## More than one repair; let even the quadratic response change

- **One fixed cyclic core, unchanged external tape:** move two arcs incident to the head per step, obtaining **4987→4989→4991** at the smallest height and the same two repairs at every symmetric height. On the second move both large contact counts decrease by 20; their signed contribution cancels. [Fixed-core route](multishift_cube.md).
- **Three repairs in one bounded material region:** keep the five-vertex core and both outer guards fixed, and move one unmarked pivot across the core between mask moves. This gives four consecutive odds **921,923,925,927** at the smallest height and four at every height. The fixed core contains only one nontransitive component. [Pivot route](signed_move.md).
- **Six consecutive odds with a changing core:** the length-six right-guard family changes all the raw count coordinates along the route while preserving the relevant projected coefficient. Coherent relabelings give at most seven arc flips on each +2 transition. [Six-count route](multishift_cube.md). [Synchronized-bank controls](signed_bezout.md) additionally give a five-count window that changes the quadratic response while using just one nontransitive component per state.
- **Grow inside the same strong core:** a port in one fixed six-vertex skeleton accepts any module P, giving a pair of readouts **312 Q3(P) + 394 Q2(P) − 17 H(P)** and **312 Q3(P) + 394 Q2(P) − 15 H(P)**. In particular, a growing transitive P gives an exact unit without parking the growth outside the core. The large endpoint and intersection terms themselves grow and cancel. [Internal growth and exact interpolation proof](MOTION.md).

For independent reproductions run each check script with Python from this directory. The C++ discovery programs cover the scoped bank searches, not the correctness of the infinite identities. No literature-priority claim is made.

## What is still missing

The shared-factor example proves that multiplication by neighboring counts need not kill the unit when the boundary is nonuniform. But the successful geometries do not yet supply an arbitrary-width window or choose the baseline from a prescribed odd target. The next search should be for a changing-drift bank on the same live seam: move h or select an isotropic direction in the mirror formula, or add a genuinely interacting second control, rather than tensoring together unmarked copies.
