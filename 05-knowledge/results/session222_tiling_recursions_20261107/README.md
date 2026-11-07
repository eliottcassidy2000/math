# Session222: tiling growth, doubling and the odd-chain plane

**COMPLETE: audited shell equivalence, physical doubling lifts, an odd-chain lattice traversal and a continuous plane-filling construction. Chan–Pak remains OPEN.** No literature-priority or unrestricted spectrum-completeness claim is made.

The user's picture has a rigorous realization once the changing quantities are kept visible. An L shell and two diagonal layers grow the same staircase. Odd chains and powers of two label a two-dimensional grid. Several real doubling laws occur, but they act on different objects.

## 1. The L and the two diagonals

The n-vertex fixed-path tiling carrier has binom(n−1,2) free cells. Both growth schemes add **2n−1 cells** and reach order n+2:

- Shift the old triangle by (1,1), adding the two-axis L.
- Keep the old triangle anchored, adding two sum-diagonal layers of lengths n−1 and n.

An explicit cell permutation identifies the two enlarged Boolean cubes. It does **not** preserve tournament incidence or Hamiltonian-path count H. The smallest decisive example is the directed triangle: zero L-growth preserves H=3 by adding a source and sink, while zero diagonal growth gives H=9. Thus vertex +2, cell +2n−1, and H-growth are separate recursions.

[Proof and recovered historical sources](shells/REPORT.md) · [L/diagonal diagram](shells/shells.svg) · [Independent shell audit](audit/SHELL_AUDIT.md).

## 2. Four genuine doubling mechanisms

- **A single tournament's defect:** inserting an interior vertex in the one-back-arc family gives H'=2H−1, so (H−1)/2 doubles. Its counts are3,5,9,17,33,... . Two diagonal layers quadruple H−1 in this family.
- **A physical count difference:** for every N>=1, a strong marked single-arc pair with gap2N has a one-vertex extension with gap4N. Both old tournaments survive as induced subtournaments, and the same old marked arc is the only disagreement. This nested lift of Session221's selector can be iterated. For example9/11 extends to19/23. The two individual H-values do not generally double.
- **Reflection-fixed information:** Mode B's n→n+2 grows the blue-code dimension by n, but the quotient B/B-perp gains exactly one bit. Its cardinality doubles. This is a literal +2→×2 connection on the stated quotient, not on H or the full tiling set.
- **Tournament order and spectral scale:** skew-Sylvester doubling sends n→2n and lambda²+1→2(lambda²+1). The singleton tower has H=1,1,3,189 at orders1,2,4,8, so it is not H-doubling.

[Proofs, constructors and examples](doubling/REPORT.md) · [Independent audit](audit/REPORT.md).

H itself is always odd. A universal H→2H operation is impossible; even H→2H+1 and H→2H−1 fail globally at3→7 and11→21. The family-specific laws above do not bypass those forbidden values.

## 3. The integers really do form an odd-chain grid

The prior decoder already uses the unique address

    n=2^k(2j+1),  j,k>=0.

Moving up doubles n. Moving right increases the odd core by2, so the actual integer increment is2^(k+1). The corresponding maps commute. Literal +2 on every row is different: its graph with doubling contains the five-cycle1–3–6–4–2–1, impossible to embed with all edges as unit edges in a square grid.

Alternating square shells of lengths1,3,5,7,... give an explicit nearest-neighbor traversal of every dyadic address. The integer labels begin

    1,2,6,3,5,10,20,12,4,8,24,40,56,28,14,7,... .

Every positive integer appears exactly once. The labels are not the visit times. The first(r+1)² visits fill the address square[0,r]²; the largest label there is2^r(2r+1), not(r+1)².

This grid has a native tournament realization recovered from earlier work: F(j,k)=R_(2j+1)[TT_(2^k)], with a fixed nested chain of odd regular cores R_q. Horizontal growth adds two whole fibers, while vertical growth substitutes TT2 into every vertex. These operations commute with compatible core/fiber labels. The construction supplies one chosen tournament at each vertex order, not every Hamiltonian-path count.

[Explicit maps, inverses and source verification](plane/REPORT.md) · [Odd-chain grid diagram](plane/oddplane.svg).

## 4. From a lattice traversal to a continuous plane-filling curve

A separate alternating-L spiral visits every vertex of Z² once. Its unit cells tile the entire plane. Transporting the odd-chain labels along this schedule preserves every scheduled path step, though not all edges of the original quadrant grid.

Joining the lattice centers with straight segments has area zero. To make the stronger continuous statement true, the report adds an explicit cell-filling rule. Interleaved Cantor digits give a continuous surjection onto a square after interpolation over deleted intervals. Connecting successive copies through shared-edge midpoints yields

    G:[0,infinity) → R²,   G([t,t+1]) = the complete cell visited at time t.

This is a genuine continuous surjection, proved without an external existence assumption. It is noninjective and adds structure beyond the original discrete path. A compact interval cannot continuously cover the unbounded whole plane.

[Self-contained construction and continuity proof](plane/REPORT.md#6-positive-continuous-realization-a-half-line-surjects-onto-the-plane) · [Independent continuity audit](audit/REPORT.md#5-subsequent-positive-continuous-extension-audit).

## 5. A precise remaining experimental question

Zero two-diagonal anchoring strictly increases H for all32,767 nonzero framed order-seven tilings and every smaller order checked. Three thousand seeded nonzero order-eight samples also increase it. These are **FINITE-EXACT** bounded checks; order eight is a sample, not an exhaustive census. The order-seven minimum positive increase is6.

**OPEN:** does every nonzero framed tiling strictly increase H under this operation at every order? The operation moves only the high endpoint of each backward arc; no common old-vertex embedding or all-order path injection has been established. The experiment gives a focused next target, not a monotonicity theorem.

## 6. Repairs and inheritance

The audit repairs THM-452's undeclared forward/backward bit gauge, its tower-order indexing, and its mistaken affine representative. The constant syndrome c is reflection-fixed, hence c lies in B; B+c=B cannot be the nonzero syndrome fiber. The correct fiber is t0+B for an actual doubled tiling t0. Nonblue statements require n>=2; the singleton double has no free tiles. The original block mechanism and properly typed forward-weight law survive.

[Corrected canon](../../../01-canon/theorems/THM-452-doubling-tiling-coset-law.md) · [Mechanism ledger](../../../01-canon/MISTAKES.md#2026-11-07-session222-repairs-the-doubling-tiling-gauge-and-affine-fiber).

Anchor: the user's tiling recursions. Niche: native odd-core/dyadic graph operations. Wildcard: positive continuous plane filling. The concept board kept cell geometry, vertex incidence, count versus difference, integer label versus visit time, and discrete versus continuous coverage separate. Inherited mechanisms are THM-291's source/sink extension, THM-477's blue code, the September25 odd-chain decoder and Session221's response selector. The principal hostile is H3→H9 under the wrong geometric embedding.

The concurrent root-change audit was read and preserved: no-singletons is not an all-K theorem because M_7(1897)=1. Neither a complete order grid nor a doubled response supplies absolute-count coverage. We used the existing meta-patterns “Type every analogy and every implication” and “Search the statement before the method”; no new universal methodology claim is promoted.

## Verification and publication

Run `python3 verify_all.py`. Only Python's standard library is needed. It checks all hashes; replays the author and independent programs normally and under `-O`; compares exact saved receipts; and keeps generated files in a temporary directory. All validation conditions stay active under optimization.

Independent universes:1,098 framed shell cases through order6;1,099 labelled parity controls through order5;1,100 framed doubling cases through order6 and tower through32;354 eligible marked selector locations through length8, including31 strong pairs; and10,000 vertices/100 complete shells in each lattice traversal. Author programs additionally check all255 apex words through length7,64 selector seed/lift examples, one-back-arc counts through order12, blue-code dimensions through order100,100,000 integer addresses and40,000 path positions. The proofs supply the all-order statements; the finite computations do not replace them.

[Checkpoint1](https://github.com/eliottcassidy2000/math/commit/e176f26bdd12eaa9944cb47e3b352de55683c2b6) published the audited shell packet; all15 files were read back exactly. The final checkpoint adds the audited doubling and plane constructions, repairs the canon, and completes this portable replay. Three concurrent commits through a6c471c45a4395bc985dbd603a9b278cef5adb02 were incorporated before publication without altering unrelated work.
