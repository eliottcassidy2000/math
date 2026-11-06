# Square-aware minimax bounds on the native order-four metagraph

Status: PROVED from finite exact graph/path certificates; independent audit PASS. No exact optimum or literature-priority claim. This strengthens both sides of the previously published order-four bounds without using the new preprint's logarithmic theorem.

## 1. Setting and square boundary information

Use the labelled path-preserving complement graph G of tournaments on four vertices, with lexicographic six-bit encoding, exactly as in session208/core-geometry. Remove its24 transitive leaves. The core has40 vertices and72 edges:48 edges A forming two subdivided cubes, and24 edges B forming six disjoint squares on the subdivision vertices. The48 vertex-relabel/converse automorphisms have exactly these two core-edge orbits. All original lengths and capacities are1 when defining the target optima S_* (maximum expected stretch) and C_* (maximum expected canonical fundamental-cut load).

Any core spanning tree T has39 edges. Let b be its number of B edges and s its number of squares in which it retains three edges. A tree retains j=0,1,2 or3 edges of each square; hence3s<=b. The number of selected A edges is39-b. Restoring the24 leaves changes neither objective optimum.

For each square Q, let J be its selected tree edges. Delete Q\J from the WHOLE core graph, leaving all other edges available. Every tree path is a path of this larger remaining graph. Thus the sum of tree distances for the four original square edges is at least the corresponding sum of shortest-path distances in that graph. This keeps boundary ports, unlike contracting the square.

Complete exact enumeration of all six squares and their15 proper edge subsets gives:

    retained j:                 0   1   2   3
    unit lengths:             20  16  12   6
    length4 onA, length5 onB:  84  68  52  30

For each j the sum is the same for every subset, including both possible shapes at j=2. The second row concerns auxiliary physical lengths used only in the dual certificate; original capacity inputs remain1. The shortest A-edge detour is5 in unit lengths and22 in auxiliary lengths4,5. check.py verifies every detour and all180 square-subset/length cases using integer shortest paths.

## 2. Improved distance lower bound

Let d_A,d_B be sums of unit-length tree distances over original A,B edges. Direct versus omitted A-edge counting gives

    d_A >=(39-b)+5(48-(39-b))=84+4b.

The first square table gives

    d_B >=120-4b-2s.

Give demand7 to each A edge and demand6 to each B edge, independently of unit lengths. Every spanning tree then satisfies

    7d_A+6d_B >=1308+4b-12s >=1308.

The total demand is7*48+6*24=480. Averaging over any tree distribution proves

    S_* >=1308/480=109/40.

This is an explicit weighted adversary for the paper's tree/edge game. It improves19/8 by retaining the four-edge square boundary instead of treating each missing edge separately. Equality would require strong simultaneous saturation; its existence is not asserted.

## 3. Improved congestion lower bound

Keep all graph capacities1, but charge a selected tree edge auxiliary length4 if it is in A and5 if it is in B. Let D_A,D_B be sums of these weighted tree distances over ORIGINAL A and B edges. The A detour certificate gives

    D_A >=4(39-b)+22(48-(39-b))=354+18b.

The second square table gives

    D_B >=504-16b-6s >=504-18b.

Therefore every spanning tree satisfies D_A+D_B>=858. By exact edge/path incidence double counting,

    D_A+D_B=4 sum_(f in A)L_T(f)+5 sum_(f in B)L_T(f),

where L_T(f) is its canonical fundamental-cut load if selected and0 otherwise. If every core edge's expected load is at most rho, the right-hand expectation is at most(4*48+5*24)rho=312rho. Hence

    C_* >=858/312=11/4.

This improves79/36. The proof uses demands on one side and physical lengths on the other, reflecting the paper's exact distance/load adjoint, without assuming its global tree theorem. No electrical-energy preservation is used.

## 4. Improved explicit upper certificates

WITNESS.json specifies four seed core spanning trees. The tuple below is (d_A,d_B,l_A,l_B), with l sums of unit-capacity canonical loads:

    distance seed1: (148,84,156,76),
    distance seed2: (172,52,132,92),
    capacity seed1: (148,84,158,74),
    capacity seed2: (172,52,134,90).

Their edge sets differ; identical distance sums do not imply identical load sums. Each seed has39 actual core edges and is connected. Uniformly apply one of the48 verified relabelling/converse automorphisms and restore all leaf edges.

For distance, mix seed1 with17/22 and seed2 with5/22. Every A-edge expected stretch is

    (17*148+5*172)/(22*48)=211/66,

and every B-edge expected stretch is

    (17*84+5*52)/(22*24)=211/66.

For capacity, use weight23/28 on seed1 and5/28 on seed2. Every A-edge expected load is

    (23*158+5*134)/(28*48)=269/84,

and every B-edge expected load is

    (23*74+5*90)/(28*24)=269/84.

Both distributions use at most96 trees. All leaves are exactly1. The checker directly verifies every edge's rational expectation across all192 seed-image trees, not only the aggregate orbit sums. Thus the current certified intervals are

    109/40 <=S_*<=211/66,
    11/4   <=C_*<=269/84.

Neither endpoint is claimed to equal the optimum. These improve267/83 and1063/330 while retaining the same support bound96; the older capacity construction with support48 remains a distinct option.

## 5. What the optimization attempt does and does not settle

Finite symmetry averaging reduces both distribution objectives to two coordinates. For any finite library of seed trees, its best rational orbit mixture can be found by two-dimensional convex hull arithmetic; a mixture of at most two orbit types suffices for this two-coordinate minimax problem. This does NOT certify that the library contains all globally useful trees. Global optimality requires a valid weighted-tree separation oracle or a universal matching lower bound.

A deterministic edge-exchange search found the four trees above. Search optimality is not asserted, and search reruns are not needed to verify the saved witnesses. A separate diagnostic multicommodity-flow MILP used72 binary edge-selection variables,10,368 directed commodity-flow variables,39 selected edges and edge-capacity coupling. The60-second solve timed out (302 branch nodes); its numeric gap and incumbent did not establish either optimum. No floating branch-and-bound lower bound is promoted here. The square table supplies the actual exact lower proof.

The supporting lines of the searched upper hull suggest candidate stronger universal inequalities4d_A+3d_B>=844 and2l_A+3l_B>=538. Both remain UNPROVED, not dependencies. The present work identifies the missing ingredient rather than declaring the small graph optimized.

Reproduce all claimed finite certificates: `python3 check.py` from any directory. Only the standard library and sibling WITNESS.json are used. No solver or floating computation is required. Sources/closest inherited mechanism: session208/core-geometry and core-improved; the main new coordinate is the selected subset of each square, not just its individual edge girths.
