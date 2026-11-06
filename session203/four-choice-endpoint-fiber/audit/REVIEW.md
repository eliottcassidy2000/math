# Independent audit of the four-choice wrapped endpoint family

6 November 2026. **PASS for the exact four-choice family at every nonnegative wrap count.** The number of certified choices remains four. No growing endpoint menu, iteration, interval coverage, or eventual odd-spectrum threshold is established.

## Reviewed inputs and physical graph

The frozen proof has SHA-256 2c6030f3f1dc2b990f99b38f40c473e1f147bf35c555a9450a35323267345315. The tensor producer has hash 6fe635693ea218abed22d38b3ad6e991db099d404bc9704e992d98b4f5092fec; its TENSORS.json has hash 5c0a8dc46574bfb7d9fe3d4c6d3cbe637fae77c04b9c8377c1d7e98270136484. The literal verifier has hash a6191e5d519c151a73766c6fd01e877a51e921a54c3be7245e5869a1cb5d4664. All seven entries of the input manifest pass and exact snapshots are retained.

The physical core is (102,0,10,83,135,30,182,47), retaining root 0 and replacing vertex 5 by a transitive spine. Only core vertex 3 has bit-dependent incidences: it beats a spine vertex exactly for bit 1. Every other incidence is inherited uniformly. The source words are 1^k followed by one of 01010,01110,10001,10101 and then 0^k. With the earliest spine bit placed least significantly, the seed integers are 10,14,17,21 and the full word is (seed<<k)+(2^k-1). The order is exactly 12+2k, and old vertex 1 remains a global sink. The construction never claims strong outputs.

Selective incidence choices are genuinely different from the earlier homogeneous-module switch. The proof derives their full colored context tensors directly, so it does not infer a selective-flip identity from the homogeneous result.

## Complete colored tensors and optional empty paths

The variables x,y,z,w count nonempty Hamilton paths with endpoint colors 00,01,10,11. Their products are subset convolution and the constant 1 is supported only on the empty body. Repeated variables remain labeled path slots.

The fixed-core gap bijection is valid for arbitrary colored bodies: each context path determines a unique core order and occupied gaps, and every permitted labeled body-path cover reconstructs that path. The rooted endpoint restrictions require both the appropriate first/last core vertex and an empty exterior gap. These restrictions are correctly included.

The producer uses merged-separator auxiliary directed cores for R,L,beta. This is legitimate, including optional empty paths, and does not invoke any tournament theorem on a bidirected graph. The independent checker instead partitions the actual core vertices into two ordered paths. It treats a bin with no core vertices as either empty or one nonempty colored body path, then imposes the root-start/root-end restrictions on the designated bin. Thus it does not reuse the auxiliary-graph construction.

Exhausting 5,760 core orders and 15,120 core partition cuts reproduces every coefficient of h,s,e,L,R,d,beta,B,V, including empty-body constants and all higher-degree zeros. The maximum degree is four. This exact finite-core identity applies at every body size; no polynomial interpolation in the spine length is used.

## Wrapping operators

A new global source of color 1 lies first in exactly one labeled path slot. If that slot has type 10, deleting the new source leaves a nonempty path of type 00 or 10. If it has type 11, deletion leaves an empty path or type 01 or 11. Choosing the labeled slot is recorded by the corresponding formal derivative. This proves

    D=(x+z) partial_z +(1+y+w) partial_w.

Similarly a new global sink of color 0 lies last in exactly one slot. Deletion from a 00 slot leaves empty, 00 or 01, and deletion from a 10 slot leaves 10 or 11. Hence

    E=(1+x+y) partial_x +(z+w) partial_z.

In particular D(1)=E(1)=0, as required because an empty cover cannot absorb a new vertex without an available slot. The constant in the coefficients is the empty-set unit, not a scalar function constantly equal to one. The operators are derivations on the cover algebra, not substitutions that multiply pointwise scalar path counts.

In coordinates Q=1+x+y+z+w, A=1+x+y, C=1+y+w, Y=y, direct computation gives

    D(Q)=Q, D(C)=C, D(A)=D(Y)=0,
    E(Q)=Q, E(A)=A, E(C)=E(Y)=0.

The affine coordinate change is invertible, so these generator identities prove D=Q partial_Q+C partial_C and E=Q partial_Q+A partial_A. The operators commute, both algebraically and because adding the exterior source/sink gives the same ordered body. A monomial Q^i A^j C^l Y^r is therefore an eigenvector for one complete wrap, with eigenvalue (i+l)(i+j).

For k=0 the iterated operator is the identity, including on eigenvectors of eigenvalue zero. For k>0 those zero-eigenvalue terms vanish. The table's delta_k is this wrap-count boundary term; it is different from the empty-set unit of the convolution algebra. Dropping it at k=0 would give incorrect h,R,beta,B,V values.

## Finite seed calculation proves the infinite family

Each induced subset of a transitive spine has one Hamilton path in its induced order. Thus Q is 1 on every subset; A allows empty subsets or first bit 0; C allows empty subsets or last bit 1; and Y allows only nonempty subsets with first bit 0 and last bit 1.

The independent checker evaluates diagonal monomials by explicitly assigning the five seed vertices to labeled slots, including empty Q/A/C slots and excluding empty Y slots. It does not use the producer's recursive convolution evaluator. It verifies the entire change-of-coordinate expansion and every grouped eigencomponent on every seed, including terms which evaluate to zero.

All eigencomponents of e,L,R,d,beta,B,V agree across the four seeds. For h and s, the only differences occur at eigenvalue 2 and are 0,8,14,22. Every baseline eigencomponent agrees with the displayed table in the proof. Therefore after k wraps the offsets are exactly {0,8,14,22}*2^k in both h and s, with the other five primary coordinates equal. B=h+L-s and V=R+beta consequently agree. The offsets are distinct for every integer k>=0.

This is a symbolic operator certificate: equality of all seed eigencomponents proves the statement for all k. The literal checks at small k are not used to extrapolate it.

## Independent controls

check_algebra.py verifies the diagonal generator formulas and commutation on all 70 monomials of degree at most four, the complete degree range of these tensors. It also verifies:

- 2,100 direct source/sink operator evaluations on all transitive colored words of length at most three.
- 10,500 direct source/sink operator evaluations on all 75 labeled colored tournaments through order three, including nontransitive bodies and the empty body. These use literal endpoint-path enumeration and explicit slot assignments.
- The complete original and diagonal tensor arrays, every seed eigencomponent, and the k=0 sum including eigenvalue zero.

The independent C++ checker constructs the spine after all seven retained core vertices rather than at the original module label. It counts paths by append-at-end subset DP with a separate root-start DP, then computes cover sums. This differs from the producer's remove-a-vertex start/end recurrence. It checks all nine coordinates directly, including beta, on all four words at k=0,1,2,3: sixteen actual tournaments through order eighteen. Normal and debug receipts agree byte-for-byte.

At k=0 the common row is (14829,0,152922,1371) and the starts are (1198,1206,1212,1220). At k=1 the row is (150887,0,2505086,13235) and starts are (12508,12524,12536,12552). The literal checks also reproduce the k=2 and k=3 rows and all offsets.

All finite counts fit 64-bit integers: ordinary Hamilton counts are at most 18!, and the two-bin quantities are bounded by 19!, below the signed 64-bit limit. Validation uses explicit failure checks in all builds. Both producer programs are replayed in their normal/optimized or optimized/debug modes and match their original receipts.

## Actual consequence and limitations

At each k this is one physical common row supporting four distinct starts. Increasing k scales the same four-point offset set while changing the row and the absolute starts. It does not increase the number of independently selectable endpoints. The common global sink remains present.

The readout 3B+2s, when invoked from its separately audited physical context, gives offsets {0,16,28,44}*2^k. They do not yield adjacent odd targets. The large starts do not themselves meet the bounded small endpoint requirements of F(8,137). There is no asserted operation for combining independently selected four-choice switches while holding the required row fixed.

## Reproduction

In this directory:

    python3 check_algebra.py > ALGEBRA_NORMAL.json
    python3 -O check_algebra.py > ALGEBRA_OPTIMIZED.json
    cmp ALGEBRA_NORMAL.json ALGEBRA_OPTIMIZED.json
    g++ -std=c++17 -O3 -Wall -Wextra check_graphs.cpp -o check_graphs
    ./check_graphs > GRAPHS_NORMAL.json
    g++ -std=c++17 -O0 -g -D_GLIBCXX_ASSERTIONS -Wall -Wextra check_graphs.cpp -o check_graphs_debug
    ./check_graphs_debug > GRAPHS_DEBUG.json
    cmp GRAPHS_NORMAL.json GRAPHS_DEBUG.json

The independent scripts use only the retained frozen tensor file and their explicit core data. Producer replays use PRODUCER_DERIVE.py and PRODUCER_VERIFIER.cpp. Publication and status promotion remain separate from this review.
