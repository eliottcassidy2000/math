# An infinite four-choice endpoint fiber from selective bundle reversals

Status: **PROOF CANDIDATE, pending independent audit.** This constructs four distinct selectable endpoints at one physical padding row, for every scale k. It does not make the number of choices grow with k and does not prove interval coverage or an eventual odd-spectrum threshold.

## 1. Actual tournaments

Use outgoing-bit rows, vertices numbered from 0. Start with

    A=(102,0,10,83,135,30,182,47), root 0.

Replace vertex 5 by a transitive spine P. For every outside vertex other than vertex 3, keep its original uniform incidence toward the replacement module. For vertex 3, use a binary word w in the increasing spine order: vertex 3 beats a spine vertex exactly when its bit is 1. All remaining old edges are unchanged. Call the resulting rooted tournament U(w). Root 0 is never replaced. The old vertex 1 remains a global sink; strong outputs are not claimed.

For each k>=0 use the four words

    1^k 01010 0^k,
    1^k 01110 0^k,
    1^k 10001 0^k,
    1^k 10101 0^k.                                  (1)

Each output has 12+2k vertices. The five-letter middle words are encoded by integers 10,14,17,21 with the first spine bit in the least significant position. Thus the complete word integer is (seed << k) + (2^k-1), on 5+2k positions.

The underlying all-or-nothing module switch comes from the pair A and (142,0,34,22,71,155,47,94): relabeling the latter by (0,1,2,5,4,3,7,6) gives (102,0,10,115,135,22,182,47), differing only by reversal of edge 5->3. Here the individual spine incidences are chosen separately. The all-order statement below is proved directly, rather than inferred from that homogeneous switch.

## 2. Profile and theorem

For a rooted tournament write h=H(U), s=start_root(U), e=end_root(U), d=H(U-root). Let R count ordered spanning two-path covers with the second path starting at the root and the first possibly empty; let L count covers with the first path ending at the root and the second possibly empty. Let beta count ordered two-path covers of U-root with either side possibly empty. Set B=h+L-s and V=R+beta. The padding row is (B,e,V,d).

**Theorem.** For every k>=0, the four tournaments in (1) have the same values of L,R,e,d,beta, hence the same padding row. Relative to the first word, their h and s offsets are both

    {0,8,14,22} * 2^k.                              (2)

In particular they give four distinct endpoint choices at one actual row.

For a completely explicit description, put delta_k=1 if k=0 and 0 otherwise. The first word has the following coefficients in the ordered basis

    (delta_k, 1, 2^k, 4^k, 6^k, 9^k, 12^k, 16^k):

| Coordinate | delta_k | 1 | 2^k | 4^k | 6^k | 9^k | 12^k | 16^k |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h | 1 | 9 | -341 | 2480 | -6600 | 11664 | 0 | 0 |
| s | 0 | 8 | -110 | 456 | -1100 | 1944 | 0 | 0 |
| e | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| L | 0 | 16 | -462 | 3352 | -7700 | 13608 | 0 | 0 |
| R | 8 | -166 | 132 | 4832 | -26400 | 71244 | -91224 | 110592 |
| d | 0 | 11 | -132 | 648 | -1100 | 1944 | 0 | 0 |
| beta | 8 | -156 | -528 | 10240 | -39600 | 94572 | -91224 | 110592 |
| B | 1 | 17 | -693 | 5376 | -13200 | 23328 | 0 | 0 |
| V | 16 | -322 | -396 | 15072 | -66000 | 165816 | -182448 | 221184 |

For each other word add the corresponding offset (2) to h and s only. No rounded or asymptotic quantities occur.

At k=0 the common row is (14829,0,152922,1371), with starts 1198,1206,1212,1220. At k=1 it is (150887,0,2505086,13235), with starts 12508,12524,12536,12552.

## 3. Colored path-cover algebra

Color each spine vertex by its bit. More generally, allow an arbitrary colored body T. Let x,y,z,w count nonempty Hamilton paths of T with endpoint colors 00,01,10,11 respectively, and let 1 be the function supported only on the empty body. Products mean subset convolution:

    (fg)(T)=sum_(S subset V(T)) f(T[S]) g(T[V(T) minus S]).

Thus a monomial counts a labeled spanning path cover with specified endpoint colors. This is not ordinary multiplication of scalar counts. The empty-set unit is essential.

Each of h,s,e,L,R,d,beta for U(T) is an exact polynomial in x,y,z,w. To obtain it, order the fixed core vertices and choose occupied gaps; each occupied gap contains one nonempty body path with permissible endpoint colors. Empty interior gaps require actual forward core edges. Restricting the first/last core vertex and end gap gives s/e; deleting the root gives d.

For R use an auxiliary separator immediately before the root, merged into it: allow all vertices to enter the merged root while retaining its original outgoing edges. For L merge a separator immediately after it: retain original incoming edges and allow all outgoing edges. For beta delete the root and insert a bidirected universal separator. Splitting at the separator gives the specified ordered two-path covers. These possibly bidirected cores are counting devices only; physical U(T) is always a tournament.

The complete coefficient arrays are recomputed by derive.py and retained in TENSORS.json. The fixed-core gap bijection proves these identities for arbitrary colored T, not only the finite initial words.

## 4. Exact wrapping operators

Let D prepend a new global source of color 1 to the body. In a labeled path cover that vertex either forms a singleton color-11 path, or is prepended to one existing path ending in color 0 or 1. Hence on any colored cover polynomial F its exact action is the derivation

    D F = (x+z) partial_z F + (1+y+w) partial_w F.    (3)

Likewise, appending a global sink of color 0 acts by

    E F = (1+x+y) partial_x F + (z+w) partial_z F.    (4)

The two derivations commute, as expected since the source and sink additions commute as ordered extensions. Their coefficients include the empty-set unit, which handles a newly added vertex occupying its own path slot.

Introduce affine coordinates

    Q=1+x+y+z+w, A=1+x+y, C=1+y+w, Y=y.

The inverse change is

    x=A-1-Y, y=Y, z=Q-A-C+1+Y, w=C-1-Y.

Direct substitution in (3)-(4) gives

    D=Q partial_Q + C partial_C,
    E=Q partial_Q + A partial_A.

Thus one source-1/sink-0 wrap acts diagonally:

    DE(Q^i A^j C^l Y^r) = (i+l)(i+j) Q^i A^j C^l Y^r. (5)

After k wraps the multiplier is ((i+l)(i+j))^k. At k=0 this is 1 also for eigenvalue 0; at positive k a zero-eigenvalue term disappears. This accounts for delta_k in the table.

## 5. Finite eigencomponent certificate of the all-order family

For each profile polynomial, change coordinates as above and group monomials by the integer eigenvalue lambda=(i+l)(i+j). Evaluate each grouped polynomial on the four five-vertex colored seed spines. On any induced subset of a transitive colored spine:

- Q is 1, including the empty subset;
- A is 1 exactly for an empty subset or one whose first color is 0;
- C is 1 exactly for an empty subset or one whose last color is 1;
- Y is 1 exactly for a nonempty subset starting in 0 and ending in 1.

These rules give a finite exact subset-convolution calculation, including the optional empty slots. derive.py performs it for every eigencomponent and preserves both the original colored polynomial and its diagonal-coordinate expansion.

The complete certificate has the following conclusions:

- Every eigencomponent of e,L,R,d,beta,B,V agrees on all four seed bodies.
- For h and s, all eigencomponents except lambda=2 agree. The lambda=2 offsets are 0,8,14,22.
- The first seed's nonzero eigencomponent values are exactly the coefficient rows in the displayed table, with lambda=0 represented by delta_k.

Equation (5) now proves the claimed identities for every k. This is a finite certificate for a symbolic operator identity, not extrapolation from counts at a few k values.

A separate literal C++ subset-DP counter checks every profile coordinate for all four words at k=0,1,2,3, including 18-vertex outputs. Its 16 graph counts agree with the closed formulas, and optimized/unoptimized builds give identical receipts. These literal checks are supplementary controls.

## 6. Scope and next missing step

All four starts have the same parity. The ordinary 3B+2s readout, when applied, gives offsets {0,16,28,44}*2^k; this family alone does not provide adjacent odd targets. Nor does increasing k increase the number of independent endpoint choices: it scales the same four-point shape while changing its common base and other coordinates.

The exact difference direction Delta h=Delta s, with L,R,e,d,beta unchanged, may interact with neutralization constructions, but no new closure, a freely growing endpoint interval, or an all-order spectrum threshold is claimed here. The counts are large and do not directly supply the small endpoint menu in the conditional bounded-profile lemma.

The core tensor and wrapping derivations above are self-contained. Historical homogeneous-module and additive-readout packets motivate the construction but are not dependencies of its main theorem.

## Reproduction

    python3 derive.py > replay.tensors.json
    python3 -O derive.py > replay.tensors.optimized.json
    cmp replay.tensors.json TENSORS.json
    cmp replay.tensors.json replay.tensors.optimized.json
    g++ -O3 -std=c++17 check.cpp -o /tmp/four_choice
    /tmp/four_choice > replay.txt
    cmp replay.txt NORMAL.txt
    g++ -O0 -std=c++17 check.cpp -o /tmp/four_choice_debug
    /tmp/four_choice_debug > replay.debug.txt
    cmp replay.debug.txt DEBUG.txt

All validation branches remain active without Python or C++ assertions.
