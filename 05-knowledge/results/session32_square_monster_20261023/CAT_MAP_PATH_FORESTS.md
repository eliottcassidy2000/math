# A precise transfer-matrix bridge: square-sum path blocks and the cat map

Status: PROVED elementary identities below. The matrix interpretation is a real preserved-structure connection, not an explanation of the tournament Hamiltonian gaps or Monster prime orders. No novelty claim against the classical path-forest literature.

Let P_m be the undirected path on m vertices, m>=1, and let

f_m(t)=det(I+tL(P_m)).

By the matrix-forest theorem, this counts rooted spanning forests by number of edges. Endpoint-aware continuant expansion gives

f_1(t)=1, f_2(t)=1+2t,
f_m(t)=(1+2t)f_{m-1}(t)−t²f_{m-2}(t), m>=3.

The actual empty graph has forest polynomial1, so the recurrence must not be extended to m=2 by assigning its empty-graph value to f_0. A formal initial value f_0=0 is useful for the recurrence, but is not the empty-graph count.

At t=1, f_m(1)=F_{2m}, with F_0=0,F_1=1. The transfer matrix is

B=[[3,−1],[1,0]].

Let the usual cat map matrix be A=[[2,1],[1,1]], and U=[[1,−1],[0,1]]. Direct multiplication gives AU=UB, so B=U^{-1}AU over SL(2,Z). Thus the path-forest recurrence and the cat map have conjugate integral transfer operators, preserving eigenvalues, powers, determinants and traces. Both have characteristic polynomial λ²−3λ+1 and spectral radius φ².

Different observables must not be conflated:

- rooted forests on P_m: (B^{m−1})_11=F_{2m};
- cat-map power trace: tr(A^m)=L_{2m};
- any extra finite-sector trace in a VOA representation is additional representation data, not part of the path graph.

In particular tr(A²)=7, whereas f_4(1)=21. These are evaluations of different matrix observables, not a count-preserving identification. The source-reported logarithmic W(2) sequence3,7,21,47,123 uses an extra3 at powers divisible by3; importing it requires the separately cited FGST decomposition. The matrix identity alone neither imports that representation nor proves a Monster or Hamiltonian-gap statement.

## The user's 3,4,5 components

At n=12 the square-sum graph has components P_3,P_4,P_5. Hence its exact rooted-forest count is

F_6 F_8 F_10 =8·21·55=9240.

The component sizes3,4,5 are intrinsic graph data. They happen to satisfy a Pythagorean equation; no Euclidean metric or right-angle structure is encoded by that fact alone. The forest recurrence, in contrast, supplies an explicit integral operator and boundary functional.

At n=25 the user's three selected induced blocks have:

A={2,7,9,16,18}, a branched tree with4 edges;
B={12,4,5,11}, P_4 with3 edges;
C={10,6,3,1,8,17}, P_6 with5 edges.

Their rooted-forest polynomials are

F_A=1+8t+20t²+18t³+5t⁴,
F_B=1+6t+10t²+4t³,
F_C=1+10t+36t²+56t³+35t⁴+6t⁵.

At t=1 these give52,21,144. Replacing the4-edge branched block A by a4-edge path would give55 instead of52. Thus the edge counts4,3,5 lose real structure even for this elementary invariant. The full G25 also contains16 boundary edges and4 core edges, so its forest count is not the product of these three values.

## Hamiltonian-path boundary

The paths in these undirected blocks and paths in tournaments have different reversal behavior. An undirected Hamiltonian path reverses inside the same graph; a directed tournament path reverses in the converse tournament. Every tournament has a Hamiltonian path, but connected square-sum graphs can fail to have one (n=14,18–22,24). Consequently no map to tournaments can preserve Hamiltonian existence on this whole family. A labelled graph/tournament bit bijection does not fix this loss.

The useful shared carrier is the specified transfer or path-counting functional with its boundary information, not a forced identification of graph objects. The companion port-state note retains those boundary endpoints and connectivity exactly.

## Primary references and scope

The rooted-forest interpretation and Fibonacci/Lucas path formulas are classical: Pavel Chebotarev, *Spanning Forests and the Golden Ratio*, https://arxiv.org/abs/math/0612792 ; Chebotarev–Shamis, *Matrix-Forest Theorems*, https://arxiv.org/abs/math/0602575 . The integral conjugacy above is checked directly, not inferred from a shared numerical growth rate. The W(2) trace remark is a source-qualified pointer to the repository's separately audited FGST-based cat-map note; no general braided equivalence is asserted here.
