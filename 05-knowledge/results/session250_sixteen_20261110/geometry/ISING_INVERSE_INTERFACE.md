# A concrete inverse interface: five-vertex switching is a signed cut-energy problem

Status: PROVED-IN-NOTE. This continues the boundary-degree theorem, not an assertion about the original 1-5-4-5-1 referent.

Let A=(a_ij) be the skew sign matrix of a five-vertex tournament. For i<j define

b_ij = sum_{k notin {i,j}} a_ik a_kj,

q_ij = sum_{(u,v,w) ordering [5] minus {i,j}} a_iu a_uv a_vw a_wj,

c_ij=(6b_ij+q_ij)/8.

Then for every switching sign vector s,

H(T^s)=15/2 + sum_{i<j} c_ij s_i s_j.                 (I)

Proof: the boundary theorem leaves only empty and two-endpoint supports. A surviving forest with endpoints i,j must be one even-length path. At order5 it has length2 or4. For a fixed length2 path, contraction leaves3 blocks; the two path directions and3! block arrangements contribute12/16=3/4 times its signed product. For a fixed length4 path, the two directions contribute2/16=1/8. Summing ordered internal vertices gives the stated formula. The constant is15/2.

Thus the spectrum **within this switching orbit** is exactly

{H(T)-2 sum_{i in U,j notin U} c_ij : U subset [5] modulo complement}.       (II)

This is an exact signed weighted-cut image, not an upper bound or analogy. No orientation is assigned to the symmetric c_ij. Missing switch-state labels can be restored by retaining U. General tournaments first require four-spin interactions at n=6; reducing those to pairwise energy without extra variables would lose information.

## Worked packet without a switch enumeration

For the regular cyclic tournament on Z/5, orient i->j when j-i is1 or2 modulo5. The c-weight is -1/2 on pentagon sides and2 on diagonals. At the unswitched state H=15.

Every switching cut has a unique representative of size0,1,2:

- Empty: H=15, one state.
- Singleton: cut weight3, H=9, five states.
- Two adjacent pentagon vertices: cut weight7, H=1, five states.
- Two nonadjacent pentagon vertices: cut weight2, H=11, five states.

This geometrically explains the packet {1:5,9:5,11:5,15:1}. In particular it is a different16-state decomposition from1-5-4-5-1, and must not be substituted for that referent.

## New inverse question and exact obstruction

A rational Clebsch eigenfunction with the right mean need not come from a tournament. The exact remaining inverse constraint is: do there exist ten skew signs a_ij in {+1,-1} satisfying the ten explicit equations8c_ij=6b_ij+q_ij? This is a finite quartic sign-feasibility interface that retains more than local adjacency means. Its constraints are necessary and sufficient, since a solution constructs a tournament and (I) recovers the entire switching function. The equivalence is constructive but is not a polynomial-time uniform algorithm or a classification.

The general-size successor is to use endpoint hyperedge coefficients from even path forests. A prospective finite-state Chan–Pak grammar could propagate this marked coefficient tensor under vertex extension; scalar H alone cannot do that. Whether a small closed sidecar exists is OPEN. For A000568 the sign feasibility still needs an isomorphism/stabilizer quotient; counting feasible coefficient tensors alone is not tournament enumeration.
