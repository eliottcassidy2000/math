# Integral directed-cut carrier for tournament ears

Status: proved elementary reduction, targeted finite check, and an obstruction; this does not prove the tournament Hamiltonian-path spectrum conjecture. It does not identify tournaments by their H value or make any claim of a numerical correspondence outside tournaments.

## Definitions and prior results

Let T be a tournament on V, |V|=n>=2, and let H be its directed Hamiltonian-path count. For S⊆V, define the child T+x_S by **x→i if and only if i∈S** and F(S)=H(T+x_S). Then F(∅)=F(V)=H. Set A_i=F({i})−H and K_ij=F({i})+F({j})−F({i,j})−H.

THM-4097 gives the exact exposed-slot formula F(S)=sum_(j in S) start_j + sum_(i outside S) end_i + sum_(i outside S,j in S) Q_ij, where Q_ij counts ordered two-path covers with boundary i,j. Its endpoint terms matter: deleting x from a child can break the parent path. THM-4118 gives F(S)=H+sum_(i in S) A_i−sum_({i,j}⊆S)K_ij with K_ij=Q_ij+Q_ji≥0. Rédei's parity theorem implies A_i and K_ij are even integers. THM-4115 already proves F(S)≥H for all S and gives the symmetric weighted cut plus field. The following flow reduction packages this known information as an integral directed cut with no field. It is not a claimed new global spectrum theorem.

## Integral directed-cut theorem

There are nonnegative integers c_ij (i≠j), computable by an integral max flow, such that for every S⊆V,

    F(S) = H + 2 sum_(i in S,j outside S) c_ij,
    c_ij+c_ji = K_ij/2,
    sum_(j≠i)c_ij = A_i/2.

The matrix c is generally nonunique and its orientation need not follow the arcs of T. When T is strongly connected of order at least three, its positive-support directed graph c is strongly connected.

Proof. Define a_i=A_i/2 and k_ij=K_ij/2. These are nonnegative integers: A_i≥0 by THM-4115 and K_ij≥0 from Q. The equation at S=V yields sum_i a_i=sum_(i<j) k_ij. For every S,

    sum_(i∈S) a_i − sum_({i,j}⊆S) k_ij = (F(S)−H)/2 ≥ 0.   (*)

Make a source, one node for each unordered pair ij, one node for each vertex i, and a sink. Give source→ij capacity k_ij, each pair→its endpoint capacity L=sum k_ij (or any larger integer), and vertex i→sink capacity a_i. Standard integral max-flow/min-cut, equivalently Hall's condition for capacitated edge-unit allocation, says all pair supply can reach the sink exactly when, for each U⊆V, the pair weights wholly inside U do not exceed the endpoint capacities in U. These are precisely (*). A full flow has each vertex sink-edge saturated. For pair ij, set c_ij to the flow from pair ij to vertex i and c_ji to its flow to vertex j. Integrality gives integer weights and the two displayed conditions. Expanding the cut sum yields sum_(i∈S)a_i−sum_({i,j}⊆S)k_ij, as required.

For strong T, every proper nonempty S has F(S)>H by the next lemma, so its c-cut is positive. A finite directed graph has positive weight on every nontrivial outgoing cut exactly when its positive-support graph is strongly connected.

## Strong-parent growth lemma

For strong T of order n≥3 and any nonempty proper S,

    F(S) ≥ H + 2 ceil((n−1)/2).

Proof. Along an old Hamiltonian path P, write b_i=1 when its ith vertex belongs to S. Inserting x before the first vertex is legal iff b_1=1; after the last iff b_n=0; and at an interior slot iff the adjacent bits are 0,1. Hence the number of valid insertions into P is

    b_1 + (1−b_n) + #(0→1 transitions) = 1 + #(1→0 transitions).

Children obtained from distinct pairs (P,slot) are distinct, since deleting x recovers P and the slot. These children form only a subset of all child Hamiltonian paths, because deletion may create an invalid old adjacency. A strong tournament has a directed Hamiltonian cycle. Its n rotations give n distinct old Hamiltonian paths. For nontrivial S, that cycle has an edge from S to its complement; this edge appears in n−1 rotations. Thus the sum of the extra 1→0 transitions over all old Hamiltonian paths is at least n−1. It follows F(S)≥H+n−1. Both counts are odd by Rédei, so their difference is even, giving the displayed bound. This argument is independent of the integer cut allocation.

## What the cut carrier does and does not prove

The carrier is concrete and stateful: feasible odd child values equal H plus twice the directed cut values of c, but c alone does not specify the full child tournament needed for the next induction. The class of all strongly connected integer weighted digraphs is a relaxation of tournament ear responses, not a characterization. Even with gcd-one cut increments and strong support, cut values need not form an interval. The code-8 strong order-five parent below realizes cut values {3,4,5,7,8,9,10,12,14,16}, so cut values 6,11,13,15 are absent. Thus neither connected support, positivity, integrality, parity, nor lattice gcd=1 suffices for a solid odd interval. This agrees with the stronger exhaustive order-five obstruction already in THM-4118 and Session151; it is a single finite certificate, not a disproof of global coverage or of a specially constructed state bank. A viable induction still needs admissible structural states, overlap or unit-move certificates, and closure/attainment; Session151's ternary bank is a separate sufficient target whose offset gadgets remain unproved.

## Targeted checks and reproducibility

`python check.py` uses literal child construction and independent subset path DP, finds an integral edge allocation by small backtracking, checks every cut and the growth bound for strong n=3 code 2 and n=5 codes 8 and 759, and prints matrices and images. Orientation convention for codes: ordered pairs i<j in lexicographic order; bit 1 means i→j, bit 0 means j→i. `out.txt` is the recorded output. The n=5 code-8 allocation is

    c = [[0,0,0,0,4], [3,0,0,0,0], [2,5,0,0,0],
         [2,5,5,0,0], [0,2,2,3,0]].

Every c-cut reproduces (F−9)/2, including singleton values, complement cuts and all 32 subsets. The arithmetic proof is all-order; the script is a targeted check, not an exhaustive proof of its premises.
