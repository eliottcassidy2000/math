# Square-sum Hamilton paths as exact matched-source inputs

Status: general identities proved below; finite path examples checked independently by code. No complete spectrum claim, no claim that all backward tournament edges have square-sum endpoints, and no novelty claim.

## A direct strong attachment transport

Let Z_j, j=0,1,2, be S161's five-vertex sequence: start with 1<2<3<4<5 and consecutively reverse 1–3, then 1–4. Add 6 and 7 with Z_j→7→6→Z_j. These are strong with H=33,51,69 and the intermediate graph literally shared. Vertex 7 has exactly one backward out-edge 7→6; deleting 7 leaves 6 as source, so no Hamilton path ends at 6. Attach labels 8,9 and reverse 9→7: 9+7=16. By S158's attachment lemma the three counts become 99,153,207, with the same old-edge flips. The new top has sole backward out-edge 9→7, and the lemma's endpoint obstruction propagates. Attaching 16→9, with interior 10..15 and 9+16=25, multiplies all three counts by 65. This produces 6435,9945,13455 on 16 vertices.

More generally take u_0=7, u_i=s_i²−u_{i−1} and u_i≥u_{i−1}+2. Add each interval of new vertices with sole backward edge u_i→u_{i−1}. Then the chain is (33Q,51Q,69Q), Q=∏(1+2^(u_i−u_{i−1}−1)). All graphs remain strong and the shared source remains literal.

## Profile of a single reversed interval

For d≥1 let I_d have vertices 0,…,d, with all edges forward except d→0. For r ordered path-cover slots, including empty slots,

  F_{I_d}(r)=r^(d+1)+r[(r+1)^(d−1)−(r−1)^(d−1)].

Every directed cycle uses d→0. It selects t interior vertices, t odd for an odd cycle; there are binom(d−1,t) such cycles. No two cycles are vertex-disjoint. The odd-cycle packing identity contributes 2r^(d−t) per cycle; summing gives the displayed formula. Set q_d=F(1), t_d=F(2). Thus q_1=1, q_d=1+2^(d−1) for d≥2, and t_d=2^(d+1)+2(3^(d−1)−1) for d≥1. In particular the d=1 reversal is transitive, and must not be treated as a factor 2.

## Transfer any square-sum Hamilton path

Take a Hamilton path p of the graph on 1,…,n in which distinct x,y are adjacent if x+y is a square. Put d_i=|p_{i+1}−p_i|. Let [x_i,y_i] be successive, pairwise disjoint intervals after the 5-vertex Z_j, satisfying y_i−x_i=d_i and x_i+y_i=s_i². One can always choose s_i sufficiently large and of the same parity as d_i, giving x_i=(s_i²−d_i)/2, y_i=(s_i²+d_i)/2. In the interval union reverse precisely y_i→x_i. All other edges in the block B, including gap labels, point forward. Let g count singleton gap labels in B, including those before its first interval. Then

  P=H(B)=∏q_{d_i},    R=F_B(2)=2^g ∏t_{d_i}.

Ordered slot profiles and H multiply under ordinal sum. S161 gives (H(Z_j),F_Zj(2))=(1+2j,32+16j). Applying the wrapper G(X)=X→a→b→X to X=Z_j⊕B gives strong tournaments whose exact chain is

  P+32R → 3P+48R → 5P+64R,   common jump 2P+16R.

The old flips 1–3 and then 1–4 are unchanged. Every selected edge in B has square-sum endpoints. The square graph path is translated into a sequence of interval gaps; the selected reversed edges are re-embeddings, not the original graph edges. Gap order can affect labels and g, so a single path does not prescribe one count without specifying the placement.

Adding an interval of gap d with h preceding singleton gaps updates (P,R) to (q_d P, 2^h t_d R). This is an exact nonautonomous two-state recurrence for both the common jump and absolute source. For paths with equal multisets of gaps, choose square-sum embeddings and pad the shorter B to a common length; both P and R agree. Distinct paths can thereby yield distinct labelled tournaments with the same three exact counts, though no non-isomorphism follows.

## Checked examples

A Hamilton path on [15] is 8,1,15,10,6,3,13,12,4,5,11,14,2,7,9. With least-square placements after labels 1..5, its interval endpoints are (9,16),(25,39),(58,63),(70,74),(83,86),(93,103),(112,113),(124,132),(144,145),(159,165),(179,182),(194,206),(218,223),(241,243); g=143. Here P=464858633025922867875 and R=2286430711049262759319343471689019920235263576776074977652877326066121533571465216. The final strong chain has 245 vertices.

Two distinct square-sum Hamilton paths on [25] have identical gap multisets:

A: 3,22,14,2,23,13,12,4,21,15,10,6,19,17,8,1,24,25,11,5,20,16,9,7,18.

B: 18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,3,22,14,2,23.

Least-square placements have last labels 846,623; pad B to 846. The interval-backward edge sets differ, so the resulting 848-vertex tournaments differ as labelled graphs, yet their three source counts and both jump counts agree exactly. Here g=592 and P=10918691908903704879820836527116367337011525929450009908515625; R is the explicit 2^592 times the t_d product printed by collision_check.py.

Controls: profile_check.py compares F at r=1,2,3,4 for d=1,..,7 against subset/endpoint plus slot-convolution DP, and directly counts the three nine-vertex strong attachments. graph_path_example.py verifies the [15] path, intervals and arithmetic. collision_check.py verifies both [25] paths, gap-multiset equality and distinct labelled edge sets. All checks remain active under Python -O. count16.py independently counts the three 16-vertex tournaments. These finite checks do not substitute for the general arguments.
