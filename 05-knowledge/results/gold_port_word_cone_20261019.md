# Strict Gold certificates for arbitrary words of two seven-thin cores

Date: 2026-10-19. Status: **PROVED restricted-class theorem**, using the realized port-word kernel and explicit finite integer certificates. The general Gold Partition and 1/3–2/3 conjectures remain **OPEN**.

## Result first

Every block in every allowed finite word of the two cores B8,B9 below supplies two fixed comparisons with

    t_0-t_1-t_2 >= (791/6015)t_0 > 0.                       (1)

The choices depend only on the block type:

    B9: first (5,7), second (1,2);
    B8: first (5,7), second (6,7).                          (2)

The second comparison is fixed in advance, not selected after the first answer. Arbitrary changes of core type and of the five permitted seam graphs are allowed. This upgrades the companion's one-comparison balance theorem to a strict Gold theorem on its entire class, including its bridgeless, ordinal-indecomposable, seven-thin subclass.

The whole two-comparison certificate in the first or last block uses vertices within nine positions of the respective endpoint in every reference linear extension.

A short proof using only the elementary raw seam cones already gives margin 57/625. Actual-law orbit symmetry improves it to 7/60; a one-neighbor buffer yields (1). All three certificate levels are recorded below. The margins are sharp on their stated outer cones, not asserted optimal over realizable whole words or over other strategies.

## 1. Cores, seams, and the exact boundary law

Use original labels. The core B9 has covers

    0<2<7<8, 3<4<5<8, 0<6, 1<6.

B8 is its induced subposet after deleting 1. Incoming ports are V={0,3}, outgoing ports U={6,8}.

Start with the ordinal sum of any nonempty word in B8,B9. Between consecutive blocks delete exactly one of these sets of comparisons, with the first coordinate in the earlier block:

    G0={(6,0),(8,3)},
    G1={(6,3),(8,0)},
    G2={(8,0),(8,3)},
    G3={(6,0),(8,0),(8,3)},
    G4={(6,3),(8,0),(8,3)}.                                (3)

All other relations remain. The deleted pairs join maximal to minimal ports, so retained transitivity does not restore them.

For an internal core extension w let beta(w) be its longest initial run in V and alpha(w) its longest terminal run in U. Use orders

    beta: empty,0,3,03,30;
    alpha: 6,8,68,86.                                     (4)

For constructing seam kernels, temporarily add the empty alpha profile before the four displayed profiles. Define K_g(alpha,beta) to count shuffles of those port words that retain every cross-comparison not in G_g. Each entry is a count of at most six shuffles on at most four letters.

The exact interleaving theorem from `port_word_cone_aperiodic_balance_20261019.md` gives the fibre of fixed core extensions w_1,...,w_k as

    product_i K_(g_i)(alpha(w_i),beta(w_(i+1))).             (5)

Its bijection is elementary: only a terminal outgoing run can mix with the next initial incoming run; all other cross-relations are retained, and the neighboring run shuffles are independent. Summing over both exterior words gives the weight of a central core extension as

    l_(beta(w)) r_(alpha(w)).                              (6)

The full left and right exterior laws are retained in l,r; no uniform induced-core marginal is assumed.

## 2. Four joint profile matrices

For the strategy in (2), let J_ij(beta,alpha) count core extensions with the indicated two answers and exact boundary profiles. Let N=sum_ij J_ij. The four global cell counts for that block are

    x_ij=l J_ij r.

The path deficit is l G_ij r, where

    G_ij=J_(1-i,0)+J_(1-i,1)-J_ij.                         (7)

Consequently a normalized margin d/e is equivalent to

    l(e G_ij-d N)r >= 0                                   (8)

for all four paths. This is bilinear in the two boundary vectors and linear in each separately. It therefore suffices to check every pair of generating rays of their cones.

Here are the exact joint matrices. Rows and columns follow (4); matrices occur in row-major answer order 00,01,10,11.

B9:

    J_00 = [0, 40, 10, 10]
           [0, 20, 4, 4]
           [0, 40, 14, 14]
           [0, 20, 6, 6]
           [0, 20, 6, 6]

    J_01 = [0, 0, 0, 0]
           [1, 10, 5, 5]
           [3, 5, 7, 7]
           [3, 10, 9, 9]
           [3, 10, 9, 9]

    J_10 = [0, 45, 10, 10]
           [0, 30, 6, 6]
           [0, 21, 6, 6]
           [0, 15, 4, 4]
           [0, 15, 4, 4]

    J_11 = [0, 0, 0, 0]
           [3, 30, 15, 15]
           [1, 3, 3, 3]
           [3, 15, 11, 11]
           [3, 15, 11, 11]


B8:

    J_00 = [0, 0, 0, 0]
           [0, 8, 0, 0]
           [0, 8, 0, 0]
           [0, 12, 0, 0]
           [0, 12, 0, 0]

    J_01 = [0, 0, 0, 0]
           [0, 0, 1, 1]
           [0, 0, 3, 3]
           [0, 0, 3, 3]
           [0, 0, 3, 3]

    J_10 = [0, 0, 0, 0]
           [0, 12, 0, 0]
           [0, 2, 0, 0]
           [0, 8, 0, 0]
           [0, 8, 0, 0]

    J_11 = [0, 0, 0, 0]
           [0, 6, 3, 3]
           [0, 1, 1, 1]
           [0, 4, 3, 3]
           [0, 4, 3, 3]

These matrices have an independent short counting derivation. Shuffle chains (0,2,7) and (3,4,5), then append 8. There are only twenty skeletons. For B8, insert 6 in every gap after 0. For B9, choose positions for 1 and 6 with 1 before 6, then reject those placing 6 before 0. Every extension arises uniquely. Read its two comparisons and boundary runs to obtain the displayed entries. The verifier separately enumerates the core posets and checks their aggregate cells by ideal dynamic programming.

## 3. A transparent 48-product strict certificate

The raw left cone is generated by all rows of all five K_g, and the right cone by all columns after deleting the unused empty-suffix coordinate. Their distinct generators are

    L0=(1,1,1,1,1), L1=(1,1,2,1,2),
    L2=(1,1,3,1,3), L3=(1,2,1,2,1),
    L4=(1,2,2,3,3), L5=(1,2,3,3,5),
    L6=(1,3,1,3,1), L7=(1,3,2,5,3);

    R0=(1,1,1,1), R1=(1,2,2,1), R2=(1,3,3,1),
    R3=(2,1,1,2), R4=(2,2,3,3), R5=(2,3,5,3).

The all-ones vectors include empty exteriors. These cones contain arbitrary nonnegative exterior laws by (5), even before exploiting symmetries of uniform exterior extensions.

For each left ray, minimize L(625 G_ij-57 N)R over the six right rays and all four paths. The resulting eight values are

    B9: (36527,25514,7703,59317,101440,102620,45632,120173);
    B8: (7360,7240,0,5300,19932,20834,3240,20478).

Every value is nonnegative, so (8) proves the uniform strict margin 57/625 by bilinearity. This is already an all-length Gold theorem with no orbit or buffer assumption.

The exact minimum on the raw product cone is 281/2946 for B9 and 57/625 for B8. Both occur at L2,R2, with tables respectively

    B9: [[1100,465],[794,587]],
    B8: [[240,88],[150,147]].

## 4. Two stronger realized cone refinements

### Initial/terminal orbit symmetry

For a specified set of maximal elements occupying a terminal run, all orders of that run occur equally often: permuting the terminal incomparable elements is an extension-preserving bijection. The initial minimal-port version is dual. Applied to whole exteriors, this allows sums of the doubleton kernel rows or columns instead of treating their orders independently.

The resulting primitive ray lists are

    L=(1,1,1,1,1), (1,1,2,1,2), (1,2,1,2,1), (1,2,2,3,3),
      (2,3,3,3,3), (2,3,3,4,4), (2,3,6,4,8), (2,6,3,8,4);

    R=(1,1,1,1), (1,2,2,1), (1,3,3,1),
      (2,1,1,2), (2,2,3,3), (3,6,8,4).

The corresponding rowwise minima of L(60 G_ij-7 N)R are

    B9: (2492,577,4202,7515,8051,10007,6626,10879);
    B8: (485,0,165,1342,1455,1827,144,1053).

This proves margin 7/60. The exact per-core cone minima are 279/2309 and 7/60, at the second displayed left ray and third right ray. Their tables are

    B9: [[840,350],[646,473]],
    B8: [[180,64],[120,116]].

### One neighboring core as a buffer

Let Ktilde_g be K_g with the empty-suffix row removed. If the left exterior is nonempty, its last core is some type s, and the incoming boundary vector has form

    l=x N^(s) Ktilde_g, x>=0.

Thus it lies in the cone of nonzero rows of all matrices N^(s)Ktilde_g. Dually use nonzero columns of all Ktilde_g N^(s) for the right cone. Include all-ones vectors for empty exteriors. Gcd normalization gives 36 left and 15 right rays. This is an exact outer-cone construction: arbitrary further context is absorbed into x and its right analogue, with no assumed convergence.

The 540 products per core satisfy (8) with d/e=791/6015. The minimum unnormalized scaled deficit over all products and paths is

    B9: 190561;
    B8: 0.

The verifier reconstructs every ray from the displayed core matrices and the four-letter shuffle definition, and records every integer inequality in its JSON certificate. Its exact per-core normalized minima are

    B9: 23365/173249;
    B8: 791/6015.

Both use rays l=(17,19,32,19,32), r=(5,13,13,5), with tables

    B9: [[61680,26524],[48276,36769]],
    B8: [[13260,4824],[9100,8906]].

This proves (1) for every allowed word. The refinement can alternatively be read as a direct Gold application of the independently established buffer-cone theorem in `port_word_one_buffer_20261019.md`.

## 5. Why changing the second comparison matters

The old triangular strategy first(5,7), then(4,7) is not stable under these two-port seam contexts. In the second block of B9,G0,B9, its exact global table is

    [[310832,0],[186720,118680]], total 616232.

The first increasing probability is 38854/77029>1/2. On that branch the second answer is forced by 4<5, and the Gold deficit is -5432. This is a genuine finite selector obstruction, not merely a bad outer-cone ray.

Changing only the second comparison to (1,2) gives

    [[206080,104752],[162128,143272]],

and all four deficits are positive. The cone theorem shows that this repair works uniformly through every permitted surrounding word, not just this example.

## 6. Structural and endpoint consequences

The companion proves that all these words have width three and incomparability degree at most seven. Their incomparability graphs are connected, hence the posets are ordinal-indecomposable. If the final block is B9 and every B8 has an outgoing seam other than the star G2, the incomparability graph is also bridgeless and the maximum degree is exactly seven. The choices of cores and seams may be aperiodic. These inherited structural conclusions do not rely on the new Gold calculation.

The entire two-comparison certificate is endpoint-local. In a first B9 block, all labels in its strategy {1,2,5,7} precede every later-block element, so they occur among the first nine positions. In a first B8 block, 5 and 7 precede every later element. Label 6 can be preceded by at most one later-block element, because every allowed seam deletes at most one comparison incident with 6; nonadjacent-block comparisons remain. Hence 6 still lies within the first 8+1=9 positions. In a final block, none of the labels used by either strategy is an incoming port, so all lie after the whole earlier exterior and within the last nine positions. These statements hold in every reference extension, while the Gold probabilities and budgets use the full global extension law.

## 7. Exact replay and scope

`gold_port_word_cone.py` uses only the Python standard library and the packet's existing counting helpers. It independently reconstructs the twenty-skeleton counts, enumerates both small cores, checks aggregate cells by ideal DP, derives all three ray sets from legal shuffles, and checks all 48+48+540 products for both cores.

It also compares the transfer-predicted full joint tables with unrestricted ideal-DP counts for every block in four fixed assembled words: 14 block checks, including all five seam types and up to 52 vertices. The old triangular-selector obstruction is independently recounted. Normal and `python3 -O` outputs agree. The published JSON receipt in `05-knowledge/results` contains the full integer certificate; the script prints a fresh receipt without modifying files.

The infinite conclusion follows from the exact seam fibre and bilinear inequalities. It does not claim that every seven-thin poset has these cores or seams, or that every finite poset has an endpoint-local Gold witness. No general conjecture is resolved.

Primary Gold reference: Marcin Peczarski, *Comments on the Golden Partition Conjecture*, Contributions to Discrete Mathematics 12 (2017), 106–109:
https://cdm.ucalgary.ca/article/download/62349/46785/178291


Reproduce from the repository root with `python3 05-knowledge/experiments/gold_port_word_cone.py`. Its standard-library helper scripts are colocated in the experiment directory.
