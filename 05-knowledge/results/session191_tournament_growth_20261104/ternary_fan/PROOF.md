# A two-port ternary fan with a two-arc nonovershooting sibling patch

**PROVED, explicit all-clock construction.** For every integer k>=1 there is a strong tournament with H=N_k=5+7·2^k and two strong tournaments with H=3N_k+2 and H=3N_k+4. They have a shared literal construction, a common old unit-start vertex, and a two-arc path between the children whose intermediate count is below both endpoints. This is one explicit infinite fan, not arbitrary-target selection or a recursively closed ternary grammar.

## 1. The response read by the exterior

For a tournament A and a physical vertex r, let F_r(A) count ordered spanning two-path covers (L,R), with L nonempty and ending at r, and R allowed empty. A one-vertex path is allowed. Let J_A(v,r) count those covers with L also starting at v. Thus sum_v J_A(v,r)=F_r(A).

For k>=1 adjoin a transitive block Z=(z_1,...,z_k) and w, with

    w -> Z -> A,     r -> w -> A-r.

Call the tournament W_k(A,r). The already established Session190 identity is

    H(W_k(A,r)) = H(A) + 2^k F_r(A).                  (1)

Its full start vector has the following particularly useful form:

    s_v(W_k) = J_A(v,r)          (v in A),
    s_w(W_k) = H(A),
    s_z_i(W_k) = 2^(k-i) F_r(A) (1<=i<=k).         (2)

**Direct proof.** A path starting w must next traverse all Z and then a Hamiltonian path of A, giving H(A). A path starting in A must reach w through r. Before w it is an A-path from v to r; after w it traverses all Z followed by the remaining A-path. This is exactly J_A(v,r). A path starting z_i has a nonempty Z-block followed by an A-path ending r, then w, a possibly empty second Z-block, and a possibly empty second A-path. Every z_j with j<i is forced into the second Z-block; each z_j with j>i independently chooses either block. The A-parts are precisely a cover counted by F_r(A), giving the last equation. Summing (2) proves (1).

If A is strong, W_k is strong: A reaches r and hence w; w reaches Z and A; and Z reaches A. Also w has sole in-neighbor r and W_k-w=Z->A, so s_r(W_k-w)=0. Reversing every arc therefore gives a strong tournament with sole-exit mark w and zero deleted-core digit. Without taking the converse, w is the corresponding sole-entry mark. This is the exact zero-response property, not digit-one readiness.

The complementary two-bin identity E_w(W_k,2)=H(W_k) follows too: for a vertex w with sole in-neighbor r, splitting covers according to whether their first path is the singleton w gives E_w-H=s_r(W_k-w), which is zero here. The k>=1 hypothesis is essential.

## 2. Five small cores on physical names

Keep r=1 and v=3 throughout. Outgoing adjacency rows are bitmasks on names 0,1,... .

| Core | Rows | H | F_r | J(v,r) | s_r |
|---|---|---:|---:|---:|---:|
| P | (10,4,9,2) | 5 | 7 | 1 | 1 |
| Q | (26,20,25,2,8) | 9 | 11 | 1 | 1 |
| A | (58,52,57,34,8,16) | 17 | 21 | 1 | 3 |
| C | (58,52,57,2,8,24) | 15 | 19 | 1 | 1 |
| B | (26,52,57,2,8,25) | 19 | 21 | 1 | 1 |

All five cores are strong. These small counts are checked independently by literal permutation/cut enumeration and predecessor-subset recurrences in `verify_fan.py`; the verifier generates `fan_receipt.json` with the full start, end, and J vectors.

The growth and patch are literal:

1. From P add x=4 with sole old out-neighbor v=3, obtaining Q.
2. Add y=5 with sole old out-neighbor x=4, obtaining A.
3. Reverse 3->5 to 5->3, obtaining C.
4. Reverse 0->5 to 5->0, obtaining B.

The old P-incidences, root r, port v, and every exterior name and arc are unchanged during growth. The two patch stages change exactly the displayed old-core arcs. Newly added core vertices are dominated by Z and w, as required by the wrapper. At every stage the physical old vertex v remains a unit-start vertex of the lifted tournament, by (2).

## 3. Exact ternary children and the nonovershooting route

Put x=2^k, k>=1. The five lifted counts, in the displayed construction order, are

    5+7x,  9+11x,  17+21x,  15+19x,  19+21x.

Consequently for N_k=5+7x the children have exactly

    H(W_k(A,r)) = 3N_k+2,
    H(W_k(B,r)) = 3N_k+4.

The two child cores use the same six names and root. Their two-arc route passes through W_k(C,r), and

    (17+21x)-(15+19x)=2+2x>0.

Thus the sole intermediate is below the lower endpoint for every permitted clock. The parent has k+5 vertices; the common growth intermediate has k+6; both children and their patch intermediate have k+7. The first four parent/children triples are

    (19;59,61), (33;101,103), (61;185,187), (117;353,355).

The output is strong at every stage and has the same named unit-start port v=3, as well as the same outer sole-entry zero-response mark w. These are exact incidence statements, not scalar substitutes for marked tournament equality.

## 4. Why this does not contradict the Session190 obstruction

Session190 required the *same inner vertex* to satisfy both s_r=1 and equal F_r at the patch endpoints. Here the lower core has s_r=3. The preserved unit-start requirement is instead carried by the second physical vertex v, via J(v,r)=1. In particular this construction does not refute Session190's seven-name minimum for its one-port criterion or its four-edit minimum for its specified endpoints.

The useful interface here is (H,F_r,J(v,r)), not the complete rooted cover response. The two child cores differ in other observables. Neither the parent family nor its children have been shown to admit unlimited repetition, and no statement is made that every attainable H has such two children. Second-generation coefficient constraints require a separate calculation.

## 5. Reproduction and credit

Run, without third-party packages:

    python3 verify_fan.py
    python3 -O verify_fan.py

Both runs produce byte-identical output and the same receipt. They independently verify all five core profiles, all literal births/flips, and 30 physical wrappers at k=1,...,6, including their full start vectors. The outer two-bin zero-response identity is independently computed for k=1,2. The all-clock statements follow from the proofs above, not extrapolation of those samples. Explicit `require` checks remain active under `-O`.

The inherited exterior is Session190's reviewed [four-arc flat-response patch](https://github.com/eliottcassidy2000/math/blob/3cb7969b5f4b25d1a4b4a5e84d33f1bb3b659081/05-knowledge/results/session190_constructive_links_20261104/spectrum/spectrum/PROOF.md). Its joint-entry viewpoint also occurs in Session189. No priority claim against the wider literature or incoming repository work is made. The new content in this packet is the explicit two-port physical fan, two-flip nonovershooting patch, and its exact marked-state distinctions.
