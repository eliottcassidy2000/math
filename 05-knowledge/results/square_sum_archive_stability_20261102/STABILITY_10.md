# Sharp weighted tournament stability at order ten

Status: **PROVED with an exhaustive finite graph certificate**. Analytic
reduction audited by an exact symbolic ledger; the finite certificate was
replayed by two different exact counting kernels on the entire universe.
Same-session audit, not an independent-agent or proof-assistant audit.

## Statement and inherited notation

Let A be any real 10 by 10 matrix with zero diagonal, entries in [0,1], and
A_ij+A_ji=1. Let H be its weighted directed Hamilton-path sum,
u=det(I+A), v=det(A), q_ij=A_ij A_ji, V=sum q_ij, and
Delta=H-u-|v|. Let C_j count directed simple j-cycles once up to cyclic
rotation, retaining their direction; put tau=C_3. Then

    Delta >= (39/64)V + (8849/18432)tau
              + C4/32 + C5/16 + C6/2 + C7 + C8/2 + C9.       (1)

In particular the sharp variance constant at order ten is **39/64**.
Equality with V>0 in Delta >= (39/64)V holds exactly for ordinal sums of
five fair two-vertex blocks, up to relabeling.

This extends the inherited sharp small-order range 3..9 to **3..10**.
It says nothing new about exclusion of raw H=7 or H=21: those exclusions
already hold at every order. Nor does it identify raw H with H/Aut.

Dependencies, all at baseline commit 66cf347030dc (full pin in ENVIRONMENT.json):

- `../session137_overnight_proof_frontiers_20261102/discriminant/GLOBAL_PROOF.md`:
  odd-cycle-family expansion of H, determinant identities and margin signs.
- `../session138_stability_and_actual_contexts_20261102/stability/TRIANGLE_VARIANCE.md`:
  tau >= W(q) and W(q)>=2V^2/n-V/4.
- The same directory, `MULTISTAR_STABILITY.md`: for a nonempty r-vertex
  remainder, u+v>=2^(1-r), u-v>=r 2^(1-r), with v=(-1)^r det A.
- The same directory, `H1_STABILITY.md`: H=1 classification and its sharp
  scalar minimization. `SHARP_NINE_GRAPH_CERTIFICATE.md` supplies the closest
  method, not the new order-ten graph inequality.

## Complete cycle ledger and the distinguished-cycle multiplicity

Write M_ab... for unordered collections of pairwise vertex-disjoint directed
cycles of those lengths. A 2-cycle is the single undirected pair with weight
q_e. No cyclic collection includes isolated vertices as parts.

For a partition t into parts at least two, with support s<=10 and e even
parts, its coefficient in H is 2^|t| if all parts are odd and zero otherwise.
Its coefficient in u is (-1)^e. Its coefficient in v, when s=10, is also
(-1)^e because the ambient order is even. Thus in G_+=H-u+v a full cover
has coefficient h, and in G_-=H-u-v it has coefficient h-2(-1)^e.
For s<10 both have coefficient h-(-1)^e. Constant terms cancel.

Define

    P = V-M22+M222-M2222,
    T = tau+3M33-M223-2M2233,
    D5 = C5-M225,
    D6 = C6-M26+M36,
    D8 = C8-2M28.

For each directed 4-cycle C let R be its complementary six-set, and define

    D4+ = sum_C w(C)(u_R-v_R) + M44,
    D4- = sum_C w(C)(u_R+v_R) + M44-2M244.

The exact ledger gives, for either sign,

    G_sign >= P+T+D4sign+D5+D6+D8+M23+C7+C9.                 (2)

Every discarded coefficient is nonnegative. `verify_reduction.py` generates
ALL partitions of every s=2..10, expands both determinants, and records both
complete ledgers and their nonnegative residuals in `reduction.json`.
This is also a short directly checkable finite polynomial identity.

Crucial convention: a collection with k cycles of length four occurs k
times when one of them is distinguished. Consequently the lifted M44 and
M244 terms have multiplicity two. The correction -2M244 in D4- is necessary;
omitting it leaves a negative residual. This draft mistake was caught before
promotion and repaired, not treated as numerical error or silently ignored.

## Paying the cycle couplings

The uniform six-vertex margin is at least 1/32 in either sign. The complement
of two disjoint 4-cycles is a single pair, so M244<=M44/4. Hence

    D4+ >= C4/32,             D4- >= C4/32.

The five-vertex complement of a 5-cycle has fifteen two-edge matchings,
each of weight at most 1/16. Therefore D5>=C5/16.

For any weighted four-vertex tournament R, choose a consistent cyclic skew
sum s_T on each of its four triples. The elementary triangle identity is

    4 tau_T = 1 + (s_T^2 - sum_(e in T) K_e^2)/2,
    K_ij=A_ij-A_ji.

Each edge belongs to two triples. Summing gives

    V_R = tau_R + 1/2 - (1/8)sum_T s_T^2.

Consequently 1-V_R+tau_R>=1/2, and grouping by a distinguished 6-cycle
proves D6>=C6/2. Similarly M28<=C8/4 gives D8>=C8/2.

For two disjoint triangles the complementary four vertices have three
two-edge matchings, giving M2233<=(3/16)M33. Thus

    T >= tau+(21/8)M33-M223
      = sum_(directed triangles C) w(C)[1+(21/16)tau_R-M22_R],

where R now has seven vertices. Here every triangle pair is counted twice.
On this remainder, the inherited variance inequalities give

    M22_R <= V_R^2/2-W_R,
    tau_R >= W_R >= 2V_R^2/7-V_R/4,

and therefore

    M22_R-(21/16)tau_R
      <= -(9/56)V_R^2+(37/64)V_R
      <= 9583/18432.

The maximum is exact by completing the square. Hence

    T >= (8849/18432)tau.                                  (3)

The numbers here are safe uniform estimates, not claimed optimal triangle
reserves. No random-matrix extrapolation enters (2) or (3).

## Graph endpoint lemma

Let P21(q) be the sum of products over a two-edge wedge together with a
vertex-disjoint edge. Summing the local triangle-variance inequality with
the disjoint edge gives M23>=P21(q). We prove

    P+P21(q) >= (39/64)V,          0<=q_e<=1/4.              (4)

The difference is separately affine in all 45 variables. Its value at
any point of the box is a convex combination of its values at box vertices.
At a box vertex let G be the undirected graph with q_e=1/4 on its edges
and zero elsewhere. Write m for its edge count, m_k for its k-edge matching
count, and p21 for its wedge-plus-disjoint-edge count. Then

    256[P+P21-(39/64)V]
       = 25m-16m2+4m3-m4+4p21.                            (5)

`graph10.cpp` checks (5) on all **12,005,168** unlabelled simple graphs
of order ten, with no connectivity, degree, density, or orientation filter.
Graph isomorphism is legitimate here: (5) has no marked vertices and all
four counts are invariant under vertex relabeling. No tournament isomorphism
class or observer gauge is inferred from this ordinary graph quotient.

The complete generation uses nauty 2.9.3 `geng -q 10`. Two exact kernels
are compared for EACH graph:

1. Subset recursion for the matching polynomial at -1/4, and a separate
   perfect-matching recursion. Set
   Z(S)=4^floor(|S|/2) sum_k (-1/4)^k m_k(G[S]). For v in S,

       Z(S) = (4 if |S| even else 1)Z(S-v)
                 - sum_(w adjacent to v in S) Z(S-v-w).

   If h is the number of perfect matchings of G, then four times (5) is
   1024-Z(V(G))-h-156m+16p21. The h correction removes the otherwise
   present fifth-matching term; its sign is essential.
2. Explicit pattern masks for all 630 two-, 3150 three-, and 4725 four-edge
   matchings in K10, plus a separate complementary-degree calculation of
   p21. This computes (5) directly.

All arithmetic is exact integer arithmetic; the largest objects counted
are bounded by their corresponding K10 pattern counts. The two kernels
agree throughout. The minimum is zero, attained only by the empty graph
and the five-edge perfect matching. The receipt contains the full edge-count
histogram as an additional universe control, not a substitute for generation.

The external complete generator is an explicit computational dependency;
this is not presented as a formally verified or generator-independent proof.
The two evaluators share the graph6 parser and generator, but not the matching
or wedge-counting formulas. Weighted-matrix controls separately compare
Hamilton-path dynamic programming and rational determinants with (1).

Equation (4), (2), (3), and the other cycle estimates prove (1) for EVERY
real weighted matrix of order ten. The finite check is a complete reduction,
not empirical evidence extrapolated to real weights or larger orders.

## Sharpness and equality

Five fair ordinal pairs have

    H=1, u=(3/4)^5=243/1024, v=-(1/4)^5=-1/1024,
    V=5/4, Delta=195/256=(39/64)V.

Conversely positive-variance equality in the unreserved bound forces tau=0
by (1). Every binary completion with positive rounding probability is then
triangle-free, hence transitive. Therefore H=1. The inherited H=1
classification and strict scalar minimization leave exactly five fair
ordinal pairs. The graph zero classification is not needed for this final
equality argument, so it does not conceal a fractional-interpolation gap.

## Reproduction

From this directory, with clang++ and nauty `geng` on PATH:

```bash
clang++ -O3 -std=c++17 graph10.cpp -o /tmp/square-stability-graph10
set -o pipefail
geng -q 10 | /tmp/square-stability-graph10 dual 12005168 > /tmp/graph10.replay.json
cmp graph10-dual.json /tmp/graph10.replay.json
python3 verify_reduction.py > /tmp/reduction.replay.json
python3 -O verify_reduction.py > /tmp/reduction.optimized.replay.json
cmp reduction.json /tmp/reduction.replay.json
cmp reduction.json /tmp/reduction.optimized.replay.json
```

The matrix controls are finite corroboration; the all-weight conclusion
uses the analytic grouping and complete endpoint certificate above.
The next unresolved small-order case is eleven. No extension through
eleven, through the remaining finite gap, or to an all-order sharp constant
is claimed here.
