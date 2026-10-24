# A marked regularization bridge, and exactly what its unmarked quotient loses

24 October 2026. **CITED / KNOWN:** the base Type-I regularization is Lemma 4.3
of Hatem–Mohamed–Müller–Siniora, [*Cofinality of Regular Tournaments*,
arXiv:2605.04846v1](https://arxiv.org/html/2605.04846v1#S4).
**PROVED:** the structural statements below. **FINITE-EXACT:** the explicitly
bounded checks and first-collision thresholds. No literature-priority claim
is made for the extracted statements, and no new formula for A000568 is claimed.

## 1. Inheritance and the connection contract

The inherited mechanism is the known two-copy Type-I regular completion.
The nearby but different construction is repository THM-4557,
`01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md`:
that skew-Hadamard doubling has a **converse** second copy, a different cross
rule, and a doubly regular input hypothesis for its automorphism theorem.
It is not the map in this packet. Its neighborhood-signature proof is a
related method, not an identification of the constructions.

Anchor: retain the entire source class while entering the regular stratum.
Niche: lift individual source arc changes to degree-preserving circulations.
Wildcard: forget the root and measure the resulting class collisions.
Concept board: rooted isomorphism; actual source arcs; four-cycle circulations;
Hamiltonian paths; native complement; converse merging; root recognition.

The map sends an arbitrary tournament T to a regular tournament D(T) with a
root z. Its root out-neighborhood recovers T. A source Hamiltonian path becomes
a specially matched alternating Hamiltonian cycle. The exact source arc
marks become four-edge blocks; a source nonpath complement reverses those
blocks, not the whole target noncycle edge set. The root is essential in
general: source-to-sink recentering produces explicit quotient collisions.
The root is intrinsically recoverable on the no-source/no-sink subdomain.

## 2. Definition and the known regularity calculation

Let T have vertices 0,...,n−1, n≥1, and adjacency matrix A, where A_ij=1
means i→j. Let U={U_i}, V={V_i}, and add z. Define D(T) by

- both induced copies U and V have orientation T;
- U_i→V_j iff i=j or j→i in T;
- z→U_i and V_i→z for every i.

This is exactly the cited Lemma 4.3, after naming its two copies U,V and its
new vertex z. In the order U,V,z, its adjacency matrix is

    [ A       I+A^t    0 ]
    [ A^t     A        1 ] .                              (1)
    [ 1^t     0^t      0 ]

Since A+A^t=J−I, the off-diagonal blocks sum with their transpose to J.
Writing d_i for the outdegree of i in T, U_i has d_i+(n−d_i)=n outgoing
arcs, V_i has (n−1−d_i)+d_i+1=n, and z has n. Thus D(T) is regular of
order 2n+1, with both z-neighborhoods isomorphic to T. These facts are
**known**, reproduced here to fix the orientation conventions.

## 3. Rooted faithfulness and the automorphism stabilizer

**MR1.** Every isomorphism f:T→S has exactly one root-preserving extension
D(T)→D(S), namely

    U_i→U_f(i), V_i→V_f(i), z→z.                         (2)

Every root-preserving isomorphism is of this form. Consequently

    Aut(D(T))_z ≅ Aut(T),
    (D(T),z) ≅ (D(S),z) iff T ≅ S.                     (3)

**Proof.** A root-preserving map takes its out-neighborhood U onto U and its
in-neighborhood V onto V; restriction to U is an isomorphism f:T→S. The
signature of V_i on U is

    {j : V_i→U_j} = N_T^−(i).                           (4)

Distinct vertices of a tournament have distinct open in-neighborhoods:
if i→j, then i belongs to N^−(j) but not N^−(i). Therefore (4) forces
V_i to map to V_f(i). Conversely (2) respects every defining arc. This proves
both existence and uniqueness. In particular, the U-to-V matching is
recoverable from the rooted tournament; it is not an additional arbitrary
matching choice. ∎

**MR2 (converse).** Interchanging U and V and fixing z gives

    D(T*) ≅ D(T)*.                                      (5)

Check the internal, matching, root, and off-diagonal arcs directly in (1).
Thus rooted faithfulness also holds modulo converse. All uses of converse
below mean reversal of **every** source arc, not native path complement.

## 4. Source-edge changes are edge-disjoint directed four-cycles

For an unordered pair {i,j}, let

    B_ij = {U_iU_j, V_iV_j, U_iV_j, U_jV_i}             (6)

be a set of four underlying target edges, retaining their actual orientations.
If i→j in T, these edges form the directed cycle

    U_i→U_j→V_i→V_j→U_i.                                (7)

**MR3.** Reversing source arc i→j reverses exactly the four arcs of B_ij in
D(T). Distinct source pairs have disjoint edge blocks. For every set E of
source edges,

    D(reverse_E(T)) = reverse_(union_{ij in E} B_ij)(D(T)).   (8)

Every such simultaneous reversal preserves all target degrees.

**Proof.** Each source adjacency entry governs its two same-copy entries and
its two cross-copy entries in (1), precisely the four entries in (6).
Their orientation is (7), and a directed cycle has one incoming and one
outgoing changed edge at every participating vertex. No target edge belongs
to two source-pair blocks. Other reversals therefore do not alter the four
arcs of this block, and each block balances independently. ∎

The remaining 3n target edges are the n fixed triangles z→U_i→V_i→z.
Hence the labelled source orientation cube embeds into regular tournaments
as a four-edge-block cube. Labelled Hamming distances multiply by four.
This statement retains the block labels. It is not an unmarked isometric
embedding of tournament-class edit graphs.

## 5. Paths, distinguished cycles, and the native complement boundary

For a source Hamiltonian path P=(i_1,...,i_n), set

    L(P)=(z,U_i_n,V_i_n,U_i_(n−1),V_i_(n−1),...,U_i_1,V_i_1).   (9)

A **distinguished cycle** means a directed Hamiltonian cycle written from z
which next visits U, then traverses each matched U_i,V_i consecutively,
alternating in this way until returning to z. This excludes other target
Hamiltonian cycles.

**MR4.** Formula (9) is a bijection between source Hamiltonian paths and
distinguished cycles. It is equivariant for Aut(T) and Aut(D(T))_z.

**Proof.** The root and matching arcs always exist. Between successive
matched pairs we need V_i_k→U_i_(k−1), equivalent by the cross rule to
 i_(k−1)→i_k. Thus (9) is directed exactly when P is a source path. Reading
the U indices in reverse order recovers P uniquely. Formula (2) commutes
with this construction. ∎

In particular, H(T) is the number of these distinguished cycles, **not**
the total number of target Hamiltonian cycles. The action on actual paths
is free, since an automorphism fixing a whole ordered Hamiltonian path fixes
every vertex. Thus H(T)/|Aut(T)| counts their root-preserving orbit classes.
The usual path-orbit mass is inherited, not a new enumerative formula.

Let C_P reverse all source edges not used by P, keeping its path arcs.
Let C_L^full reverse every target edge outside the cycle L(P).

**MR5.** The source operation lifts exactly to reversal of the B_ij for
nonpath pairs {i,j}. It keeps L(P), is involutive, and gives D(C_P(T)).
For n≥2 it is never C_L^full on the labelled target:

    lifted source changed edges = 4[binom(n,2)−(n−1)]
                                = 2(n−1)(n−2),
    full target changed edges   = binom(2n+1,2)−(2n+1)
                                = (2n+1)(n−1).           (10)

The operations differ on exactly 5(n−1) target edges. At n=1 both are the
identity. Indeed no nonpath block meets L(P), so the first support is a
subset of the second, and subtraction gives the exact discrepancy.

This is an all-size operation bridge **with actual arc and block marks**.
It does not identify the source's native C-lines with native full-target
C-lines, does not transport port clocks, and does not recover an unmarked
metagraph after throwing away those marks.

## 6. A canonical root on a natural large subdomain

Write c(x,y)=|N^+(x)∩N^+(y)| in D(T). Directly from the arc rules,

    c(U_i,V_i)=0,
    c(z,U_i)=d_i,
    c(z,V_i)=n−1−d_i.                                  (11)

**MR6.** If T has neither a source nor a sink, z is the unique vertex x
such that c(x,y)>0 for every y≠x. Consequently every automorphism fixes z,

    Aut(D(T)) ≅ Aut(T),                                (12)

and D is faithful on unmarked isomorphism classes of tournaments with
neither a source nor a sink. It is also faithful modulo converse on that
domain.

**Proof.** The root satisfies the condition exactly when 0<d_i<n−1 for
every i. Every other vertex has its matched partner with c=0. An isomorphism
preserves this intrinsic condition, hence preserves z, and MR1 applies.
For converse, apply MR2 and note that the no-source/no-sink condition is
converse-invariant. ∎

When T has a source or sink, **no** vertex satisfies this particular
condition. This does not imply that D(T) has an automorphism moving z, or
that every such source has a class collision. For example 1+C3 has root
stabilizer and full automorphism group both of order three. The theorem's
hypothesis is sufficient, not a characterization of all root-faithful cases.

## 7. Source-to-sink recentering and explicit quotient obstructions

For tournaments A,B, let A+B denote their order join: every A vertex beats
every B vertex. Let s be the source of T=1+Q, and t the sink of S=Q+1.
Use the same labels q for the vertices of Q in both source tournaments.

**MR7.** There is an isomorphism D(1+Q)→D(Q+1) given by

    old U_s → new z,
    old V_s → new U_t,
    old z   → new V_t,
    old U_q → new U_q, old V_q → new V_q.              (13)

**Proof.** At the old U_s, the out-neighborhood is all U_q followed by V_s,
which induces Q+1. Its in-neighborhood is all V_q followed by z, also Q+1.
Under the asserted matching U_q↔V_q and V_s↔z, all diagonal arcs point from
the new out-copy to the new in-copy. The unchanged Q-to-Q blocks already
satisfy the defining rule. The remaining arcs are

    z→U_q, V_s→V_q, U_q→V_s, V_q→z,

which become exactly the four types involving the new sink t. The root
arcs hold by the two neighborhood descriptions. This covers every pair. ∎

**MR8 (ordinary loss).** Take Q=C3. Then 1+C3 and C3+1 are not isomorphic
(their source/sink structures differ), but their regularizations of order
nine are isomorphic by (13). The two source tournaments are converse, so
this example alone does **not** disprove faithfulness on converse-merged
classes.

**MR9 (merged loss, smaller witness).** Take Q=1+C3. The order-five sources

    A=1+1+C3,             B=1+C3+1                       (14)

are neither isomorphic nor converse-isomorphic, while D(A)≅D(B) has order
11. Their strongly connected component sizes in dominance order are
(1,1,3) and (1,3,1). Converse reverses this sequence; B's sequence remains
(1,3,1). An isomorphism must preserve the ordered SCC decomposition, so
neither equality is possible. Formula (13) supplies the target isomorphism.

For a second, larger witness with no singleton inside Q, take Q=C3+S4,
where S4 has arcs 0→1→2→3→0, 0→2 and 1→3. This S4 is strongly connected.
Then 1+Q and Q+1 have SCC sequences (1,3,4) and (3,4,1); the converse of
the latter has (1,4,3). They are different ordinary and merged classes,
yet their order-17 regularizations are isomorphic by (13).

The result file includes complete source and target adjacency matrices and
explicit old-to-new permutations for all three witnesses. No numerical
invariant alone is being used as evidence of target isomorphism.

## 8. Exhaustive checks, threshold scope, and hostile controls

`python3 check.py --output normal.json` and its `python3 -O` counterpart
exhaust **all 1,099 labelled tournaments of source orders 1 through 5**.
Every source is independently constructed by matrix blocks and by the
published edge recipe. Target automorphisms and isomorphisms are found by a
generic exact adjacency backtracker with no knowledge of the construction.
Source automorphisms use direct permutation enumeration. Every source
Hamiltonian path is found both by permutation testing and subset DP; every
distinguished cycle is separately scanned over all matched-pair orders.

The finite class counts are:

| source order | ordinary sources | unmarked D images | merged sources | merged D images |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 |
| 2 | 1 | 1 | 1 | 1 |
| 3 | 2 | 2 | 2 | 2 |
| 4 | 4 | 3 | 3 | 3 |
| 5 | 12 | 9 | 10 | 9 |

Here “merged” quotients by isomorphism and full converse. The finite-exact
absence of a smaller collision, together with MR8/MR9, establishes the
first ordinary collision at source order four and the first merged collision
at source order five. These are thresholds for **this map**, not an A000568
threshold or a count of all regular tournaments.

The checks include every single source-edge flip (10,650 total), every source
Hamiltonian path (7,887 total), all root-stabilizer claims and no-extreme full
automorphism claims, all small source class pairs, all 1,100 recenterings for
Q of orders zero through five, and the explicit order-17 witness. Sixty
fixed-seed examples of source orders 6 through 10 are additional sample
controls, not exhaustive censuses.

Mutation controls reject: a converse second copy; an unreversed cross rule;
reversed matching diagonals; a forward rather than reversed path; a lone
target-edge flip; replacing lifted source C by full target C; omitting the
root-recognition hypothesis; claiming unrestricted full automorphism
equality; and counting all target Hamiltonian cycles as distinguished ones.
For transitive T of order two, D(T) has two target Hamiltonian cycles rooted
at z but only one distinguished cycle, and |Aut(D(T))|=5 while |Aut(T)|=1.

The implementation uses explicit exceptions, never Python assertions, so all
validity checks survive optimized mode. Normal and optimized JSON outputs
are byte-identical. `verify_manifest.py` checks SHA-256 values of the portable
packet. No third-party Python packages, repository checkout, network access,
canonical-labelling program, or hidden atlas is required to rerun it.

## 9. An elementary class-count corollary

Let A_n count ordinary tournament isomorphism classes (A000568), and let S_n
count ordinary classes that are isomorphic to their converse. Here S means
**self-converse**, not strongly connected. Set A_0=S_0=1 for the empty
tournament. For n≥2, let N_n count ordinary classes with neither source nor
sink, and let J_n count self-converse classes with neither source nor sink.
Then

    N_n = A_n − 2 A_(n−1) + A_(n−2),
    J_n = S_n − S_(n−2).                               (15)

**Proof.** A tournament has at most one source and at most one sink. Removing
the unique source bijects source-bearing classes with all (n−1)-vertex
classes; the same holds for a sink. Having both bijects with (n−2)-vertex
classes by removing those two distinct vertices. Inclusion-exclusion proves
the first formula. A self-converse tournament with a source also has a sink,
and vice versa. Removing both leaves a self-converse tournament, and
adjoining a source and a sink to any self-converse tournament reverses this
operation on ordinary classes. This proves the second formula. ∎

Write R_m for ordinary regular tournament classes, S^reg_m for their
self-converse classes, and M^reg_m for their converse-merged classes. MR6
and its converse-compatible version give

    R_(2n+1) ≥ N_n,
    S^reg_(2n+1) ≥ J_n,
    M^reg_(2n+1) ≥ (N_n+J_n)/2.                         (16)

The last numerator counts converse orbits by fixed points. These inequalities
are elementary comparisons with regular tournaments at a **larger order**.
They give a weak construction-based lower bound, with no claim of asymptotic
sharpness, a fast evaluation of A_n, or new enumerative priority. In particular,
(15) takes the source class counts as inputs; it does not compute A000568
from regularity. The checker verifies both formulas and the merged orbit
count through source order five against the independently enumerated classes.

## 10. What this settles and what stays open

The faithful object is a rooted regular tournament, with a recoverable paired
copy and, for operation transport, actual source-edge blocks. On the
no-source/no-sink domain the root becomes intrinsic. On unrestricted sources,
source-to-sink recentering loses even converse-merged class information.
The fixed-path operation still has an exact regular circulation lift, but
it is a proper sub-operation of target noncycle complement. These are the
precise surviving bridges; no port-clock theorem, general count formula,
or full metagraph reconstruction follows.
