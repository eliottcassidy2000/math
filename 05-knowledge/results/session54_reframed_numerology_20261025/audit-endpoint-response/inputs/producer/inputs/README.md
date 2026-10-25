# Full Paley-11 minimum-blocker classification: exactly121

Status: PROVED with small exact rank certificates; independent audit requested. P11 has beta_HP=10 and exactly121 minimum blockers:55 two-source deletions,55 two-sink deletions, and11 isolated-vertex stars. The same argument recovers all63 P7 minima, including its two additional seven-element Hall families. No enumeration of binomial(55,10) deletion sets is used.

This separate packet supersedes the historical uncertainty about the P11 total in the earlier frozen packets. Those packets are not edited. Both the1023 count conjecture and the specified F1024× extrapolation are now false for the full P11 minimum-blocker family.

## 1. Arc-orbit averaging forces constant cycle intersection

Let T be an arc-transitive tournament on n=2r+1 vertices with at least one directed Hamiltonian cycle. Let B destroy every Hamiltonian path. Every Hamiltonian cycle C must meet B in at least two arcs: zero or one deletion from C leaves a spanning directed path.

For any such C, average |B intersect gC| over the automorphism group. Arc transitivity makes every arc equally frequent, so the mean is

    |B|*n / |E(T)| = |B|/r.

It is at least2. Hence |B|>=2r=n-1. Isolating any vertex gives the matching upper bound. Thus beta_HP=n-1 under these hypotheses.

If |B|=2r, the average is exactly2 and every summand is at least2. Therefore |B intersect gC|=2 for every g. Since the argument starts with ANY Hamiltonian cycle C, every Hamiltonian cycle meets every minimum blocker in exactly two arcs. This step does not assume all cycles lie in a single orbit.

For P7 and P11, the explicitly checked affine group x->ax+b (a a nonzero square) is arc-transitive, and the rank certificate supplies a Hamiltonian cycle. Thus the argument proves their blocking numbers independently of the earlier cycle-packing supplement.

## 2. A small exact rank certificate forces row-column potentials

Let m=|E(T)|. View Hamiltonian-cycle incidence vectors in Q^m. Their differences are orthogonal to every vector

    b_(u,v)=a_u-d_v,

because each Hamiltonian cycle uses each vertex once as a tail and once as a head. The bipartite support graph, with tail vertices and head vertices and one edge per tournament arc, is connected in both P7 and P11; an explicit spanning tree is recorded. Therefore the row-column potential space has dimension2n-1, and the cycle-difference rank is at most m-2n+1.

The supplied small certificates attain this bound:

- P7: one base cycle plus8 further cycles; an8 by8 minor has determinant1.
- P11: one base cycle plus34 further cycles; a34 by34 minor has determinant-2.

All cycle vertices/arcs and both determinants are checked by exact integer arithmetic (Bareiss elimination). Consequently the orthogonal complement of the certified cycle differences is EXACTLY the potential space. By section1, every minimum blocker incidence vector lies in that space.

This is a rational linear-span conclusion. In particular, the P11 determinant-2 is NOT silently treated as a unimodular or saturated integral cycle-lattice certificate.

## 3. Integrality comes from the support graph, not the cycle minor

Suppose an integer arc vector has a rational representation b_uv=a_u-d_v. Fix one potential at0 by adding the same constant to every a and d. Along a support-tree edge, the equation a_u-d_v=b_uv propagates an integer value from one endpoint to the other. Connectivity reaches all2n potentials. The pre-existing rational solution guarantees consistency, so all potentials can be chosen integral.

For a binary blocker vector, every allowed arc satisfies

    0 <= a_u-d_v <= 1.

Shift all potentials by their minimum, making their common global minimum0; write M for their maximum. For each k=1,...,M define

    S_k={u:a_u>=k},   T_k={v:d_v>=k}.

The lower inequality implies N^-(T_k) is contained in S_k. The cut D_k consisting of arcs from S_k to V\T_k has size

    |D_k|=r*(|S_k|-|T_k|).

Moreover

    b_uv = sum_k 1_{(u,v) in D_k}.

The upper inequality says these cuts are arc-disjoint. Their deficiencies delta_k=|S_k|-|T_k| are nonnegative integers.

No active level has delta_k=0. If it did, the cut D_k would be empty, while compatibility already forbids edges from the complementary tail set into T_k. The bipartite support graph would split into two components unless S_k,T_k were simultaneously full or simultaneously empty. Neither occurs at k between1 andM with global minimum0 and maximumM. Thus every active delta_k is positive.

Since the blocker has size2r,

    sum_k delta_k=2.

There are only two possibilities: one deficiency-two cut, or two nested deficiency-one cuts. This is a complete decomposition, not a heuristic binary relaxation.

## 4. Doubly regular expansion classifies deficiency-one cuts

For a doubly regular tournament, let r=2lambda+1, n=2r+1. For a set T of t vertices, put U=N^-(T), the union of their in-neighborhoods. Counting indegrees and common in-neighbors, followed by Cauchy-Schwarz, gives

    |U| >= r^2*t/(lambda*t+r-lambda).

A compatible deficiency-one cut requires |U|<=t+1. For1<=t<=2r-1,

    r^2*t-(t+1)*(lambda*t+r-lambda)
    =r*(r-2)+lambda*(t-1)*(2r-1-t)>0

when r>=3. Hence the only deficiency-one cuts are t=0, a full outgoing row at one vertex, and t=n-1, a full incoming column at one vertex.

The level sets are nested. Two row cuts would have identical singleton tail sets and overlap; two column cuts would similarly coincide and overlap. Both contradict binary arc values. The only two-level possibility is a full incoming column at w at the lower level and a full outgoing row at v at the upper level.

They are disjoint iff v=w or the tournament arc is w->v. If v=w, their union is the star at v. If w->v, their union is a mixed source/sink deletion. Any Hamiltonian path starting at w and ending at v avoids all deleted arcs: it uses neither an incoming arc at its first vertex nor an outgoing arc at its last.

The supplied forward endpoint witnesses are

- P7, from0 to1: 0,2,3,4,5,6,1.
- P11, from0 to1: 0,3,4,2,5,6,7,10,8,9,1.

Arc transitivity transports these to every ordered pair w->v. Thus every non-star mixed deletion retains an HP and is NOT a blocker. Reverse-orientation endpoint witnesses are also supplied and checked, although the mixed-cut exclusion needs only the forward orbit.

## 5. Deficiency-two cuts and the final counts

A one-level potential is exactly a compatible deficiency-two Hall cut. The analogous expansion comparison is

    r^2*t-(t+2)*(lambda*t+r-lambda)
    =r*(r-3)+lambda*(t-1)*(2r-2-t).

For P11, r=5, the right side is strictly positive for1<=t<=8. Only t=0 and9 remain: two outgoing rows (two sinks) or two incoming columns (two sources). There are55 of each. Together with11 stars, all of which certainly block every HP, this gives exactly121 minimum blockers.

For P7, r=3, the right side vanishes at t=1 and4. At t=1, compatibility forces S=N^-(c), T={c}, yielding the seven common-out-neighbor blockers. At t=4, writing S=V\{c}, compatibility forces T={c} union N^-(c), yielding the seven common-in-neighbor blockers. The endpoints t=0,5 give21 two-sink and21 two-source blockers. Adding seven stars gives63.

The categories are distinct, either from their residual degree/neighborhood profiles or the explicit certificate check. Every Hall cut and every star has2r arcs and blocks every HP. Sections1-4 show there can be no further minimum blockers.

## 6. Exact finite certificate and hostiles

verify.py does NOT enumerate all P11 deletion sets or all Hamiltonian cycles. It checks the9/35 listed cycles, the exact full-rank minors, support connectivity, regularity and double regularity, arc transitivity, endpoint paths, and an explicit difference-cycle decomposition as an additional independent beta certificate.

It independently enumerates compatible deficiency-one/two cuts through in-neighborhood unions. Combined with the rigorously exhaustive two-level alternative, the complete binary potential candidate sets have84 elements at P7 and176 at P11. Endpoint-subset HP-DP leaves exactly63 and121 blockers; all21 and55 rejected mixed candidates receive explicit surviving paths. All candidate words meet every listed Hamiltonian cycle exactly twice, so those mixed examples are a concrete hostile against the FALSE implication “constant HC intersection2 implies HP blocking.”

The actual classification proof excludes those mixed words uniformly using the endpoint-orbit witness; the DP is an independent finite check, not a hidden assumption about unenumerated deletion sets.

## Recovery, portability and scope

Targeted searches of current canon/results/hypotheses for Hamiltonian row-column potentials, blocker potentials, Hamiltonian affine hulls, and cycle-difference ranks found no prior version of this classification route. Earlier P7 counts and P11 Hall counts are recovered, not claimed new. No literature-priority claim is made.

The general criterion proved here is conditional on explicit hypotheses: an arc-transitive Hamiltonian tournament with connected bipartite support and maximal Hamilton-cycle difference rank has minimum blockers in the binary row-column potential space; the doubly regular expansion and forward endpoint witness then give the stated Hall/star classification. No maximal-rank claim for all tournaments or all Paley orders is made.

Run replay.py --self-test from any directory. It verifies the exact manifest, normal/-O output identity, and file-add/content-mutation rejection. All rank and endpoint input data are included; only Python's standard library is required. Earlier frozen packets remain unchanged.
