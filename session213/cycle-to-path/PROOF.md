# The exact path response of the universal cycle comb

**Status: PROVED elementary identities and a scoped no-go theorem; FINITE-EXACT family spectra. No claim of a Chan–Pak completeness proof or literature priority.**

## 1. Inheritance, target, and conventions

H(T) counts directed Hamiltonian vertex words. C(T) counts directed Hamiltonian cycles modulo cyclic rotation. Paths of one vertex are allowed. Two-path covers below are ordered and both paths are nonempty.

The inherited Session212 binary-comb theorem, `session212/cycle-bridge/PROOF.md`, realizes every positive C on strong tournaments. Its endpoint and two-vertex gate laws do not turn that C into an unrestricted H. The inherited Session210 structure atlas explains why endpoint counts, cycle counts, odd-cycle packings and total H must stay distinct. Its one-apex and fixed-bank range claims are not claimed here. Sources are pinned in `SOURCES.json` to main `25fde962779e6155a0581c1680062de45d63bb37`.

The closest proved mechanism is the ordered two-path-cover decomposition of a path after removing a marked vertex. In particular `session203/uniform-incidence-extraction/uniform/PROOF.md`, Section5, already proves the stronger positive weighted-cut identity and strict negative pair coefficients for every attachment. Its exposed-step counts Q_ij are precisely the K_ij used here. The present contribution is the explicit comb-specific evaluation of all required coefficients, not general quadraticity or submodularity. The hostile example is H=35, which exists in the ambient spectrum but will fail in this precise comb family. The corrected near miss is mistaking a binary cycle code for a binary total-path code. The least-used sidecar here is the full pairwise two-cover response to the coding apex.

Concept board: marked cycle code / ordered two-path covers / attachment cut function / induced-path monotonicity / exact family inverse. Anchor: a faithful cycle-to-path interface. Niche: an exact quadratic readout with an all-length inverse. Wildcard: which core-preserving compilers are impossible?

## 2. A general attachment identity and its lost coordinate

Let U be a nonempty tournament. Let s_i and e_i count its Hamiltonian paths starting and ending at i. Let K_ij count its ordered spanning two-path covers (P,Q) with P ending at i and Q starting at j. Thus K_ii=0. Add b, setting x_i=1 when i->b and x_i=0 when b->i.

**Theorem 1 (recovered attachment identity).** The exact path count is

    H(U+b) = sum_i s_i(1-x_i) + sum_i e_i x_i
             + sum_(i != j) K_ij x_i(1-x_j).                    (1)

Consequently the multilinear quadratic coefficient of x_i x_j is

    -(K_ij+K_ji),                                               (2)

and all second discrete derivatives are nonpositive. In particular attachment path count is a submodular set function of the in-neighborhood of b.

**Proof.** Partition a Hamiltonian path according as b is first, last, or internal. The first two cases give the endpoint terms. In the third, deleting b leaves precisely an ordered two-path cover with endpoint i and initial vertex j, and exactly the two displayed incidence conditions are required. Conversely each counted object inserts b uniquely. The partition is disjoint and exhaustive. Expanding its indicator factors proves (2), hence submodularity. QED.

This identity is included as the recovered elementary mechanism, not a new general theorem or priority claim. Session203 also proves K_ij+K_ji>=2H(U-{i,j})>=2 for distinct vertices, so strict curvature is inherited as well. It identifies the actual obstruction to a scalar cycle readout: the response contains a whole nonnegative two-cover matrix. A cycle count alone does not encode it.

## 3. Definition of the inherited comb

Fix m>=1, q=sum_(i=0)^(m-1) epsilon_i 2^i, 0<=q<2^m, tau in {0,1}, and v in {0,1}. Start with the transitive chain x0,...,x(m+1). Every earlier chain vertex beats every later one. Add a,b with

- a->x1 and xj->a for every j!=1;
- b->x0;
- x1->b exactly when tau=1;
- x(i+2)->b exactly when epsilon_i=1;
- b->a exactly when v=1 (so v=0 is Session212's forward B+).

Call the result B(m,q,tau,v). It has m+4 vertices. The inherited exact cycle identity is

    C(B)=q+tau+1-v.                                            (3)

The following theorem computes a different observable.

## 4. Exact quadratic path formula

Write f_i=[i=0] and l_i=[i=m-1], where brackets denote 0/1 indicators. Put

    U_i = 2^(m-1-i) + 2^i,
    D_i = 2^(m+1-i) 3^i + 2^(i+2) - 2 l_i,
    L_i = D_i - (1-v)(U_i-f_i+l_i) - tau(U_i+f_i-l_i),

    B0 = 2^m+1 + (1-v)(4*2^m-2)
                   + tau(2^(m+1)+2) - tau(1-v)2^m,

    W_d = 2^(m-d)(3^(d-1)+1) + 2^(d-1) + [d=1]
          for 1<=d<=m-1.

**Theorem 2.** For all the above parameters,

    H(B(m,q,tau,v))
      = B0 + sum_i L_i epsilon_i
           - sum_(i<j) W_(j-i) epsilon_i epsilon_j.             (4)

The pair interactions are strictly negative and depend only on bit-distance. Formula (4) needs O(m^2) exact integer operations, whereas it counts paths in an actual tournament on m+4 vertices. It does not assert monotonicity as the integer q increases.

### Proof by two-cover counts

Delete b and call the remaining tournament U. Vertex x0 is a universal source. On x1,...,x(m+1),a the tournament is a transitive order with just its extremal arc reversed, a->x1. The standard two-increasing-run count gives

    H(U)=2^m+1,
    e_(x(i+2))=2^i, e_(x1)=e_a=1,
    s_p=0 for p!=x0.                                         (5)

For completeness, a path in the nontrivial core either is the increasing order ending at a, or uses a->x1 and has two increasing runs. The path ending at xj has 2^(j-2) choices for j>=2, and there is one ending at x1. This proves (5) by summing the geometric series.

In (1), take the baseline b to be a universal source. For a vertex p of U, its linear coefficient at that baseline is

    A_p=e_p-s_p+sum_(j != p) K_pj.                              (6)

Here are the required cover counts. In the table p=x(i+2), r=x(j+2), i<j, and d=j-i:

    A_p = D_i,
    A_a = 4*2^m-2,
    A_(x1) = 2^(m+1)+2,
    K_(a,x1)+K_(x1,a) = 2^m,
    K_(p,a)+K_(a,p) = U_i-f_i+l_i,
    K_(p,x1)+K_(x1,p) = U_i+f_i-l_i,
    K_(p,r) = 2^(m-d)+[d=1],
    K_(r,p) = 2^(m-d)3^(d-1)+2^(d-1).                       (7)

We supply the counts rather than infer them from numerical fitting. Every chain segment in a path is increasing, and every nonempty segment immediately after a must start with x1.

For K_(p,r), put a in the first or second cover path. If it is in the first, that path has the form A,a,B ending at p, with B beginning at x1. The second chain path C begins at r. Vertex x0 is forced into A; the i vertices x2,...,x(p-1) each choose A or B; the m-1-j vertices above r each choose A or C; the vertices strictly between p and r are forced into A. This gives 2^(m-d-1). If a is in the second path and its post-a segment is nonempty, the same count results. If that segment is empty, there is one extra cover exactly when p and r are consecutive. This gives K_(p,r).

For K_(r,p), the corresponding two nonempty post-a cases again have factors 2^(m-d-1) each, but the d-1 vertices between p and r can now choose any of the three chain segments. They give 2^(m-d)3^(d-1) in total. The empty post-a case has 2^(d-1) choices, giving the final line of (7).

For the special endpoints, direct two-run allocation gives

    K_(p,a)=2^i+l_i,             K_(a,p)=2^(m-1-i)-f_i,
    K_(p,x1)=2^i-l_i,            K_(x1,p)=2^(m-1-i)+f_i,
    K_(x1,a)=0,                 K_(a,x1)=2^m.                 (8)

For example, in K_(p,a) the path beginning at a either is the singleton a, possible only when p is the last chain vertex, or continues at x1. The latter case leaves i intermediate vertices free between the two paths. In K_(a,p), all vertices below p lie in the path ending at a and vertices above p can choose either path; exclude the one case where p=x2 and the first path would enter a from x1. The other formulas follow from the same forced-start rule; (8) also states the endpoint corrections explicitly.

To check A_p without hiding a summation, count all ordered covers with their first path ending at p=x(i+2). If a lies in that first path, assign chain vertices to its pre-a and post-a segments and to the second path. The post-a segment starts at x1 and ends at p. There are

    2^(m-i)3^i - 2^i

assignments, subtracting those with an empty second path. If a lies in the second path and its post-a segment is nonempty, there are 2^(m-i)3^i assignments. If that segment is empty, there are 2^(i+2)-2l_i assignments; at the last endpoint exclude exactly the two chains entering a from x1. Adding e_p=2^i as in (6) gives D_i.

For A_a, splitting the chain into the path before a and the other nonempty path gives 2^(m+2)-1 assignments, excluding the two whose first path enters a from x1. Thus sum_j K_(a,j)=2^(m+2)-3, and adding e_a=1 gives A_a. For A_(x1), placing a in the first path gives 2^(m+1)-1 covers and placing it in the second gives exactly two. Adding e_(x1)=1 proves its displayed coefficient.

Finally set the attachment bits x_a=1-v, x_(x1)=tau, x_(x0)=0, and x_p=epsilon_i in (1). The constant, linear and pairwise coefficients are exactly (4) by (5)-(8). QED.

### A concrete failure of scalar monotonicity

At m=3,tau=0,v=0, the codes q=6 and q=7 give (C,H)=(7,101) and (8,97), respectively. The cycle code increases while the path readout decreases. The negative pair terms explain the failure; no rounding or numerical approximation is involved.

## 5. A complete finite inverse with no hidden length cutoff

**Lemma 3 (inherited insertion monotonicity, elementary proof).** If U is an induced subtournament of T, then H(U)<=H(T).

Insert each new vertex z into an old Hamiltonian path p1,...,pn. If z->p1, put z first. Otherwise take the first j with z->pj; if it exists, then p(j-1)->z->pj and insert there. If none exists, put z last. The old path remains the subsequence on its original labels, so this deterministic insertion is injective. Repeat for all new vertices. QED.

Every comb contains the U from (5), hence

    H(B(m,q,tau,v)) >= 2^m+1.                                 (9)

**Corollary 4.** To decide whether a target h>=3 is realized by any member of this comb family, at any order, it is enough to evaluate (4) for

    1<=m<=floor(log2(h-1)), all q, tau, v.                    (10)

The same bound with h replaced by R gives the entire all-length family support through R. This is an exact terminating inverse for this family, not for arbitrary tournaments. A Gray-code implementation updates one coordinate and the m first derivatives at each step; all values through R cost O(R log R) integer operations and O(R) storage if one witness per attained value is retained.

`RESULT.json` reports the exact computation through R=65535. Its universe is all four parameter choices and all bit words at every permitted m, not isomorphism classes. It also provides the complete support through35:

    {3,5,9,13,15,17,19,23,25,27,29,31,33}.

In particular this universal cycle family does not itself realize H=11 or H=35 at any order. The smallest absent odd target at least23 is35. This is not an ambient exclusion: the general tournament path-spectrum problem remains open, and35 is an inherited attained value.

## 6. Why arbitrary induced-core-retaining cycle compilers cannot work

**Theorem 5.** There is no operation F on all strong tournaments such that

1. T occurs as an induced subtournament of F(T); and
2. H(F(T)) is bounded by a finite function of C(T) alone.

In particular no such operation can satisfy H(F(T))=alpha C(T)+beta with fixed constants.

**Proof.** Let T_n, n>=3, be the transitive tournament on 0,...,n-1 with the single extremal arc reversed to (n-1)->0. Every spanning cycle uses this arc, and the remainder is the unique increasing chain, so C(T_n)=1 and T_n is strong. The two-increasing-run count used in (5) gives H(T_n)=2^(n-2)+1. Lemma3 forces H(F(T_n))>=2^(n-2)+1, unbounded at the fixed cycle count1. QED.

**Corollary 6 (the obstruction already occurs inside the marked binary family).** Fix any k>=1. The strong combs B(m,k,0,1), for m>=max(1,bit_length(k)), all have C=k, while H>=2^m+1 is unbounded. Thus a compiler retaining all such padded marked inputs as induced subtournaments cannot make its output H depend only on k either. The binary construction itself, without any arbitrary-input generalization, already has unbounded H on every positive C fiber.

This excludes only compilers retaining the stated entire input classes as induced subtournaments. It says nothing against deletion, reorientation, input-size restrictions, specially controlled marked families, or a construction starting afresh from the integer C. In particular choosing one canonical minimal-length input for each k avoids this obstruction and remains a possible route. It does not resolve Chan–Pak. The useful next obligation is to control or replace the two-cover matrix rather than transport just the cycle scalar.

## 7. Verification and current frontier

Run `python3 verify.py > RESULT.json` and `python3 -O verify.py > RESULT.optimized.json`, then compare them. The checker uses explicit exceptions, not assertions. It contains separately implemented subset-DP endpoints, literal permutation controls, and a cover-matrix computation by all ordered nonempty subsets. It checks the general attachment identity on every labeled core through order four and every attachment; it checks (4) on every comb with m=1,...,8, all bits and both apex choices; and it checks the coefficient identities, finite inverse and negative controls. The all-order proof does not depend on extrapolating any of those finite checks.

The global question remains: does every odd h>=23 occur as H(T)? The cycle comb settles the cycle statistic and gives a useful exact total-path test language, but that language is insufficient. No unbounded coverage conclusion is claimed here.
