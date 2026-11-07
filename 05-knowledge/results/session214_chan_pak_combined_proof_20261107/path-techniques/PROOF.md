# Fixed endpoint defect has a finite dyadic classification at a sole exit

Status: PROVED elementary statements, pending independent proof audit. Computational controls are separate. Chan–Pak completeness remains OPEN. No literature-priority claim.

## 1. Scope and inheritance

Let H(T) count directed Hamiltonian vertex words, with H(empty)=1. Let e_a(T) count those ending at a. A sole exit means N+(a)={r}. We study the positive odd defect

    d(T,a)=H(T)-2e_a(T)

at such a mark, assuming e_a(T)>0. This is a restriction on the *output mark*, not on arbitrary tournaments with a specified total H.

Session213's `cycle-arithmetic/EQUALITY.md` classifies d=1 at every positive non-sink endpoint: the near-transitive extremal-arc-reversal family, with transitive prefix, and e_a a power of two. Its universal classification is inherited. Section 3 below recovers its sole-exit subcase and proceeds to arbitrary fixed d. Session203's `uniform-incidence-extraction/uniform/PROOF.md`, Section 5, supplies the strict pair capacity inequality; its short block-insertion proof is recalled below. SCC multiplicativity and the elementary lower bound H(C)>=|C| for a nontrivial strong tournament are inherited. The latter follows by rotating one Hamiltonian cycle (Camion). The weighted one-apex word count is inherited; here it also follows directly from two path covers.

The closest mechanism is the endpoint injection, the hostile is mistaking a single marked endpoint count for unrestricted H, and the additional coordinate is an ordered spanning two-path-cover capacity. Anchor: exact odd total-path counts. Niche: stability of endpoint equality. Wildcard: compress fixed-defect tournaments rather than search bounded orders without a bound. Live concepts are endpoint deficit / sole exits / ordered covers / strong-component order / dyadic padding / finite inverse.

## 2. Exact defect and endpoint identities

Put U=T-{a,r} and O=N+_U(r). For distinct x,y in U, let K_xy(U) count ordered pairs (P,Q) of nonempty disjoint directed paths covering U, with P ending at x and Q starting at y. Put

    K(U,O)=sum_(x,y in O, x!=y) K_xy(U),
    L(U,O)=#{(P,Q): P possibly empty, Q nonempty,
                    V(P) disjoint V(Q), V(P) union V(Q)=U,
                    Q starts in O}.

**Theorem 1.** For any sole-exit mark a->r, with U nonempty,

    e_a(T)=L(U,O)-K(U,O),
    H(T)-2e_a(T)=H(U)+K(U,O).                         (1)

In particular

    d >= H(U)+2 sum_{{x,y} subset O} H(U-{x,y})
      >= H(U)+|O|(|O|-1).                            (2)

Proof. A path ending at a must use r before its end; r must be followed by a nonempty path Q starting in O. It has unique form P,r,Q,a. Every vertex of U beats a. Its only possibly invalid junction is from the end of P to r: it fails exactly when P ends in O. Hence e_a=L-K.

For the second identity use the inherited injection from paths ending at a: move a from the end to immediately before r. Its image consists of paths ending in U. Such a path in U has form P,a,r,Q with Q nonempty. It lies in the image if and only if P is empty or its last vertex beats r: precisely the condition that removing a and appending it gives a path. The complement of this image is therefore counted by K(U,O). Paths ending at r have forced final segment a,r and are in bijection with paths of U, so there are H(U) of them. Partitioning paths by their last vertex gives H=e_a+(e_a+K)+H(U), proving (1).

For completeness the inherited inequality K_xy+K_yx>=2H(U-{x,y}) has a short injection proof. In each Hamiltonian word z_1,...,z_k of U-{x,y}, insert either adjacent block (x,y) or (y,x), treating its internal arc as exposed. Cutting that exposed step gives the required two-cover. If k>=1, put A_l=[z_l->x], B_l=[z_l->y]. Counting both blocks over every gap gives

    2+sum_(l=1)^(k-1) [(1-B_l)A_(l+1)+(1-A_l)B_(l+1)] >=2.

For k=0 both exposed words work. Removing the block recovers the old word. Sum over unordered pairs in O to obtain (2). This is the prior capacity argument, not a new independent inequality. QED.

The nonnegative K counts the paths missed by the equality injection. In particular its vanishing is equivalent to |O|<=1. All path counts of nonempty tournaments are positive and odd, so K=d-H(U) is even.

## 3. A finite skeleton for every fixed defect

Assume first that T is strong. List the SCCs of U in dominance order

    C_1 => C_2 => ... => C_t.

Then O meets C_1. Indeed, any path from r to C_1 must first enter U in C_1: arcs between components only point forward. Conversely O meeting C_1 makes T strong: r reaches all of U, all of U reaches a, and a->r. Let j be the last component meeting O.

Delete every singleton component C_i with i>j. Denote the remaining U by U_0, with the inherited set O, and the corresponding tournament by T_0. Let s be the number deleted, h=H(U)=H(U_0), k=K(U,O), and A=L(U_0,O).

**Theorem 2 (exact compression and converse expansion).** One has

    K(U_0,O)=k,
    L(U,O)=2^s A,
    e_a(T)=2^s A-k,
    H(T)=2^(s+1) A+2h-d,        d=h+k.                (3)

The marked tournament T_0 is strong and has the same defect d. Conversely, adding any s singleton SCCs anywhere after the last O-containing component of U_0, with r losing to them, produces a strong sole-exit tournament satisfying (3).

Proof. Remove one such singleton z. In every two-cover counted by K, the first path ends in O, so z cannot occur there: its component is later than every component meeting O, and no path travels backwards in the component order. It must lie in Q. Removal and insertion in the component order are inverse, so K is unchanged.

In a cover counted by L, Q starts in O. Vertex z can lie in either P or Q, independently. Removing z cannot empty Q because z is not its O-start. Conversely insert it in either path at the unique position prescribed by component order. All cross-component arcs are forward. Even when one path has no vertices before or after z, or P was empty, this gives one valid insertion in each bin. Thus L doubles. H(U) is unchanged by singleton-component removal, by SCC multiplicativity. The last component meeting O survives and O still meets the first remaining component, so strongness persists. Iteration and (1) prove (3), as does the reverse construction. QED.

**Theorem 3 (effective skeleton bound).** For positive odd d, every T_0 in Theorem 2 has order at most

    B(d)=d+floor(log2(max(1,d-1)))+4.                  (4)

This bound is deliberately coarse, not optimal.

Proof. The product of H(C_i) over the nontrivial components is h<=d. Each nontrivial strong component has H(C_i)>=|C_i|>=3. Consequently their total number of vertices is at most h<=d: for numbers at least three, their sum is at most their product if there are at least two; the one-component case is immediate. If there are none, this total is zero.

We bound the singleton components retained in U_0. If j=1 there is at most one, namely C_1 itself when it is a singleton. Suppose j>1. Pick x in O intersect C_1 and y in O intersect C_j. If b singleton components lie strictly between C_1 and C_j, then

    K_yx(U) >= 2^b.                                  (5)

To construct these covers, put a Hamiltonian path of C_1 starting at x in Q, a Hamiltonian path of C_j ending at y in P, and put a fixed Hamiltonian path from every nontrivial component between them in P. Each of the b singleton components between them may independently go to P or Q. Put every component after C_j in Q, using a fixed Hamiltonian path per component. Concatenate within each path in component order. Both paths are nonempty and have the required boundaries. Their allocations are distinct. Thus 2^b<=k=d-h<=d-1, and b<=floor(log2(d-1)). Including C_1 and C_j contributes at most two more singletons. No singleton remains after C_j by definition. Adding the nontrivial vertices and a,r proves (4). When j=1, the separate bound d+3 is smaller. QED.

**Corollary 4 (finite affine-dyadic spectrum).** At any fixed positive odd defect d, the endpoint counts of strong sole-exit tournaments are a finite union of rays

    {A 2^s-k : s>=0},                                 (6)

where every (A,k) is obtained from an actual strong marked tournament of order at most B(d), and k=d-h with h a positive odd integer. Their total-path counts form the corresponding finite union

    {2A 2^s+2h-d : s>=0}.                             (7)

One may enumerate all strong sole-exit marked tournaments of orders at most B(d), retain those of defect d, compute their h,k,A, and output these rays. This is an exact finite algorithm, not an order-cutoff conjecture. It need not be computationally practical; duplicate rays may occur. Since A>0, each ray contributes O_d(log X) values at most X. Thus the spectrum at fixed d has cardinality O_d(log X).

No assumption of a finite-state numerical quotient or a local scalar growth law is used. The finite state is the bounded physical skeleton, and the unbounded coordinate is the number of singleton components removed.

## 4. Arbitrary tournaments and bounded-defect conversion barriers

Let T be any tournament with a sole-exit mark a and e_a(T)>0. The mark lies in its final SCC C. That component is nontrivial, and the sole exit r lies in C. Write M for the product of the path counts of the earlier components. Then

    H(T)=M H(C),   e_a(T)=M e_a(C),   d(T,a)=M d(C,a).   (8)

Thus M is a positive divisor of d. For every fixed d there are only finitely many possibilities for M and d(C,a), so Corollary 4 and scaling by M again give a finite union of affine dyadic rays and O_d(log X) endpoint or total-path values. Only those M realized by an actual prefix are retained; (8) is not permission to introduce arbitrary integer factors. A finite exact prefix test is available: strip its singleton components, leaving at most M vertices by the same product argument as in Theorem 3.

**Corollary 5.** Fix D. The total-path spectrum of tournaments admitting a positive sole-exit endpoint a with

    1 <= H(T)-2e_a(T) <= D

has O_D(log X) values at most X. In particular it cannot contain every sufficiently large odd integer. Equivalently a proposed construction that preserves a cycle register as such a sole-exit endpoint and keeps its additive excess over twice the register uniformly bounded cannot prove an odd tail.

The sole-exit condition on the output endpoint is essential to the claimed scope. Session213 excludes exact excess one for every endpoint, but neither that theorem nor this result rules out unbounded excess, output endpoints with several exits, changing the register, or unrelated constructions. The original deletion of a sole-exit vertex in a cycle tournament need not leave a *sole-exit marked endpoint* in the output, so Corollary 5 applies only when that additional condition is met.

## 5. Complete defect-three endpoint spectrum at a sole exit

For a transitive chain U=(y_0,...,y_(ell-1)), define

    w(1)=1,   w(t)=2^(t-2) for t>=2.

For i<j, the exact two-cover capacity is

    K_(y_i,y_j)+K_(y_j,y_i)=2w(j-i).                  (9)

Indeed K_(y_j,y_i)=2^(j-i-1): all chain vertices before i are forced into P, all after j into Q, and each between i and j chooses either bin. The opposite orientation has one cover if j=i+1 and none otherwise. Consequently

    d=1+2 sum_{{i,j}: y_i,y_j in O, i<j} w(j-i).      (10)

This is a direct pair-cover identity; the appearance of the familiar one-apex weights is inherited combinatorics.

**Theorem 6.** Over all finite tournaments with a sole-exit mark a, e_a>0, and H=2e_a+3, the set of endpoint counts is exactly

    {1}
    union {3*2^q : q>=0}
    union {3*2^q-2 : q>=1}
    union {5*2^q-2 : q>=1}.                           (11)

Equivalently the total-path spectrum is exactly

    {5}
    union {6*2^q+3 : q>=0}
    union {6*2^q-1 : q>=1}
    union {10*2^q-1 : q>=1}.                          (12)

Proof. First suppose T strong. By (1), either (h,k)=(1,2) or (3,0).

If h=1, U is a transitive chain. Strongness requires y_0 in O. Formula (10) forces precisely two vertices in O, at distance one or two. Thus O={y_0,y_1} or O={y_0,y_2}.

- O={y_0,y_1}, ell=2: e_a=1. For ell>=3, one has e_a=3*2^(ell-2)-2.
- O={y_0,y_2}, ell=3: e_a=3. For ell>=4, one has e_a=5*2^(ell-3)-2.

For a direct count, L(U,O)=sum_(y_i in O) 2^(ell-1-i): a cover whose Q begins at y_i forces earlier vertices into P and assigns every later vertex freely. Here k=2 and e_a=L-2. This proves both formulas, including their short boundaries.

If h=3, the SCC decomposition of U is

    TT_p => C3 => TT_q,       p,q>=0.                 (13)

To see that no other components are possible, their H values multiply to three, and any nontrivial strong component has order at most its H. The unique nontrivial component must therefore be C3. Since k=0, capacity positivity forces O={y}; strongness places y in the first SCC.

- If p=0, y is a triangle vertex and e_a=L=3*2^q. There are exactly three two-covers of C3 with the second path starting at y: its vertex set is {y}, {y,z} for the unique outneighbor z of y, or all of C3. Each following singleton chooses either path.
- If p>=1, y is the first chain vertex and e_a=L=12*2^(p+q-1)=3*2^(p+q+1). The first vertex is forced into Q, each other chain vertex chooses either path, and C3 has twelve ordered two-covers allowing empty bins (six with one empty bin and six with both nonempty).

These endpoints are all in {3*2^q:q>=0}, and that entire ray already occurs when p=0. This proves strong necessity and provides every listed value by explicit construction.

For nonstrong T, (8) gives M=1 or M=3. M=1 only adds a transitive prefix. If M=3, the final component has defect one; Session213's equality classification, or (1) and (10) with h=1,k=0, gives e_a(C)=2^q. Hence e_a(T)=3*2^q, already listed. Conversely these have the strong realizers above. This proves (11)-(12). QED.

## 6. What this says about the original target

The exact equality H=2c+1 only permits power-of-two endpoint registers c, as Session213 proved. Allowing excess three at a sole-exit output mark creates actual additional infinite families, including H=23,47,95,... and H=19,39,79,..., but still only logarithmically many values up to X. Theorem 3 proves that *every fixed excess* has this same sparse form under the sole-exit condition, so increasing a fixed error allowance cannot repair the approach.

For instance c=5 is not in (11), so no sole-exit marked tournament has e_a=5 and H=13. This is not an exclusion of H=13 as an unmarked tournament count; such tournaments exist. The finite verifier includes that hostile explicitly.

This session does not extend the inherited coverage bank, prove a positive-density path spectrum, classify arbitrary endpoint defects, or settle Chan–Pak. It supplies a bounded physical skeleton, an exact unbounded inverse at fixed defect, a complete next-defect spectrum, and a precise route that must be escaped.

## 7. Reproduction

`python3 verify.py` and `python3 -O verify.py` produce identical JSON. The checker separately validates the exact cover identity, SCC compression, the coarse bound, the transitive pair weights, explicit infinite-ray instances over a finite parameter grid, the defect-three classifier, and the H=13 hostile. It uses explicit exceptions. See RESULT.json for exact finite universes. No finite tests substitute for the proofs above.
