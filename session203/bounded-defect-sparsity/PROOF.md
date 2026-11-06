# All bounded physical Z defects have sparse profile support

Status: PROVED / independently audited in the stated scope.

## 1. Statement and exact scope

For a nontrivial rooted tournament U with positive root-start count s, use h,e,L,R,d,b and B=h+L-s as in the frozen dependencies. Its physical Z defect is

    D=B(ZU)-2B=R+b-e-B >= s >= 1.

**Theorem.** For every fixed positive integer K, there is a finite constant A_K such that, for every real X>=1,

    #{B<=X : some actual positive-start nontrivial rooted U has D<=K}
      <= 1 + A_K (floor(log2 X)+2K+5)^(2K+1).       (1)

Consequently all bounded-absolute-defect profiles have zero density among odd B. The bound concerns all physical sources, not a selected construction language. No claim of uniformity in K is made. F(8,137) remains open: its witnesses are allowed to have unbounded D.

The proof uses the nonnegative defect decomposition from SHARP_DEFECT_PROOF.md, the two-terminal cover bound from COVER_EQUALITY_PROOF.md, elementary tournament SCC endpoint facts, and the inherited physical edge extension E(U) with H(E(U))=B. The edge extension adds a vertex z with v->z and z->U-v. The original singleton source is excluded from the definition of physical D here.

## 2. Removing at most two exceptional vertices

Let k be the root outdegree. For positive start, k>=1. The exact outgoing-cover identity and nonnegative defect decomposition give

    D-s = f+2p+(b-L-d),
    b-L = sum_{a: v->a} L_a(U-v),
    f>=0, p>=0, L_a(U-v)>=d=H(U-v).

If k>=2, these imply

    D>=s+(k-1)H(U-v),
    H(U-v)<=K,  k<=K.                            (2)

Thus U is obtained from a tournament T with H(T)<=K by adding v, with at most K exceptional outgoing neighbors; all other old vertices beat v.

Suppose k=1, with unique outneighbor a. The rooted edge U=(v->a) has (B,D)=(1,3); reserve its single B value as the first term in (1). Otherwise T=U-{v,a} is nonempty. Let I={x in T: x->a} and put t=H(T). Write K_T(I) for the number of ordered covers of T by two nonempty paths (Q,P) with last(Q) and first(P) both in I. The endpoints are distinct.

The exact bad-cover decomposition in the dependency gives

    L_a(U-v)-H(U-v)
      = sum_{x in I} start_x(T)+K_T(I)
      = t-s+K_T(I).

For unique root outdegree one, f=0 and p=end_a(U-v)=sum_{x in I}end_x(T). Hence exactly

    D=t+2 sum_{x in I} end_x(T)+K_T(I).             (3)

The two-terminal cover lemma gives K_T(I)>=binom(|I|,2). Since t>=1, (3) yields t<=K and |I|<=K: for |I|>=2, binom(|I|,2)>=|I|-1, and the cases |I|<=1 are immediate. Thus U is obtained from T with H(T)<=K by adding v,a, with v->a, T->v, and a dominating all but at most K marked vertices of T.

The nonempty-T qualification is necessary in (3). The rooted edge was explicitly removed before applying it.

## 3. Finitely many fixed cores and linearly ordered chain modules

If H(T)<=K, its strong-component Hamilton counts multiply to at most K. Every nontrivial strong component has order at most its Hamilton count and count at least three. The sum of the counts of the nontrivial components is at most their product, unless there is just one, when equality holds. Hence the total number of vertices in all nontrivial strong components is at most K. There are only finitely many possible ordered lists of these components.

Retain all those vertices and all marked exceptional vertices from Section 2. At most K additional retained vertices are singleton strong components. After deletion of these retained blocks, the remaining singleton components of T fall into at most 2K+1 ordered runs. Each run is a transitive tournament. Every vertex outside a run has the same orientation to all its vertices. The runs themselves are totally ordered: all edges go from an earlier run to a later one.

For fixed K there are only finitely many retained patterns, including marks and their orientations. U is a member of one of these finite patterns with arbitrary nonnegative lengths n_1,...,n_r of its transitive runs, where r<=2K+1. Adding the inherited auxiliary vertex to form E(U) preserves all these properties. Its fixed core has at most

    c<=2K+3

vertices. Its only unbounded parts are the r linearly ordered transitive runs. We may allow length zero even when a particular original description required a nonempty run, since this only enlarges the set used for an upper bound.

## 4. Exact positive counting on a fixed pattern

Fix any tournament pattern with c fixed core vertices and r transitive variable runs, linearly ordered among themselves and homogeneous toward every core vertex.

A Hamilton path orders the c core vertices. In each of its c+1 gaps (before, between, or after these core vertices), it visits a list of variable runs in increasing run order, with at most one nonempty segment from each run. There are finitely many such skeletons. Discard any skeleton whose consecutive edges contradict the pattern.

For a valid skeleton let q_i be the number of occurrences of run i. Then 0<=q_i<=c+1. Assigning its n_i labeled linearly ordered vertices to these q_i nonempty occurrences determines the path uniquely: each occurrence lists its assigned vertices in increasing order. The number of assignments is

    Surj(n_i,q_i)=q_i! S(n_i,q_i),

the number of onto maps to q_i labeled bins, with Surj(0,0)=1 and Surj(n,0)=0 for n>0. Therefore the exact Hamilton count is the finite nonnegative sum

    H = sum over valid skeletons product_i Surj(n_i,q_i).   (4)

This is a literal decomposition, not an inclusion-exclusion estimate. Distinct skeletons and assignments give distinct Hamilton paths, and every Hamilton path supplies exactly one of each.

If n>=q>=2, assigning the first q vertices one to each bin and every remaining vertex freely to the first two bins proves

    Surj(n,q)>=2^(n-q).                           (5)

Suppose the count in (4) is at most X and some n_i exceeds N=floor(log2 X)+c+1. Any active skeleton with q_i>=2 would contribute at least 2^(n_i-q_i)>X by (5), a contradiction. An active skeleton also cannot have q_i=0 because n_i>0. Thus every active skeleton has q_i=1, whose factor is exactly one.

Reduce this n_i to N. Since N>=c+1, feasibility of every q_i<=c+1 is unchanged: an inactive skeleton remains inactive, and all active factors at this run remain one. The total Hamilton count is unchanged. Repeat independently for every oversized run. Every count <=X therefore has a representative with

    0<=n_i<=floor(log2 X)+c+1.

One fixed pattern supplies at most (floor(log2 X)+c+2)^r different counts up to X. This includes runs whose length is completely invisible to the count; those are normalized rather than incorrectly bounded in the original source.

## 5. Finish and implication for inverse methods

Apply Section 4 to the finitely many E(U) patterns from Section 3. Since B=H(E(U)), c<=2K+3, and r<=2K+1, summing their bounds and restoring the single rooted-edge value proves (1).

For any fixed staggered width J, the finitely many shifts X=B+2j, 0<=j<=J, of these bounded-defect supports still have density zero. Thus a proof that every sufficiently large odd base has bounded-start staggered witnesses cannot use only sources whose absolute Z defect stays below a fixed constant. It must allow unbounded defect along its target family. This is a limitation on that proposed mechanism, not a Hamilton-spectrum gap and not a refutation of F(8,137).
