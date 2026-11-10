---
id: THM-5121
title: "Bounded-deletion dichotomy, binary shadowing and physical ratio control"
status: "PROVED, independently audited; exact scalar coverage remains OPEN."
session: Session246, 2026-11-10
source: 05-knowledge/results/session246_unifying_principles_20261110/descent/PROOF.md
audit: 05-knowledge/results/session246_unifying_principles_20261110/audit/AUDIT.md
---

# Real-state control does not erase arithmetic gaps

Use the inherited physical O/Z/J source grammar from t=(1,0,2,1).
States (D,E,U,d) have 0<=E<=d, D>=d, U>=D and U>=2d.

**Deletion dichotomy.** On any infinite physical trajectory the deletion
coordinate d is bounded iff the word has finitely many root-changing J's.
If d_i is the deletion count immediately after its i-th J, then

    d_(i+2)>=d_(i+1)+d_i.

Indeed O preserves D,d; Z increases D by U−E>=d; J makes d the old D and
has the same positive D increment. This gives Fibonacci growth along the
root-change subsequence. This is about a supplied word, not an invariant
minimum number of J's over all presentations.

**Uniform shadow theorem.** On an O/Z tail let d be fixed, a=U_0−d>0,
x_n=D_n/U_n and epsilon_n=0 for O,1 for Z. Starting y_0=x_0, put

    y_(n+1)=(y_n+epsilon_n)/2.

For every tail and n>=0,

    U_n>=d+2^n a,
    |x_n−y_n|<=n d/(2^n a),
    y_n=2^(−n)x_0+sum_(j=0)^(n−1) epsilon_j 2^(j−n).

The O step is exact. The Z error is
[d(D+U)−2EU]/[2U(2U−d)], of magnitude at most d/(2U−d).
The displayed bound follows by contracting and summing the errors.

For t in [0,1] and r>=1, read the r binary digits of
min(2^r−1,floor(2^r t)) least-significant first. The resulting actual suffix
has |D_r/U_r−t|<=(1+r d/a)/2^r. This is constructive ratio control, not a
selector for a prescribed integer Hamilton count.

The ratio converges iff the tail letters are eventually constant: to zero
for eventual O and one for eventual Z. A repeated binary block of length l>=1
has endpoint limit A/(2^l−1), where A=sum epsilon_j 2^j; intermediate phases
follow the binary affine recursion. Mixed blocks have nonconstant cycles.

There exists a single deletion-one physical trajectory with ratio limit set
[0,1]: concatenate all finite binary words repeatedly in increasing length.
Nevertheless its balanced counts G=D+2U belong to the inherited incomplete
deletion-one image. Session199 hidden-structure Theorem3 bounds that image
through N by cN+O(log N), with

    c=sum_(j>=0)1/(2^(j+2)+1)<9/20.

The density bound is credited to that source, not reproved here. The new
conjunction shows explicitly that full continuum limit-set coverage and
uniform physical ratio control coexist with positive-density integer gaps.
The original word and full integer state are essential arithmetic sidecars.

[Proof, positive controls, source recovery and standalone exact verification](../../05-knowledge/results/session246_unifying_principles_20261110/descent/PROOF.md).
