# Independent cycle-side audit of the root's determinant inequality

**Status: proof accepted, exact finite checks below.** This is a cross-lane audit of the root's proposed global cycle-expansion proof, requested after the carrier packet was completed. The carrier worker supplies an independent sign/empty-sidecar review, a weaker-induction simplification, and the equality classification. No external novelty claim.

Let A be the loopless adjacency matrix of a tournament T, and put

    H=H(T),  u=det(I+A),  v=(-1)^n det(A).

Set H(empty)=u(empty)=v(empty)=1. The full inequality is deliberately **restricted to nonempty tournaments**: the empty tournament has1<1+1.

## 1. Determinant kernel positivity

For nonempty T, B=I+A has B+B^t=I+J. The matrices I+sA, 0<=s<=1, have positive definite symmetric parts, so are nonsingular; continuity from s=0 gives u>0. Put z=B^(-1)1 and t=1^t z. Then

    2t = z^t(B+B^t)z = ||z||^2+t^2.

Since z is nonzero, 0<t<2. Tournament complementation gives `-A=B^t-J`; the determinant lemma therefore gives

    v=u(1-t),  hence u>|v|.

For the empty sidecar, u±v are0 and2. Thus u±v>=0 for every sidecar, including empty. This proves the exact positivity needed by the combinatorial induction without a numerical eigenvalue argument.

## 2. Independent sign audit

Adjoin x with out-neighborhood S in the old tournament T. Every new directed cycle through x corresponds uniquely to a simple directed old path P starting in S and ending outside S. Such a path has at least two vertices. Write R=T-V(P), and split paths into E/O by even/odd cardinality |V(P)|.

A distinguished child-cycle has length |V(P)|+1, and its determinant sign in det(I+A_child) is (-1)^|V(P)|. A full-cover cycle contributes a factor-1 in v regardless of its length. The literal cycle expansions are therefore

    u_child = u + sum_E u_R - sum_O u_R,
    v_child = - sum_(E union O) v_R,
    H_child = H + 2 sum_E H_R.

The third is inherited OCF. Consequently, for epsilon=+1 or-1,

    H_child-u_child+epsilon v_child
      = (H-u)
        + sum_E [2(H_R-u_R)+(u_R-epsilon v_R)]
        + sum_O (u_R-epsilon v_R).                       (1)

No orientation reversal, cyclic-rotation multiplicity, loop, or missing fixed point changes these signs.

## 3. Weaker induction suffices; exact equality boundary

First prove H>=u for all T, including empty, by induction using

    H_child-u_child
      = (H-u)+sum_E [2(H_R-u_R)+u_R]+sum_O u_R.

Every term is nonnegative, since u_R>=1 and all smaller weak gaps are nonnegative. Then every term of(1) is nonnegative by kernel positivity. For **every nonempty tournament** this proves

    H(T) >= det(I+A)+|det(A)|.

This avoids incorrectly assuming the stronger inequality for the empty sidecar.

Moreover, H-u is positive whenever T is cyclic: choose x on a directed triangle; the corresponding two-vertex path gives an E term at least u_R>=1. Conversely a transitive tournament has H=u=1.

If cyclic T has n>=4, choose x outside a fixed directed triangle. Its old parent remains cyclic, so the parent term H-u in(1) is at least1. Both signed gaps are then at least1. A directed triangle itself has(H,u,v)=(3,2,-1), giving equality. Therefore

    H=u+|v| iff T is nonempty transitive or T is a directed triangle.

For cyclic tournaments of order at least4, the strengthened gap H-u-|v|>=1 is sharp by padding a directed triangle with singleton SCCs.

## 4. Exact check

`determinant_audit.py` uses this lane's independently written exact class canonicalizer and Hamiltonian-path DP. It separately implements integer Bareiss determinants, cross-checks them by literal determinant permutation sums through order5, and literally sums all cut-selected simple paths for every parent class through order5. It checks both signed recurrences, every empty sidecar, strict u>|v|, the full inequality and the equality classification on all7,412 nonempty tournament classes through order8.

Reproduce:

```bash
python3 -B determinant_audit.py > determinant.normal.json
python3 -B -O determinant_audit.py > determinant.optimized.json
cmp determinant.normal.json determinant.optimized.json
```

An initial audit harness incorrectly inherited a nonempty-only HP helper's empty value0. The explicit sidecar H(empty)=1 wrapper corrected that implementation before any check was accepted; the mathematical recurrence did not change. This is separate from the genuine failure of the full inequality at empty, which is retained as a hostile control in the output.

## 5. Stronger unweighted defect: exact packing-one boundary

The root subsequently proposed the following refinement. The cycle-side audit accepts its proof, including the sharp strictness mechanism.

Let C_all(T) count all literal directed simple cycles of length at least3, and let C_Ham(T) count those using all tournament vertices. Put C_proper=C_all-C_Ham. These are unweighted cycle counts, not OCF-weighted packings. For every **nonempty** tournament,

    H(T) >= det(I+A)+|det(A)|+C_proper(T).                (2)

Equality in(2) holds **iff T has no two vertex-disjoint directed triangles**. If it has two such triangles, the gap above the right side of(2) is at least2. Packing of arbitrary directed cycles and of directed triangles have the same maximum size: each directed cycle's induced support contains a directed triangle.

### Weak gap pays for every cycle

The weaker bound from Section3 strengthens to

    H(T)-u(T) >= C_all(T)                               (3)

for every T, including empty. Induct using the unsigned weak-gap recurrence there. Every new directed cycle through x corresponds to one path P; its increment is either `2(H_R-u_R)+u_R` or `u_R`. Both are at least1, since the earlier weak inequality applies to R and u_R is a positive integer, including u_empty=1. The old gap pays for all cycles avoiding x. Thus it pays once for every cycle of the child, proving(3).

### Both signed margins pay for every proper cycle

Use signed recurrence(1). The old-parent term is at least C_all(parent) by(3), and these are precisely the child cycles avoiding x; all are proper in the child. A path P has nonempty remainder R iff its new cycle through x is proper. For these paths, strict determinant domination makes `u_R-epsilon v_R` a positive integer, hence at least1. The additional even-path term `2(H_R-u_R)` is nonnegative. Empty-R paths correspond to Hamiltonian cycles and contribute nonnegative terms, but need not contribute1. Consequently

    H_child-u_child+epsilon v_child >= C_proper(child)

for both epsilon values. Taking their minimum proves(2).

### Disjoint triangles force a surplus of two

If C1,C2 are disjoint directed triangles, select x in C1 and use the two-vertex old path completing C1. This is an even path; its nonempty remainder R contains C2, so by(3), `H_R-u_R>=1`. Its summand in either signed recurrence is therefore at least3, whereas counting that proper cycle needs only1. All other assigned terms remain valid lower bounds. Both signed margins, and hence their minimum, exceed C_proper by at least2.

Conversely, suppose no two directed triangles are disjoint. Then no two directed cycles are disjoint. Write O/E for the numbers of odd/even cycles. The determinant expansions contain no multi-cycle terms, the full-cover expansion has at most one cycle, and OCF has only single odd cycles. Therefore

    u=1+O-E,  v=-C_Ham,  H=1+2O,
    H-u-|v|=O+E-C_Ham=C_proper.

These identities also cover nonempty transitive tournaments (all cycle counts zero), orders1 and2, and a cyclic SCC padded with singleton SCCs. They prove the converse equality direction. Two dominance-ordered directed triangles give a sharp strictness example: `(H,u,v,C_proper)=(9,4,1,2)`, leaving surplus2.

### Extended exact check

The same audit now enumerates all directed cycles, not only odd ones, on every one of its7,412 complete classes through order8. It checks(3), (2), the equality iff, the surplus>=2 implication, and the packing-one determinant expansions independently of the proof. The original470 literal insertion checks and permutation-determinant controls remain active. The final normal/optimized output hashes, which supersede the earlier shorter determinant receipt, are recorded below after both executions complete. The carrier classification, its source pins and its primary output are unchanged.

Final extended paired receipt (normal and optimized byte-identical):

    SHA256 f98d25d324026aa32ebadfe207b295acce380aaaf1284e23c95d5b9d92be28db

The packing-at-least-two class counts at orders6,7,8 are25,386,6734; the corresponding minimum strict surpluses are2,3,3. All remaining classes satisfy the exact packing-one equality. These finite minima are not asserted at larger orders.
