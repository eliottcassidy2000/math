# A one-birth ray step needs defect at least72

**PROVED corollary of inherited identities and the exact finite census
below.** For every k>=1, start with any actual tournament of Hamilton
count h_k=405*2^k+71, permit any preparation preserving that full count,
and then add one vertex y beaten only by the selected source root r.
If the newborn graph has the desired next count h_(k+1)=2h_k-71,
its actual indegree-one readiness defect is at least72. In particular
this operation cannot give a unit-ready next source, even after an
arbitrary count-preserving change of core.

The elementary bound gives at least71. The finite H71 head rules out
equality for the entire ray. This is an operation obstruction, not a
missing Hamilton count, not a restriction on multi-birth returns, and
not a Chan--Pak proof.

## 1. The inherited inequality and the exact source budget

Write the actual source profile as(h,s,q,F), where s counts Hamilton
paths starting at r, q=H(T-r), and F counts ordered spanning two-path
covers with the first path nonempty and ending at r. The second bin
may be empty. The inherited newborn identity is

    H_new=h+F-s,       readiness defect_new=s.

Cutting each Hamilton path immediately after r injects its h paths
into the F covers, so F>=h. This bound and its newborn consequence
H_new>=2h-s are already explicit in
[ready-return obstruction, O2](../../ready_return_pentagon_20261103/OBSTRUCTION.md).
They are not new here.

More precisely, let d be the indegree of the selected source root and
let Delta=F-h-q+s. The inherited
[weighted-cut identity](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse)
gives

    q>=s,       Delta>=d(d-1),       Delta=0 iff d<=1.

For a requested one-birth count2h-c, the exact source budget is

    F-h=s-c,
    q+Delta=2s-c,
    s>=c+d(d-1).                                      (1)

The first two statements are equivalent to the target equation; the
last is a necessary consequence. They retain the start channel lost
by merely marking an indegree-one output.

For our ray c71, a unit-ready newborn would require s1 and would exceed
the requested count by at least70. No same-H preparation changes this
inequality. The possible sharp boundary s71 would require F=h,
q=s71 and Delta0. Thus d<=1. If d0, the source root is a global source,
so s=h; but h_k>=881 differs from71. Therefore the boundary would force

    d1,       source profile(h_k,71,71,h_k).            (2)

The equality F=h is possible for strong sources: the directed triangle
has profile(3,1,1,3). It must not be replaced by a false strict inequality.

## 2. Why the sharp boundary reduces to strong H71 cores

Suppose(2) exists. Delete its root y and call the remaining tournament
U. Let a be the sole in-neighbor of y. Then the source is exactly the
newborn construction from(U,a), with

    H(U)=71,       S_U(a)=0,       F_U(a)=h_k-71=405*2^k. (3)

These statements follow directly from the newborn profile map:
s_y=H(U)-S_U(a), F_y=H(U)+F_U(a).

Hamilton counts multiply over ordered strongly connected components.
Since71 is prime, U has exactly one nontrivial strong component S with
H(S)=71; every other component has count1 and hence is a singleton.
Thus U=TT_A -> S -> TT_B for some A,B>=0.

The inherited SCC cover factorization
[O1](../../ready_return_pentagon_20261103/OBSTRUCTION.md)
now gives three possibilities for F_U(a):

* a lies before S: 71 times a nonnegative power of2;
* a lies inside S: F_S(a) times2^A;
* a lies after S: P_S(2) times a nonnegative power of2.

The first has odd part71 and cannot equal405*2^k. In the other two
cases a necessary condition is that some F_S(a), or P_S(2), has odd
part405. This argument even allows more weak-parent placements than
strongness of the final source permits, so it is a safe exclusion.

The inherited published strong-order bound is H>=75 at order9, and
is nondecreasing thereafter; see the citation in O2 above. Therefore
every possible S has order at most8.

## 3. Exact finite head, with its completeness boundary

The input is the inherited complete6880-class order-eight census
[tournaments8.txt](../../../../session212/cycle-spectrum/tournaments8.txt),
whose provenance and checks are in
[Session212 PROOF.md](../../../../session212/cycle-spectrum/PROOF.md).
This lane does not regenerate that isomorphism census. It independently
counts all its Hamilton paths, selects its18 H71 presentations, and
extracts their nontrivial SCCs. Every smaller strong H71 core appears:
append enough transitive sink vertices to obtain an order-eight class,
then recover its nontrivial component.

The result is eleven strong core presentations, seven of order7 and
four of order8. Their complete outgoing masks and all81 marked profiles
are in [birth_budget_receipt.json](birth_budget_receipt.json). The
root-cover lists are:

|order|P_S(2)|F_S at every root, in the receipt's vertex order|
|--:|--:|---|
|7|1012|287,167,105,118,102,90,72|
|7|1012|241,184,115,97,149,77,78|
|7|1012|273,139,167,103,95,73,91|
|7|1012|291,137,109,143,100,89,72|
|7|1072|292,147,113,172,79,119,79|
|7|1012|201,193,137,104,152,79,75|
|7|1012|207,189,112,168,110,74,81|
|8|1976|535,243,517,167,112,75,182,74|
|8|1976|531,287,505,190,128,103,87,74|
|8|1976|523,349,236,153,404,85,84,71|
|8|1976|595,404,268,217,177,86,71,87|

No F equals405, and every F is at most595<810. Hence no F has odd
part405. The possible P2 values1012,1072,1976 have odd parts253,67,247,
respectively. None has odd part405 either. This rules out(3), hence(2),
at every k>=1. Combining with(1), the newborn defect must be at least72.
QED.

## 4. Complete strong balanced-t71 templates and their count support

There is a positive classification behind the exclusion. Let F71 be the
finite set of root-cover values in the table and P71={1012,1072,1976}.
The complete set of full counts of strong marked tournaments with
profile(H,71,71,H) is exactly

    {71+2^a f : f in F71, a>=1}
       union {71+2^c p : p in P71, c>=0}.              (4)

These are all actual graph families, not just necessary numerical forms.
The receipt compresses their overlaps into odd-part/minimum-exponent
pairs, each representing all amplitudes H-71=oddpart*2^e with e at
least its listed minimum. This gives a disjoint exact dyadic-ray
description of the count support.

For the proof, delete the marked root as in(3), now allowing arbitrary H.
Again U=TT_A -> S -> TT_B with H(S)=71, and the sole predecessor a
of the deleted root has S_U(a)=0. The newborn tournament is strong
exactly when a lies in the final SCC of U: every old vertex must reach
a before it can enter the newborn, and this condition also suffices
because the newborn beats every old vertex except a.

If B0, a lies in S. All81 source-core start counts are positive, so
S_U(a)=0 forces A>=1. The amplitude is2^A F_S(a). Conversely prepend
TT_A to S and birth at any selected core root; the newborn graph is
strong and has the required balanced profile.

If B>=1, a is the final singleton of that tail and the amplitude is
2^(A+B-1)P_S(2). Every exponent c>=0 is attained by taking A0 and
B=c+1. These sources too are strong and balanced. This proves both
directions of(4), and gives the full source-template description
(S,A,B,a) up to isomorphism of the marked core. Duplicate count values
do not identify distinct templates.

## 5. Reproduction and the useful next question

```text
python -B 05-knowledge/results/changing_core_returns_20261113/returns/birth_budget.py
python -O -B 05-knowledge/results/changing_core_returns_20261113/returns/birth_budget.py
```

The receipt pins the inherited census by SHA256, records all source
masks, and checks the weighted-cut inequalities at every selected root.
Triangle and H71 equality cases are positive controls for F=h; the two
F404 entries are near misses, not405 certificates. Normal and optimized
receipts agree exactly.

There is an actual repeatable neighboring construction. Take either saved
H71 core/root with F404. Prepend a global source, retaining that marked
core vertex. The resulting weak parent has H71,start0,F808. A newborn
birth therefore gives a strong balanced source(879,71,71,879). Repeated
keep-root births R construct

    (404*2^k+71,71,71,404*2^k+71),     every k>=1.

The receipt stores and independently counts the first four literal
graphs. The all-k statement follows from the exact R map. The gap from
the requested405-ray is exactly2^k; the missing single core cover at405
is amplified by the same clock. These outputs have defect0, so this is
a repeatable neighboring scalar family, not a unit-ready return.

The attempted positive construction was meaningful: an actual balanced
profile(h_k,71,71,h_k) would let the physical keep-root birth R repeat
the entire scalar ray, since R sends(h,t,t,h) to(2h-t,t,t,2h-t).
The finite head proves why that particularly economical proposal fails.
The surviving route must pay a larger start/defect budget, use a
different multi-birth operation, or change the selected target equation.
This does not yet give the missing repeatable unit-ready return.
