# Quantitative root-response defect via prescribed pair covers

**INHERITED bound; alternative proof independently audited.** Session244, 10 November 2026.
No full-count coverage theorem is claimed.

The bound is already implied by the stronger weighted inequality in
THM-4875, `THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, part(A1),
which credits THM-4873 and `two_entry_funnels_20261103/READY.md`, equation(7).
That inequality has the extra weighted term
`2 sum_{{a,b} subset O} H(T-{r,a,b})`. The following SCC pair-cover proof
and parity decomposition are independent recovery, not a new stronger bound.
The saturation characterization is likewise inherited from Session231
`descent/ROUND2.md`, section2. These are retained as verified dependencies
for the new synchronized compiler and root-profile analysis.

Let T be any nonempty tournament rooted at r. Write h=H(T), U for ordered
spanning two-path covers (P,Q) with P starting at r and Q possibly empty,
O=N+(r), k=|O|. Let e_O(T-r) count Hamilton paths of T-r ending in O,
with value zero for empty O. Then

    U-h >= k(k-1)+e_O(T-r).                         (1)

Thus a nontrivial strong tournament with U=h has root outdegree one.
The singleton is a genuine exception to that outdegree conclusion.

## 1. Prescribed unordered marks

**Lemma.** For any tournament A with distinct marks a,b, a spanning pair
of nonempty paths exists with one path starting at one mark and the other
ending at the other mark.

Write the ordered strong components of A-a as C1=>...=>Cs, with b in Ci.
Every vertex in a strong component can start or end a Hamilton path: cut
a Hamilton cycle there; singleton components are immediate.

If i=s, use (a) and a Hamilton path of A-a ending b. Otherwise let C=Ci+1.
If a beats some c in C, use P=(a), then a Hamilton path of C starting c,
then all components after C; the other path runs through C1,...,Ci and
ends b. If a beats no vertex of C, then C=>a. Use one path starting b,
through Ci and every component after C, and the other through C1,...,Ci-1,
then C, then a. Empty prefixes/suffixes are omitted. This proves the lemma.
Strong-component ordering and strong-tournament Hamiltonicity are inherited
standard facts; no literature-priority claim is made.

## 2. Parity strengthens existence

Let K_ab count spanning pairs in A whose first path starts a and second
path ends b. Neither count K_ab,K_ba depends on the arc between a and b,
since the two marks lie on different paths. The lemma gives K_ab+K_ba>0.

Compare the tournaments differing only in the arc on {a,b}. Paths not
using that pair consecutively coincide. Joining the second path to the
first identifies K_ab with paths using b->a, and K_ba with paths using
a->b. The difference of Hamilton counts is therefore K_ab-K_ba up to sign.
Both full counts are odd by Redei's theorem, so their difference, and also
K_ab+K_ba, is even. Consequently

    K_ab+K_ba >= 2.                                  (2)

## 3. Exact defect decomposition

Split each Hamilton path immediately before r. Its distinguished suffix
starts r. A cover is in the image iff its complementary path Q is empty
or ends in N-(r). Thus U-h counts covers with Q nonempty ending in O.
The distinguished singleton path(r) contributes e_O(T-r).
Every remaining distinguished path begins r,a for a in O; Q ends at a
distinct b in O. Removing r gives exactly the pair covers counted above.
Different unordered pairs give disjoint classes, and they exhaust them:

    U-h = e_O(T-r) + sum_{{a,b} subset O}(K_ab+K_ba).

Equation(2) proves(1). When O={a}, U=h iff no Hamilton path of T-r ends a,
equivalently a is outside its terminal strong component. Strongness forces
a in the first component. This recovers the exact saturated-state structure
of ../composition/STRUCTURE.md, whose reachability proof is independent.

The bound is attained by the directed triangle and by rooted bitmask graphs
[8,1,3,6], root3, with(h,U,k,e_O)=(5,7,2,0), and
[16,1,3,7,14], root4, with(h,U,k,e_O)=(9,15,3,0).
No claim of sharpness for every k is made. The pair-indexed even-positive
contributions are the physical sidecar lost by a scalar response inequality.

## Reproduction

Compile verify.cpp with C++17, -Wall -Wextra -Werror, then run it. It checks
all 33,867 labelled tournaments of orders1..6 and all202,013 rooted cases,
including every outneighbor-pair coefficient in the exact identity. A
separate audit uses permutations through order5. The --fail flag deliberately
throws. Normal and -DNDEBUG builds must produce identical JSON. Proofs above
supply the all-order quantifiers; the finite runs certify the implementation.
