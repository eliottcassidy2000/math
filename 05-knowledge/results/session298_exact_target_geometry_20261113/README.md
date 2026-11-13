# Session 298: exact-target marked preparation

## Proved: neutral preparation at every target 2^m+1

Use the physical directed triangle with parts C=TT_m, X={x}, Y={y}, arcs C→X→Y→C. Fix a named root in C; change only the total order inside C. This preserves the whole unrooted tournament type, hence H=2^m+1, but changes the rooted response. This is the all-order extension of the marked preparation in [the H17 controller](../eight_nine_twentyseven_20261112/constructors/SIX_RETURN.md), not a claim that all targets already have such sources.

With l predecessors and r successors of the root in C (l+r=m−1), the marked state (h,s,q,F) is
```
h=1+2^m, s=2^r, q=1+2^(m−1), F=2^(l+1)+2^(r+1)3^l−[r=0].
```
Here m≥2, s is the root-start count, q counts deletion paths, and F counts ordered spanning pairs (P,Q) with P nonempty ending at the root, Q possibly empty. For h the odd-cycle packing formula gives h=1+2(2^(m−1)): an odd cycle must use x,y and an odd nonempty resource subset. Deletion applies the same formula. Starting at the root, choose any subset of its successors before x,y, then use remaining resources; when the subset is all successors and l=0, the second resource run is empty. This gives s.

For F, if neither x nor y is in P, choose a subset of the l predecessors for P; summing H on its complement yields 2^l+2^r3^l−[r=0] (empty resource complement has H=1, not 2). If only y is in P it must precede its increasing run, giving 2^l; x alone is impossible. If both are in P they occur consecutively between two increasing runs, each predecessor has three assignments and each successor two, giving 2^r3^l.

Let a=h−s, d=F−q and δ=d−a. If l=0 then δ=0. Move the root one place to l=1 with r=m−2≥1, holding its vertex name. Writing z=2^r gives d=3+4z, q=1+2z, a=1+3z and δ=2+z. The inherited physical commands R (keep root) and N (mark newborn) both add root→newborn→all other old vertices. For k≥1 their word R^k N R returns (b,t,t,b), where
```
t=2^k d+q−δ,   b=3t+2δ.
```
For reference, applying R once replaces s by q, F by 2F−q and leaves q fixed. After k steps h=a+q+(2^k−1)d, s=q, F=q+2^kd. Applying N then R gives the displayed formulas. The output is strong since each birth adjoins a vertex reachable from the old root and reaching the other old vertices.

The return is primitive precisely when gcd(t,2δ)=1. The unprepared phase has b=3t with t>1 for every k≥1, so no primitive return. In the prepared phase t is odd. If k is a positive multiple of the multiplicative order of 2 modulo the odd part of δ (period 1 if that odd part is 1), then t≡d+q=4+6z≡−8 modulo that odd part. Thus there are infinitely many primitive branches at each exact initial target. For odd m, r is odd and z≡−1 mod3, so taking k even gives t≡0 mod3 and b≠0 mod3. These outputs meet the [native balanced controller](../source_selection_closure_20261103/NATIVE_SELECTION.md) entry guard. This proves an infinite family of physical same-target preparations with renewable exits, not universal same-target normalization. The original H17 preparation moves two positions and gets a stronger all-k branch theorem.

## Design interpretation and proof targets

A real relaxation of the repair equation b−3t=2δ is nearly vacuous; integer feasibility adds parity and coprimality, and a physical lift still requires (h,s,q,F) to come from an actual tournament with a marked root. The distinction mirrors [the baby Hodge wall](../../../01-canon/theorems/THM-509-the-two-layer-baby-hodge-dichotomy.md): a coarser observer and an integral realizability question should not be conflated. This is not a statement about algebraic cycles in the Hodge conjecture. The rational affine formula is exact, but its extension to real points is not constructive. Parity plus the odd-part period supplies actual integer parameters without an appeal to an ambient continuous approximation.

See [finite legal repair geometry](LATTICE.md) for a precise conductor theorem and connections. Chan–Pak remains open: missing is an exact-source supply for targets outside the constructed family, or a target-decreasing legal repair that returns with its required interface.
