# Independent adversarial audit: cyclic physical incidence kernel

2026-10-25, Session 62. **PASS**, with the coefficient-field/weight scope caveat below. No repository files were edited by this auditor. All files here are audit scratch artifacts.

## Scope and exact geometric reduction

Let p>3 be prime and 1<=L<=p. The actual point set is

    T_L={(-a-3b,j-9b) mod p: a,b in {0,1}, 0<=j<L}.

Its four x-coordinates are distinct and its columns have no repetitions, so |T_L|=4L. For p=2,3 this point-indexing assertion fails; at p=3, w=(3,9) is zero and is not a direction. Thus those primes are outside the theorem, not untested cases. For L>=p the point set stabilizes at T_p; repeated j labels must not be counted as extra points.

The physical geometry is over F_p, but line-sum coefficients may lie in ANY field K. Put I={0,...,L-1} in Z/p, and use zero extension outside I. With column arrays A=f00,B=f01,C=f10,D=f11, set

    h=A, P=A+B, R=A+C, S=B+D,
    (tau_d r)_j=r_(j-d).

The horizontal line equations give S=-tau_9 R. Thus R=r must be supported in E=I intersect (I-9). The w-line label y-3x=j+3a gives

    (1-tau_3)P=-tau_3(1-tau_9)r.

The identity 1-tau_9=(1-tau_3)(1+tau_3+tau_6) implies

    P=-(tau_3+tau_6+tau_9)r+c*1.

No division by 3 is being used over K. Since p>3, translation by 3 is one p-cycle, so the kernel of 1-tau_3 consists exactly of constants in EVERY characteristic, including p.

The vertical sums are equivalent to

    sum h=0, sum r=0, p*c=0.

Requiring P=0 outside I is the only remaining support constraint. Reconstruction is unique:

    A=h, B=P-h, C=r-h, D=-tau_9*r-P+h.

This independently verifies the proposed cyclic reduction, including the sign and shift conventions.

## Independent proof of the p>=17 closed formula

For the wrap range, put m=p-L<=8 and rotate the missing interval G=I^c to {0,...,m-1}. Then E excludes G and G-9. For p>=17 these intervals are disjoint. At j in G, the term r_(j-9) vanishes. The remaining equation is

    r_(j-3)+r_(j-6)=c,

where the first term can survive only for max(0,m-6)<=j<=2, and the second only for max(0,m-3)<=j<=5, always also intersected with 0<=j<m.

The ranks of these boundary equations on r, with c fixed, are:

| m | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| boundary rank | 0 | 1 | 2 | 3 | 3 | 3 | 3 | 2 | 1 |
| c permitted by boundary | yes | yes | yes | yes | yes | yes | yes | no | no |

For m=7,8 an identically zero left-hand row forces c=0. For m<=6, every c is possible at the boundary. The sum-r equation is independent because at least one unconstrained r coordinate remains, EXCEPT p=17,m=8: E then has just one coordinate, already forced to zero. That is precisely the edge handled by the positive part below.

Together with the already proved nonwrap formula, for p>=17 and 1<=L<=p this proves

    dim J=(L-1)+(L-10)_+ +(L-p+6)_+ +(L-p+3)_+
          + 1_(char K=p)*1_(L>=p-6).

The modular constant mode enters at L=p-6. New geometric dimension terms enter at L=p-5 and L=p-2. At L=p the formula is 2p-2+1_(char K=p).

## Audit of the producer's all-prime path-forest theorem

For L<p, let Gamma have all vertices Z/p and edges {j-6,j-3} for j in G=I^c. Anchor all vertices outside E to zero. Gamma is a proper subgraph of the step-3 p-cycle, hence a disjoint union of paths and isolated vertices. Its equations are r_x+r_y=c on every edge.

Let f be the number of components with no anchors, o=1 if at least one such component has an odd number of vertices (otherwise 0), and b=1 if in every anchored component all anchors occupy the same bipartition class (otherwise 0). Isolated vertices count as components.

For c=0, each unanchored path contributes one alternating vector. Its total sum is zero for an even path and +/-1 for an odd path, over every field, including characteristic 2. Anchored components contribute nothing. Therefore imposing sum r=0 leaves dimension f-o.

A c=1 solution exists before the sum condition exactly when b=1. If b=1, choose an integer 0/1 solution by alternating on every path, with anchored vertices zero. Because L<p there is at least one edge, so its integer sum beta satisfies 1<=beta<=|E|<=p-1. If o=0 the sum cannot be changed by a homogeneous solution, and beta is nonzero in F_p, ruling out a c-mode. If o=1, subtract beta times a suitably signed odd free-path vector to make the sum zero. Consequently

    dim J=L-1+f-o+1_(char K=p)*b*o,   L<p.

The full-cycle L=p case is separate and equals 2p-2+1_(char K=p).

This proof passes the adversarial audit, including char 2/3, multiple anchors, isolated free/anchored vertices, tiny physical primes, and the possible sum obstruction. It is a theorem, not an extrapolation from finite tests.

The small-prime dimension tables below list the generic value for L=1,...,p; add one in characteristic p from the specified threshold onward:

| p | generic dim J for L=1,...,p | extra char-p mode from |
|---|---|---|
| 5 | 0,1,2,4,8 | L=4 |
| 7 | 0,1,2,3,4,8,12 | L=5 |
| 11 | 0,1,2,3,4,6,8,10,12,16,20 | L=9 |
| 13 | 0,1,2,3,4,5,6,8,10,12,16,20,24 | L=10 |

These are exceptional shift-overlap patterns, not collisions of the four physical columns.

## Independent finite-exact controls

`direct_audit.py` constructs the actual distinct point set, groups points into physical lines by determinants with each of the directions (1,0),(0,1),(3,9), and forms the FULL 0/1 line-incidence matrix. It imports no producer code and does not use the cyclic reduction to build this matrix.

- Physical primes: 5,7,11,13,17,19,23,29,31,37,43.
- Every integer 1<=L<=p, without quotienting, pruning, or sampling L.
- Coefficient characteristics: 2,3,5,7,11,13,17,19,23,29,31,37,43,47,101.
- **3,525 direct finite-field rank tests; zero mismatches.**
- Every direct rank agrees with the cyclic reduction. All p>=17 direct ranks agree with the displayed closed formula.

`exact_smith_audit.py` goes further: it computes exact Smith normal forms of all 235 FULL physical incidence matrices in the same (p,L) universe. This is not a Smith computation on the reduced equations. Every nonzero Smith invariant is 1 except for a single p exactly when the path-forest theorem predicts the extra characteristic-p mode. Thus these 235 finite geometries are independently verified over EVERY coefficient field, including characteristic zero. No hidden characteristic-2, characteristic-3, or other torsion occurs in this finite universe.

Results are in `direct_results.json`, `direct_cases.json`, and `full_smith_results.json`. Reproduce with:

    python direct_audit.py
    python exact_smith_audit.py

These controls are FINITE-EXACT evidence. The all-range conclusion rests on the proofs above, not on their extrapolation.

## Essential weighted-coupling caveat

The common-kernel theorem is valid over arbitrary K. The universal coupling identity dim ker M=2 dim J additionally needs THREE NONZERO weights in K.

For the physical coupling over F_p with p>3, weights (-3,-9,1) are indeed nonzero, so the coupled conclusion applies. It must not be silently extended with those literal numerical weights to characteristic 3: then the first two weights vanish.

A concrete hostile control is physical p=17,L=2 with coefficient field F_3. Here dim J=1 and the three pairwise-zero-line-sum spaces each have dimension 2. Literal weights (-3,-9,1)=(0,0,1) produce coupling nullity 4, not 2 dim J=2.

No claim is made about realizable binary errors, actual sections, or improvements of norm-graph failure bounds.

## Final consolidated-note readback and harness freeze

The consolidated producer note `cyclic_common_kernel.md` was independently reread through its arbitrary-index, complete pairwise-space, explicit-basis, and coupling conclusions:

- The cyclic reduction and path theorem hold for every nonempty subset I of F_p, not only intervals. For full I the graph has p isolated free vertices, so the same path formula gives the full-length answer.
- The explicit interval x-bases satisfy all listed boundary equations, have the listed independent dimensions, and become zero-sum after subtracting their sums at the low-block anchor. That anchor exists whenever it is used. The characteristic-p generator with r=e_(p-9)+e_(p-8)+e_(p-7)-3e_0 and c=1 is valid throughout m<=6.
- For n=|I|, eta_s=|I intersect (I+s)|, c_s=2 if eta_s=0 and 1 otherwise, and kappa the number of components on I under differences +/-3 or +/-9, the stated formulas are correct:

      dim W_u=2n+eta_3-4+c_3,
      dim W_v=eta_3+eta_9+kappa,
      dim W_w=2n+eta_9-4+c_9.

  The proof via the rank of an oriented bipartite graph incidence matrix works over every field. Each fixed-j u/w incidence subgraph is a four-edge square; different squares meet precisely through index differences +/-3 or +/-9. For p>=17 and interval I, the two wrap-9 joins occur at L=p-8 and L=p-7 and successively merge the three ordinary residue classes, giving the claimed kappa formula.
- Coupled nullity and rank statements are correct with the explicit nonzero-weight hypothesis; the characteristic-3 literal-weight boundary is correctly separated.

**Exact scope:** these conclusions have an independent proof readback. This auditor's finite matrix/SNF universe remains the 235 interval geometries stated above; it must not be relabelled as an exhaustive arbitrary-subset computation. The producer has a separate broader checker.

The two audit scripts now use explicit always-active `require` gates, with no bare Python asserts. Mathematical JSON and stdout contain no runtime measurements. Normal and `python -O` runs produced byte-identical three JSON certificates and byte-identical stdout for each harness. The deterministic comparison and hashes are recorded in `harness_receipt.json`.

Each harness has a `--negative-control` option that deliberately injects one false dimension claim. It must terminate nonzero with the intended explicit gate, both normally and under `-O`. These four rejection outcomes are also recorded in the receipt. Use a separate output directory for negative controls so they cannot replace successful certificates:

    python direct_audit.py --negative-control --output-dir /tmp/direct-negative
    python -O direct_audit.py --negative-control --output-dir /tmp/direct-negative-O
    python exact_smith_audit.py --negative-control --output-dir /tmp/smith-negative
    python -O exact_smith_audit.py --negative-control --output-dir /tmp/smith-negative-O

Dependencies: Python 3, NumPy, and SymPy. The exact Smith harness imports only `incidence` and `require` from this auditor's `direct_audit.py`, not producer code.
