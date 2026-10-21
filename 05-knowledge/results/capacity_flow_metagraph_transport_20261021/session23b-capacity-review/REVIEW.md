# Independent review: capacity reorientation and the fixed-path face

Date: 2026-10-21. Verdict: **PASS**, with classical imports and finite controls separately identified. The audited packet is `../session23b-capacity/PROOF.md` and its accompanying program. This review and its independently written tests live outside the repository; no Git operation was performed.

## 1. Scope and source inheritance

I read the repository predecessor `01-canon/theorems/THM-646-line-metagraph-score-complement-skeleton.md`. Its part (i) already proves the pointwise outdegree law for the nonpath-complement operation, with base path directed from its last label toward its first. The present first-to-last convention reverses the endpoint constants, yielding `s+s'=(n-1)1-p` with `p=e_last-e_first`. Therefore the score law is inherited, not newly established here. The new packet correctly identifies its added content as the affine integral capacity/demand transport and the resulting bounded-flow enumerators. Neither this audit nor the packet establishes literature priority for general reorientation transformations.

Primary-source spot checks confirmed the relevant lattice-polytope imports: Stanley, *Combinatorial Reciprocity Theorems*, Section 6 (https://math.mit.edu/~rstan/pubs/pubfiles/22.pdf), and Beck–Sottile, *Irrational proofs for three theorems of Stanley*, Theorem 2 and Corollary 4 (https://arxiv.org/html/math/0501359v4). The packet also supplies an elementary integrality proof, so its special-case integrality does not depend on merely citing total unimodularity.

## 2. Exact transport and signs

With incidence column `e_head-e_tail`, reversing an edge and replacing `f` by `M-f` changes its divergence contribution from `B_e f_e` to `B_e f_e-M B_e`. Hence

    B' Phi_D(f) = B f - M B_D 1.

The map is an ambient lattice affine isometry and is its own inverse after changing the orientation. For tournament outdegree vectors, `B_D1=s'-s`, so the invariant is `Ms+b`. Equal labelled score vectors imply `B_D1=0` and therefore a true bounded-circulation isomorphism.

For the actual operation `D=E\P`, `B_D1=B1-p`. A regular source has `B1=0`, giving target demand `Mp`, not zero. The first vertex is the source and the last is the sink under the inflow-minus-outflow convention. The path itself, its orientation, and both labelled endpoint ports remain fixed. The stated triangle and five-vertex counterexamples to circulation-to-circulation preservation are appropriate.

The multivariate weight formula and total-use correction are also correct. No total-use, tournament-isomorphism, or automorphism-action invariance follows from the scalar lattice equivalence.

## 3. Full regular polytopes

For odd `n>=3`, `1/2` on every arc is a strictly feasible circulation, so the affine hull is the incidence kernel of dimension `d=m-n+1`. The vertex integrality argument is valid: nonintegral edges cannot have a degree-one vertex, so a nonintegral support would contain an undirected cycle permitting a two-sided signed-circulation perturbation.

Interior integer points at dilation `M>=2` are exactly `1+` the points at dilation `M-2`. Thus reciprocity gives `N(-M)=(-1)^d N(M-2)`, and series reciprocity gives a palindromic numerator of exact degree `d-1`. The n=1 exception is necessary and explicitly retained. The same-score-fibre conclusion is not a binary-code assertion: a regular five-vertex fibre has size24.

The exact five-vertex polynomial and numerator are correct:

    N5(M)=(M+1)^2(3(M+1)^4+5(M+1)^2+4)/12,
    h5*(z)=1+17z+72z^2+72z^3+17z^4+z^5.

The proven degree-six bound makes seven exact values a proof by finite determination, not an extrapolation. Four additional capacities independently reproduced the formula.

## 4. Fixed-path face and its presentation interpretation

If two regular bases contain the same labelled directed Hamiltonian path, their disagreement set avoids that path. Their affine map therefore preserves every equation `f_P=0`. At capacity1, reversing the support bijects such circulations with all regular tournaments containing the fixed path. Summing `H(T)/|Aut(T)|` is justified because the automorphism action on ordered Hamiltonian words is free.

For any proper nonempty vertex set of size `a`, a regular tournament has exactly `a(n-a)/2` outgoing cut edges. The path removes at most `min(a,n-a)`. For odd `n>=5`, the difference is positive. This proves `T\P` is strong, not merely weakly connected. Summing directed cycles through its arcs yields a strictly positive circulation and proves exact dimension

    delta=m-2n+2=(n-1)(n-4)/2.

The remaining incidence matrix still has rank `n-1`. Integrality follows either by the face argument or the elementary incidence proof. The denominator exponent is consequently exactly `delta+1`.

For the standard cyclic five-vertex base, write the nonpath arcs as `02,30,40,13,41,24`. Their circulation coordinates have the form

    (x,y,z,y,y,x),  x=y+z.

The capacity conditions reduce to `y,z>=0` and `y+z<=M`. This is a lattice-equivalent standard two-simplex and proves

    R5(M)=binomial(M+2,2),   sum R5(M)z^M=1/(1-z)^3.

At n=3 the remaining graph has no cycle and `R3(M)=1`; the dimension and strong-deletion assertions correctly exclude it.

## 5. Additional refinement discovered during this audit

For every odd `n>=5`, the fixed-path face has **exact codegree3**. There is no interior integer point at dilation1. At dilation2, every remaining coordinate would have to be1, but `B_D1=-p` is nonzero. Since `T\P` is strong, choose a directed first-to-last path `Q` in it. Then `1_D+1_Q` has zero divergence and all coordinates1 or2, giving an interior point at dilation3.

Therefore the numerator in denominator `(1-z)^(delta+1)` has exact degree `delta-2`. This does not imply palindromicity. Indeed, the independent program finds precisely34 binary flows on the nonpath graph of cyclic R7 with demand `p`; translation by `1_D` bijects these with the interior points at dilation3. Hence the restricted n7 numerator has constant coefficient1 and degree-seven leading coefficient34: it is **not palindromic**. The corresponding n5 count is1, consistent with its unimodular simplex.

This codegree proof was sent to the producer before it was incorporated into the packet. The producer independently verified the numerical nonpalindromicity witness and incorporated it into the final frozen proof. A separate direct orientation census in supplement.py also confirms its endpoint-shift score-fibre interpretation:1 presentation at n5 and34 at n7.

## 6. Centered-current interpretation

For `j=f-(M/2)1`, reorientation flips exactly the coordinates indexed by D. Thus `B'j'=Bj`, and each fixed-edge quadratic energy `sum r_e j_e^2` is preserved. This is an orientation-coordinate equivalence. It does not impose Ohm's law, vertex-potential representability, or energy minimization. The packet now states these limitations explicitly. Uncentered `sum r_e f_e^2` and total flow use are not claimed invariant.

## 7. Independent computational evidence

`audit.py` imports no producer module. Its full-flow route assigns nonpath coordinates and solves the Hamiltonian-path spanning-tree coordinates by prefix-cut equations, unlike the producer's star-tree elimination and score-coefficient recurrence. All checks use explicit exceptions, not removable assertions.

Verified in normal mode and Python `-O` with identical deterministic JSON:

- Full regular n5 counts at M=0,...,10: 1,24,219,1136,4175,12216,30429,67264,135621,254200,449031.
- Exact degree-six finite differences and the claimed numerator.
- Direct enumeration of path-normalized regular tournaments gives n3/n5/n7 counts1/3/91.
- All three regular n5 bases have the same restricted point sets under transport at M=0,...,6, with counts1,3,6,10,15,21,28.
- Restricted nonpath-complement transport agrees as complete point sets at n3,n5 capacities1,2,3 and n7 capacity1, retaining zero path coordinates and changed demand.
- 11,556 cut checks across all3 n5 and91 n7 regular path presentations; every deleted-path graph has positive outgoing cut size. Minimum remaining cut size is1 at n5 and2 at n7.
- All8,281 images of91 restricted binary flows across91 n7 regular bases have zero divergence and unchanged zero path coordinates.
- Exact interior-at-dilation3 counts:1 at n5,34 at n7.
- 1,536 unequal-capacity transport cases on K3 with capacities(1,2,3), arbitrary induced demands, centered-current divergence, and nonconstant fixed-edge energy weights.
- Independent permutation counts give H=175 and189 and automorphism orders7 and21 for cyclic and Paley R7; both have2640 unrestricted binary circulations, while their q^4 total-use coefficients are28 and21. This confirms genuine information loss.

The final frozen additions were separately checked in supplement.py, also in normal and optimized modes with identical output. All three n5 path-face total-use polynomials are exactly1+q^3+q^4,1+q^3+q^5,1+q^4+q^5. The cyclic representative's joint capacity/use series follows immediately by writing M=a+b+c and summing three geometric series, giving1/[(1-z)(1-zq^3)(1-zq^5)]. Independent Hamiltonian-word enumeration gives the n7 presentation histogram{171:57,175:25,189:9}; representative automorphism orders are3,7,21, and171/3+175/7+189/21=91.

No claim in this review classifies the full original metagraph, its port wiring, all tournaments, or A000568. The audit proves/checks the named transport and count invariants only.

## 8. Reproduction

    python3 audit.py > normal.json
    python3 -O audit.py > optimized.json
    cmp normal.json optimized.json
    python3 supplement.py > supplement-normal.json
    python3 -O supplement.py > supplement-optimized.json
    cmp supplement-normal.json supplement-optimized.json

The adjacent manifest records the independently reviewed producer snapshot and review artifact hashes.
