# Independent audit: finite simplicial gauge and code-spectrum bridge

Verdict: PASS for the stated finite claims. No substantive mathematical correction required. This is not an audit of any attached manuscript, and does not establish continuum Yang–Mills theory, a physical gauge identification, or a uniform continuum mass gap.

## Scope and independence

Reviewed `session34-gauge/SIMPLICIAL_WILSON_FANO.md`, the full rainbow histogram, and `session34-root/FINITE_CODE_SPECTRAL_BRIDGE.md`. Constructed a separate exact checker rather than rerunning the producer checks. `independent_audit.py` and its normal/optimized JSON certificates are in this directory. Python normal and `-O` outputs agree; failures use explicit exceptions.

## Independent spectral derivation

Let r_0=0 and r_s be the single bit indexed by nonzero syndrome s. Enumerate the 16 code characters by their values directly on C. For each character chi, the real-space syndrome quotient has adjacency

    A_chi(s,s+a) = chi(r_s + e_a + r_(s+a)),  a != 0.

The parenthesized expression has syndrome zero. Each A is symmetric, has zero diagonal and seven signed off-diagonal entries per row. The sector Hamiltonian is

    t(7I-A_chi) + diag(0,1,1,1,1,1,1,1),  t=h/J.

Exact 8×8 determinants reproduce four types with multiplicities 1,7,7,1. This construction is independent of the producer's Fourier-coset/weight-enumerator construction. Full characteristic polynomials, including repeated-weight factors, are recorded in the certificate.

The cubic coefficients can also be obtained without the characteristic polynomial. Write A in blocks (0,v^T;v,Q). The low eigenvalue obeys the Schur-complement equation

    z = 7t - t² v^T [I+t(7I-Q)-zI]^-1 v.

Because v^T v=7 and z=7t+O(t²), expansion gives

    z = 7t - 7t² - (v^T Qv)t³ + O(t⁴).

Independent exact evaluation gives cubic coefficients −42,−6,+6,+42 with multiplicities 1,7,7,1. Substitution into the independently obtained determinants verifies these terms.

Consequences:
- There is precisely one low branch in each sector near t=0; the other seven approach 1. The four distinct cubic coefficients rigorously order low branches for sufficiently small positive t.
- The first excitation has degeneracy seven in that regime. It is not an invariant-sector excitation.
- Restoring units multiplies dimensionless eigenvalues by J, hence Δ_full=36h³/J²+O(h⁴/J³). The dimensions are energy, as required when h and J both have energy units.
- The invariant 2×2 determinant is z²−(1+8t)z+7t. Its exact gap is sqrt((J−6h)²+28h²), with six other levels at J+8h. At J=0 its two upper level descriptions merge into the sevenfold eigenvalue 8h; no extra distinct level should be counted.
- For J>0,h>0, the ground eigenvalue is positive: equality in the two nonnegative terms would require both a constant cube function and support only on C. Perron–Frobenius proves uniqueness and code invariance for finite parameters. At J=0 the ground energy is zero and unique. At h=0 the zero eigenspace has dimension16, which is excluded by the hypothesis h>0.
- The assertion that the full ground lies in the invariant sector holds for every finite J>=0, not only asymptotically. The asserted first excited sector ordering is only asymptotic, correctly so.

## Independent simplicial and Fano checks

The symmetric link convention u_ji=u_ij is essential and correctly stated. Using the skew tournament matrix instead changes directed traversal signs and invalidates unadjusted Wilson formulas. Vertex switching acts on links by a coboundary, not by line-code translation.

For all 1024 K5 edge assignments the independent checker verifies every tetrahedral Bianchi relation, root-gauge reconstruction, and unchanged triangle curvature; it finds exactly64 curvature assignments. This agrees with the general proof: root triangles parametrize the quotient, and tetrahedral equations supply all remaining triangles.

Independent boundary enumeration verifies the K5 normalized triangular partition polynomial

    1+5t⁴+10t⁶,

and the fixed C4 numerator

    2t²+5t⁴+8t⁶+t⁸.

The proof that the first nonempty closed surface has four faces and is a tetrahedral boundary is valid. The nonabelian transported identity also checks algebraically: both ordered sides reduce to U_ij U_jk U_kl U_li. Neither transport nor factor ordering may be omitted.

Reconstructing all105 undirected K7 four-cycles yields21 rainbow cycles, grouped into seven supports of multiplicity3. The supports match the specified line normals (4,1,3,7,6,5,2). Direct Hamiltonian-path dynamic programming, independent of the claimed Fourier formula, verifies on all128 Fano orientations

    4H = 693 + 3 Σ_C Wbar_C.

The histogram is H=189 on16 and H=171 on112. The normalized Wilson sum is21 on16 and −3 on112. Switching vertex0 preserves every C4 Wilson loop but changes the baseline count189 to135. Therefore the invariant extension on the saturated slice is not the original tournament Hamiltonian count away from its chosen representative slice.

A separate binary-rank calculation gives rank(L+B)=13. Thus the128 slice points lie in distinct64-element gauge orbits and the saturated field count8192 is correct. The partition-function factors64 and16/112 are therefore correct for unnormalized counting measure.

## Independent full-rainbow extension check

Rather than enumerating32768 root-gauge representatives, the auditor enumerated the image of the linear rainbow-loop evaluation map. This image has2048 elements, hence rank11 and kernel dimension10. The cut vectors together with line-code vectors span dimension10 and annihilate all21 rainbow loops, proving the kernel is exactly B plus the embedded four-dimensional line code.

The full 2²¹-field histogram agrees exactly:

    R=-7:282624; -3:637952; 1:688128; 5:322560;
      9:143360; 13:21504; 21:1024.

The maximum21 requires every normalized loop to equal1. Its1024 fields are precisely the gauge saturation of the16 Fano code orientations. This conclusion concerns the selected normalized rainbow action, not the unconstrained triangular Wilson action.

## Interpretation guardrails

The notes correctly distinguish the full edge-switching group, Fano line flips, and Hamming-code translations. Projection onto code invariants is an exact algebraic operation, but only an independent physical argument could identify that sector as a physical gauge-invariant Hilbert space. A gap diverging in that sector and a full gap tending to zero are compatible.

The limit J/h→infinity is a finite128-state strong-penalty limit. It is neither a thermodynamic limit nor a continuum lattice-spacing limit. Four edge factors, four restricted line-bit factors, four surface-face factors, third-order tunneling and inverse-square penalty scaling refer to different gradings. No renormalization, reconstruction, reflection positivity, continuum existence, or physical mass-gap conclusion follows from these finite certificates alone.
