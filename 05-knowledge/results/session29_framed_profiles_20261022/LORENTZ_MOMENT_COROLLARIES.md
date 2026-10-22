# A genuine indefinite observer in tournament wiring

These are exact consequences of the independently audited three-port transfer in `session29-profiles/ADDENDUM.md`. They concern one infinite source family, with multiplicity on edges, not all tournament classes. The combinatorial proof is in that addendum; this note isolates its spectral consequences.

Let

K = [[172,-24,408], [24,-4,68], [408,-68,706]],
J = diag(1,-1,1), e = (1,0,0)^T,
g_r = e^T K^r e.

For r>=2, g_r is the sum of the Hamiltonian counts of strong target tournaments over the family's 18^r outward metagraph edges. Parallel lines are counted separately. The extension g_0=1 and g_1=172 is convenient but not part of that outward-edge assertion.

## 1. The scalar observable requires one negative direction

Direct multiplication gives JK=K^T J. The Krylov matrix B=[e,Ke,K^2e] has determinant -4406400, and det K=50400. Therefore for every integer s>=0,

H_s := (g_(s+i+j))_(0<=i,j<=2) = B^T J K^s B,

det H_s = -19416360960000 * 50400^s < 0.

All g_s are positive. Hence each H_s has signature (2,1): its determinant is negative, so it has an odd number of negative eigenvalues, while its positive top-left entry rules out negative definiteness. None is singular.

This proves a sharp, observation-specific obstruction. No tail (g_(s+r))_(r>=0) is the moment sequence of a positive measure on the real line, nor can it equal <u,A^r u> for a self-adjoint Hilbert-space operator with the requisite moments defined. In either proposed realization H_s would be a Gram matrix and positive semidefinite. This is impossible. Any finite linear realization needs at least three states because H_0 has rank three. The displayed three-state indefinite realization achieves dimension three and exactly one negative direction, so both bounds are sharp.

An explicit negative square witness is

q(x)=x^2-874x-45144,
q(K)e=10800*(0,1,0)^T.

For the linear functional L(x^r)=g_r,

L(q^2)=(q(K)e)^T J(q(K)e)=-116640000.

This is not a claim that K itself cannot be made self-adjoint in a positive inner product. It has three distinct real eigenvalues and can. The obstruction concerns this scalar count observable with the same vector on both sides, not the bare eigenvalue set. Real poles alone do not certify a positive spectral measure.

## 2. Exact rational resolvent and sharp growth

The full scalar series is

sum_(r>=0) g_r z^r
 = (1-702z+1800z^2)/(1-874z-43344z^2-50400z^3).

Its denominator has degree three with no cancellation, also certified by the nonsingular three-state Hankel matrix. The characteristic polynomial is

p(x)=x^3-874x^2-43344x-50400.

Its roots lie respectively in (-46,-45), (-2,-1), and

(921.115395614954,921.115395614956).

The last root rho is strictly dominant in modulus. If the other two roots are alpha and beta, the exact three-atom signed spectral expansion is

g_r = sum_(lambda in {alpha,beta,rho})
      [(lambda^2-702lambda+1800)/p'(lambda)] * lambda^r.

The two outer roots have positive weights and the middle root has negative weight: the numerator is positive at all three roots, and p' has signs +,-,+ there. This is the spectral form of the same one-negative-direction obstruction.

In particular,

g_r ~ c rho^r,
c=(rho^2-702rho+1800)/p'(rho)>0.

The average Hamiltonian count at the strong endpoint of a uniformly chosen outward edge is therefore

g_r/18^r ~ c (rho/18)^r.

This is an algebraic asymptotic for an explicitly specified edge-weighted family. It is not an asymptotic for A000568, for the maximum Hamiltonian count M_k, or for distinct target classes. The exact checker records rational enclosing intervals for rho, c, and rho/18.

## 3. What the Lorentzian analogy preserves

For each marked tile, mirror changes its transfer matrix to J K_t^T J. An even mirror-fixed seam word has target count H^2-S^2+R^2, where the prefix's H,S,R are actual path and path-cover counts. This preserves composition, reflection, and the target Hamiltonian observable exactly. It is stronger than noticing a repeated constant or dimension.

The negative term comes from a forbidden port connection in the path-counting problem, not from a physical timelike direction. No embedding into II_(25,1), no identification with the Leech cusp, and no map to a polyhedral Dehn invariant has been constructed. The useful shared principle is that composition requires marked interface data: scalar counts, unframed angular spectra, and unmarked incidence decks forget information needed when pieces are joined.

## Certificate

`check_lorentz_moments.py` uses only Python's standard library and exact integers/rationals. It verifies the adjoint relation, determinants, negative square, recurrence, 40 shifted Hankel determinants, and three disjoint root brackets. Normal and optimized-mode outputs are identical. The determinant identity and spectral proof establish the infinite statements; finite checks are reproducible controls, not the induction premise.
