# Independent audit: the Lorentzian tournament count observer

Verdict: ACCEPTED, with no corrections required to the reviewed report.

Reviewed: `LORENTZ_MOMENT_COROLLARIES.md`, `check_lorentz_moments.py`, and its exact output. The underlying all-word combinatorial transfer was independently accepted in `session29-wiring-audit/AUDIT.md`; this audit checks the new spectral and asymptotic consequences, not the entire combinatorial construction again.

## Proof audit

The identity JK=K^T J and Je=e imply

(g_(s+i+j))_(i,j=0..2) = B^T J K^s B,

where B=[e,Ke,K^2e]. Exact determinants are det B=-4406400 and det K=50400, giving

det H_s=-19416360960000*50400^s

for every s>=0. These matrices are real symmetric and nonsingular. Their negative determinant permits either one or three negative eigenvalues; their positive first diagonal entry rules out three. Thus every shift has inertia (2,1). Positivity of every g_s follows independently from g_0,g_1,g_2>0 and the positive-coefficient cubic recurrence, as well as from the audited counting interpretation.

A positive measure would make its moment Hankel matrices positive semidefinite. A same-vector self-adjoint Hilbert-space representation would make them Gram matrices, under the stated moment/domain assumptions. Therefore the no-positive-moment assertion is correct for every tail, not merely an artifact of the artificial empty-word term. Hankel rank three forces linear realization dimension at least three. J realizes the original sequence with exactly one negative direction. If desired, the same sharpness for each shifted sequence follows using the symmetric form JK^s, which is congruent to H_s through B and hence also has inertia (2,1).

The report correctly distinguishes this obstruction from whether K itself admits some positive inner product making it self-adjoint. Its three distinct real eigenvalues do permit such an inner product. The negative result concerns this specified scalar observable with the same vector on both sides.

The polynomial witness q(x)=x^2-874x-45144 satisfies q(K)e=(0,10800,0)^T and L(q^2)=-116640000 exactly.

## Resolvent, residues, and growth

Direct symbolic inversion gives

sum_(r>=0)g_r z^r=(1-702z+1800z^2)/(1-874z-43344z^2-50400z^3).

The numerator and denominator are coprime. The characteristic polynomial has one root in each of (-46,-45), (-2,-1), and (921.115395614954,921.115395614956), certified by exact endpoint signs. Since it is cubic, these are all its roots, distinct and real. The positive root is strictly dominant in absolute value.

For each root lambda the residue weight is

(lambda^2-702lambda+1800)/p'(lambda).

The numerator is positive at each root. The derivative signs, ordered by the roots, are +,-,+. Thus the signed spectral expansion has precisely one negative weight. The independent checker verifies its first 80 moments by taking polynomial traces in Q[x]/(p), without evaluating numerical eigenvalues or importing producer code.

The coefficient c of the dominant term is positive. Rational interval arithmetic encloses it between approximately0.228309188369572646 and0.228309188369577137; the exact rational endpoints are in `evidence.json`. Likewise rho/18 is between approximately51.1730775341641111 and51.1730775341642222. Decimal displays here are orientation aids; the artifact records exact rational bounds.

Therefore g_r~c*rho^r and the mean H at the target of a uniformly sampled outward edge is asymptotic to c*(rho/18)^r. Dividing by18^r is justified by the separately audited edge count, including parallel-edge multiplicity. This is not an asymptotic for unweighted target classes, A000568, or the maximum tournament Hamiltonian count.

## Independent controls

`independent.py` uses SymPy with exact integers and rational numbers. It does not import or invoke the producer implementation. It checks:

- characteristic polynomial, symbolic resolvent and numerator/denominator coprimality;
- adjoint relation, Krylov determinant and exact negative polynomial vector;
- 80 moments derived from the rational series, rather than repeated matrix multiplication;
- 76 shifted Hankel determinant and full matrix identities;
- 80 signed spectral moments through quotient-algebra traces;
- exact disjoint root brackets and valid rational enclosures for c.

Normal and optimized runs give byte-identical JSON. The original standard-library checker was also replayed and its output matches the submitted output exactly. These finite controls support the reviewed algebraic identities and proofs; the infinite statements do not rely on extrapolation from them.

## Preserved scope

The Lorentzian form describes an exact marked-interface transfer and mirror-adjoint relation. Neither this audit nor the report constructs a map to a Lorentzian lattice, the Leech cusp, a polyhedral Dehn invariant, or a physical time coordinate. That restraint is mathematically necessary and is correctly retained.
