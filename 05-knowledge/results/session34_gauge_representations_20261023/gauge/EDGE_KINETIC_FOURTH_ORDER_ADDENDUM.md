# Actual vertex-gauge quotient: matching tunneling and fourth-order cancellation

Status: PROVED finite-dimensional small-parameter expansion, with exact combinatorial coefficients. This is a32768-state binary lattice model. Its spectral gap is NOT a continuum Yang–Mills mass gap. It is also NOT asserted to be the gap on the unrestricted2^21-dimensional edge Hilbert space, which retains additional vertex-gauge sectors. The root's128-state line-flip model is a different operator.

Use relative edge fields x=u/u^0, and let X_e flip the single edge e. Work in the Hilbert space of functions invariant under vertex-cut gauge transformations B. Its32768 basis states are gauge orbits; X_e acts as translation by [e]. The21 translation vectors are distinct. Define

H_edge(h,J)=h sum_(21 edges)(I-X_e)+J v,
v(x)=(21-R(x))/8,  J>0, h>=0.                               (K1)

R is the normalized rainbow Wilson sum. By the full census its distinct v-values are0,1,3/2,2,5/2,3,7/2; the16 ground orbits have v=0 and the next unperturbed energy is J. Half-integer costs MUST NOT be rounded. All21 elementary edge flips from any ground state have cost J, since each edge lies on four rainbow quadrilaterals.

Let A_crown be adjacency of the16-ground-sector crown graph K_(8,8) minus a perfect matching, whose seven generators are the matching representatives proved in LOGICAL_EDGE_DISTANCE_ADDENDUM.md. Its spectrum is7 (once),1 (seven),-1 (seven),-7 (once).

**Theorem (through fourth order).** The sixteen low eigenvalues of(K1), counted with multiplicity, have expansions

E_mu(h,J)=21h-21h^2/J-6 mu h^3/J^2-77h^4/J^3
          +O(h^5/J^4),                                     (K2)

where mu ranges over that crown spectrum. Consequently the lowest spectral gap is

Delta_edge(h,J)=36h^3/J^2+O(h^5/J^4)                         (K3)

as h/J->0+, with no fourth-order term. The error constants concern this fixed finite model. They assert no volume-uniform, continuum, or thermodynamic limit.

## Schur-complement derivation

Remove the common21h shift and divide by J. Put t=h/J, A=sum_e X_e, D=v, and P the16-dimensional zero space of D; Q=I-P. We study D-tA. The cut quotient is an abelian translation space. Ground translations K/B leave D invariant, so every diagonal effective matrix entry is identical across ground sectors.

PAP=0: no one-edge subset belongs to the kernel modulo cuts. For a low eigenvalue lambda, eliminating Q gives the exact equation on its ground component

lambda psi = -t^2 P A Q (D_Q-t Q A Q-lambda)^(-1) Q A P psi.  (K4)

Since D_Q>=I and ||A||<=21, the inverse has a convergent Neumann expansion whenever21|t|+|lambda|<1. Standard finite-matrix spectral separation, or direct Weyl bounds, isolates these16 levels for sufficiently small t (for example|t|<1/100 suffices for separation). Equation(K4) also gives lambda=O(t^2).

Expand the inverse, using lambda=O(t^2):

lambda psi = [-t^2 B0-t^3 B1-t^4 B2-lambda t^2 B_lambda
              +O(t^5)] psi,                                (K5)

B0=P A Q D_Q^-1 Q A P,
B1=P A Q D_Q^-1 Q A Q D_Q^-1 Q A P,
B2=P A Q D_Q^-1 Q A Q D_Q^-1 Q A Q D_Q^-1 Q A P,
B_lambda=P A Q D_Q^-2 Q A P.

All remainder bounds are operator-norm bounds on this fixed finite space, locally uniform in the low spectral window. The required matrices are:

B0=B_lambda=21I; B1=6 A_crown; B2=518I.                      (K6)

Proof of the coefficients:
- Two flips return to a ground orbit only by repeating the same edge. There are21 choices, each intermediate cost1, giving B0 and B_lambda.
- Three flips reach a different ground orbit only via one of the seven unique3-edge matching representatives. Each matching has6 orderings. The first intermediate edge has cost1; the two-edge partial state has the same rainbow syndrome as the remaining third edge because the full matching lies in the kernel, so its cost is also1. This gives B1. No three-flip sequence closes to the same gauge orbit because a nonzero cut has at least6 edges.
- At four flips there is no off-diagonal matrix element between ground orbits: a final even edge subset of size at most4 in the kernel must be empty, since the kernel weight enumerator has no2 or4. Thus B2 is diagonal, hence scalar by ground-translation symmetry. For any ordered pair of distinct edges(a,b), exactly two allowed closed four-flip words survive Q at every intermediate step: abab and abba. Their denominator is v(a+b), since first and third partial states are single-edge syndromes. Among420 ordered distinct pairs,42 have cost1,168 have cost3/2, and210 have cost2. Therefore

B2=2[42+168/(3/2)+210/2] I=518I.

The pair-cost counts are exact finite incidence data checked by the companion. They can also be read from common-rainbow-face counts k: v(a+b)=2-k/2. Here k=2 for42 pairs, k=1 for168, k=0 for210; the total common-face count252 equals21*4*3, as a consistency check.

Substitute lambda=-21t^2+O(t^3) in the last term of(K5). The fourth-order scalar becomes -518+441=-77. Thus the effective equation is

lambda psi=[-21t^2 I-6t^3 A_crown-77t^4 I+O(t^5)]psi.

Diagonalizing the explicit crown matrix yields(K2). Multiplicity counting can be obtained by the isolated finite spectral projection and the same Schur complement; no infinite-dimensional spectral theorem is needed. For t>0 the crown eigenvalue7 gives the lowest branch, eigenvalue1 the next seven, separated by36t^3+O(t^5). The other excited levels remain order J away. This proves(K3).

The finite Hamiltonian has a unique ground state for h>0: its off-diagonal entries are nonpositive, and the21 edge translations generate the entire quotient, so the connected finite Perron–Frobenius argument applies. This uniqueness is consistent with the small-t lifting of the16 classical ground sectors.

## Comparison with the line-flip model

The leading cubic matching tunneling has the same crown adjacency and coefficient6 as the root's line-code model, but this alone does not identify the models. Here the single flips are actual edges, the Hilbert dimension is32768 after ordinary gauge projection, the first two scalar shifts are21h and-21h^2/J, and all quartic splitting cancels. In the separately audited root model, line flips yield different fourth-order band corrections. Any comparison should cite that model's own derivation, not infer its coefficients from this note.

The term 'tunneling' here refers only to a finite degenerate perturbation expansion. Gauge projection is the six-dimensional vertex-cut action B, NOT quotienting away the four-dimensional code symmetry C. Further quotienting by C would remove the logical sectors responsible for this small gap, and is not the Hilbert space in(K1).

## Reproduction

edge_perturbation.py uses exact Fraction costs, verifies32768 rooted-gauge states and16 ground orbits, all21 single-edge costs, all441 two-step words, all9261 three-step words, and all420 distinct-pair costs. Normal and optimized outputs agree. The quartic coefficient is obtained from the explicit path classification and exact rational count, not a numerical fit to eigenvalues or a finite-difference extrapolation.
