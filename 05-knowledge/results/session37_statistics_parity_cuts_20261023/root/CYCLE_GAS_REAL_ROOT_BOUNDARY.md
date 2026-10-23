# Recovering the nine-vertex cycle-gas boundary with a transparent strong witness

Status: PROVED elementary boundary theorem with an explicit strong witness. This concerns the cycle-fugacity polynomial below, not every possible fermionic representation of tournament counting.

## 1. Inheritance and the precise partition function

The threshold itself is OLD repository work: THM-020 `01-canon/theorems/THM-020-real-roots.md` gives real-rootedness through order8, and THM-025 `01-canon/theorems/THM-025-real-roots-disproved.md` already supplies an order9 counterexample with polynomial1+94z+10z²+z³. Session37 does not claim a new threshold. It supplies the particularly simple three-triangle/one-reversed-arc witness1+49z+7z²+z³, an elementary degree-two/Mantel proof of the lower boundary, and the precisely scoped positive-free-fermion obstruction.

For a tournament T, let O(T) be its directed simple odd cycles, each counted modulo cyclic starting point, not modulo reversal. Its conflict graph joins cycles sharing a vertex. Define

I_T(z)=sum_(vertex-disjoint odd-cycle families F) z^|F|.

The empty family contributes1. The inherited odd-cycle formula gives H(T)=I_T(2). Sources: THM-499 `01-canon/theorems/THM-499-the-spectral-reframe-boundary-H-is-finer-than-the-spectrum.md`, THM-505 and the primary Irving–Omar/Grinberg–Stanley path identities. The present result does not claim this representation as new.

This is a hard-core gas of cycles. Its exclusion is imposed on shared vertex support. It is not automatically a free fermion gas, a physical quantum Hamiltonian, or a determinant of a positive one-particle operator.

A finite positive free-fermion occupation polynomial has the form

F(z)=det(I+zK)=product_j(1+z lambda_j), lambda_j>=0,

where K is Hermitian positive semidefinite. Every zero is negative real. This necessary property supplies a sharp test.

## 2. All tournaments through order eight pass the root test

THEOREM. If |T|<=8, every zero of I_T(z) is negative real. If I_T is constant there are no zeros.

Proof. Each cycle uses at least three vertices, so I_T has degree at most2. Write I_T=1+mz+ez². Here m is the number of odd cycles, and e is the number of unordered disjoint cycle pairs. Construct the compatibility graph on the m odd cycles, with an edge for disjointness. Three pairwise disjoint cycles would require at least9 vertices, so this graph is triangle-free.

For each compatibility edge uv, triangle-freeness gives d(u)+d(v)<=m. Summing gives sum_v d(v)²<=me. Cauchy–Schwarz gives sum_v d(v)²>=4e²/m when m>0. Hence e<=m²/4. Thus the quadratic discriminant m²−4e is nonnegative. If e>0 the sum of its two roots is−m/e<0 and their product is1/e>0, so both are negative. If e=0<m the lone root is−1/m. The empty case is immediate. QED.

This is a root criterion, not an intrinsic one-particle realization. Algebraically the roots permit a diagonal positive realization when degree<=2, but it need not intertwine any tournament operation or assign one mode per vertex.

## 3. A minimal-order strong counterexample

Start with three directed triangles A,B,C, with all cross-arcs A->B->C and A->C. Choose a in A and c in C, and reverse just a->c into c->a. Call the tournament U. It is strong: all three internally strong blocks lie on the directed block cycle A->B->C->A.

Every cycle not internal to a block must use c->a, since this is the only backward cross-arc. Hence all new cycles intersect one another.

For a new cycle avoiding B, its A segment starts at a and its C segment ends at c. Each segment has a unique possible directed path of any length1,2,3 in its triangle. The total cycle length is odd precisely for the four pairs of segment lengths(1,2),(2,1),(2,3),(3,2). Each of these four new odd cycles is disjoint from the internal B triangle, and from neither other internal triangle.

For a new cycle meeting B, choose segment lengths r,s,t in{1,2,3}. The A and C segments are unique; a B segment of each length has three choices. Exactly14 of the27 triples have odd sum, so there are42 such odd cycles. Each meets every original triangle and hence is incompatible with all three.

Therefore

I_U(z)=(1+z)^3+4z(1+z)+42z
      =1+49z+7z²+z³.

Its cubic discriminant is−348172<0. Consequently it has one real zero and a nonreal conjugate pair. Positive coefficients imply that the real zero is negative. In particular H(U)=I_U(2)=135.

Together with §2, order9 is the exact smallest tournament order at which this odd-cycle partition polynomial fails real-rootedness.

## 4. What the obstruction means

No Hermitian positive semidefinite finite K can satisfy I_U(z)=det(I+zK). More generally a fixed positive free-mode product cannot represent every tournament odd-cycle gas. The witness exposes overlapping-cycle interactions that do not reduce to independent occupations.

This does NOT rule out an interacting fermionic model, a signed Grassmann integral, a non-Hermitian companion-matrix determinant, or a determinant ratio in the vertex-labelled square-zero algebra. Those constructions have different hypotheses and retain additional data. At z=2 it also does not preclude representing the single scalar135 as a determinant: the theorem concerns the entire fugacity polynomial.

## 5. Exact certificate

`probe_cycle_gas.py` independently enumerates all directed odd cycles of U by fixing the minimum cycle vertex and permuting the remainder, then computes disjoint-family coefficients by a lowest-vertex include/exclude recurrence. The explicit out-neighbour masks are

[250,508,505,464,480,456,128,256,65].

The individual odd-cycle counts are c3=8,c5=20,c7=18,c9=3, totaling49. Direct Hamiltonian subset DP independently checks135. Normal and optimized Python use explicit exception checks and produce identical JSON. The all-order minimality proof is the triangle-free compatibility argument, not a finite enumeration extrapolation.
