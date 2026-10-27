# Session 77: typed color contraction, alignment memory, and Petersen

## Scope and conventions

These are finite-dimensional SU(3) representation calculations over C, not a theorem about QCD dynamics, confinement, continuum construction, or mass gap. V=C^3 carries the fundamental representation, V* its dual. The ordered orthonormal basis is indexed 0,1,2; epsilon_012=epsilon^012=+1. The invariant state epsilon^{ijk} belongs to V^tensor3; the dual volume form epsilon_ijk belongs to (V*)^tensor3. Their numerical arrays agree in this basis, but V and V* are not thereby identified as complex SU(3) modules. The B_A below are invariant states in V^tensor6; their Gram pairing is Hermitian, using adjoints for the dual contractions. Tensor coordinates here are ordinary commuting numbers, not quark-field Grassmann variables, so no additional fermionic sign rule is asserted. Named tensor ports retain their order. All contractions below use these fixed conventions. The elementary algebra and ETF/Petersen structures are classical; the report makes no priority claim.

## 1. Why an untyped ternary color tree cannot work

For zeta a primitive cube root of unity, zeta I acts on V^tensor3 as 1 and on V as zeta. Equivariance therefore forces every map V^tensor3 -> V to vanish. Epsilon is instead a map V^tensor3 -> C. Its scalar output cannot be fed into another fundamental port. Moreover, a rooted full ternary node has three child ports AND one parent port; an epsilon web vertex has three incident ports in total.

A valid recursively composable mixed type is V tensor V* tensor V -> V. Its intertwiner space has dimension two, with basis

F(a,phi,c)=phi(a)c,  G(a,phi,c)=phi(c)a.

Indeed this Hom space is isomorphic to End_SU3(V tensor V), and V tensor V=Sym^2 V plus Lambda^2 V is a sum of two inequivalent irreducibles. At a=e0, c=e1, phi=e0*, F=e1 but G=0. Triality alone forgets this channel distinction.

Important qualification: for the FULL unprojected space V^tensor r tensor (V*)^tensor s, triality r-s=0 mod3 is necessary AND sufficient for some singlet to exist. Sufficiency follows by pairing min(r,s) ports with deltas and grouping the remainder into epsilon triples. It does not determine the singlet dimension, the value of a chosen contraction, or whether a selected irreducible channel contains a singlet: the irreducibles 8 and 10 have zero triality and no invariant vector.

End(V)=1 plus 8 is a nine-dimensional composition-closed matrix carrier: (x phi)(y psi)=phi(y)x psi. Nine is minimal only among linear matrix carriers containing every rank-one operator. This does not replace the full tensor category or justify an absolute minimality claim.

## 2. An exact C3/A2 connection inside the two octets

The standard decomposition is V^tensor3=1 plus 8 plus 8 plus 10. In the highest-weight subspace with a=|001>, b=|010>, c=|100>, the octet multiplicity space is the plane of coefficient sum zero. Raising the single 1 to 0 sends all three basis vectors to |000>, so this plane is exactly the highest-weight kernel of weight (1,1).

Pair-12 symmetric/antisymmetric orthonormal vectors are
u=(2a-b-c)/sqrt6, v=(b-c)/sqrt2.
Pair-23 vectors are u'=(a+b-2c)/sqrt6, v'=(a-b)/sqrt2.
Their overlap matrix is

    [ 1/2       sqrt3/2 ]
    [ sqrt3/2  -1/2     ].

Thus re-pairing mixes the two octet channels; retaining only the label `8` loses observable data. Cyclic permutation of a,b,c acts on this centered A2 plane. The operator C-C^-1 has square -3 there. This provides an actual representation-theoretic bridge to the centered triangular/Eisenstein C3 operator, rather than an identification of epsilon with an additive scalar charge. It does not identify the surrounding tournament and SU(3) composition algebras.

## 3. Six named color ports give a five-dimensional signed carrier

For each pair i<j in {1,...,5}, put A=(0,i,j) and B=its sorted complement. On a color assignment c in {0,1,2}^6 define the integer tensor

B_A(c)=epsilon(c_A) epsilon(c_B).

There are ten such unordered baryon partitions, with the block containing 0 listed first. In lexicographic pair order, the exact Gram matrix is determined entrywise by this displayed contraction formula. Its diagonal is 36, every off-diagonal is +12 or -12, and

G^2=72G, rank(G)=5, spectrum(G)={72 multiplicity5, 0 multiplicity5}.

The normalized vectors B_A/6 therefore form a real (10,5) equiangular tight frame: absolute overlap 1/3 and squared overlap 1/9, frame bound 2. An independent infinitesimal invariance computation on all 90 zero-weight words gives a 360-by-90 system with kernel dimension 5. Hence these ten tensors span ALL six-fundamental invariants, not merely a five-dimensional subspace selected in advance.

A pointwise relation on all 729 color assignments is

-B_012|345 + B_013|245 - B_014|235 + B_015|234 = 0.

The orbit of this one relation under all 720 named-port permutations spans the entire five-dimensional coefficient kernel, as the checker verifies exactly. Thus these port-translated four-term relations present the six-port quotient. Every three-channel normalized Gram minor has determinant 1-3/9 +/- 2/27, namely 16/27 or 20/27, and every pair is independent. The displayed four-channel relation therefore has minimum possible dependent support: this dictionary has spark 4. This is a six-port statement, not an all-size lower bound.

Dropping minus signs destroys this relation. Declaring distinct baryon partitions orthogonal also fails: reclustering overlap has absolute value 1/3. The exact source-preserving norm of sum a_A B_A is a*G a, not a sum of independent coefficient norms. The Gram matrix is the gluing object; rank and absolute overlaps are invariant under independently changing tensor signs, while matrix entries and coordinates transform covariantly.

This five-dimensional color space is not by itself a prediction of any physical state's amplitudes or hidden-color percentage.

## 4. Petersen appears after an explicit sign gauge

Use the tensor ordering above and D=diag(1,1,1,1,1,1,-1,-1,-1,1). Set

S=D(G-36I)D/12,  A_P=(J-I-S)/2.

The checker verifies that A_P is a simple cubic adjacency matrix satisfying

A_P^2=2I-A_P+J,

with spectrum 3^1, 1^5, (-2)^4. It also finds and checks an explicit isomorphism to KG(5,2), the graph on the ten edges of K5 with adjacency when those edges are disjoint. Thus A_P is the Petersen graph. An explicit map from the ten ordered partitions to edges of K5 is (01,12,23,03,24,02,14,34,13,04), with K5 vertex labels 0,...,4.

This is a sign-switched Seidel/Gram realization, NOT a canonical gauge-invariant physical graph and NOT a Petersen-to-K5 minor contraction. The canonical information here is the signed Gram data up to port/sign conventions and its associated switching class. The link with K5 uses edges and disjointness. Classical equiangular-frame/Seidel graph correspondence is documented by Waldron below; we do not claim this construction is new.

## 5. Closed composition requires typed webs, not trees alone

With the stated conventions the exact local identities are

sum_k epsilon_ijk epsilon^lmk = delta_i^l delta_j^m - delta_i^m delta_j^l,
sum_ab epsilon_iab epsilon^jab = 2 delta_i^j,
sum_i delta_i^i=3.

The first is an oriented antisymmetrizer resolution, the second a digon factor 2, the last a closed loop factor 3. Altered leg order can change signs; these tensor equations, not an unspecified drawing convention, are the checks used here. Every index instance is verified exactly.

The established composition-closed framework is the oriented A2/SL3 web category at q=1: fundamental/dual boundary words, delta/epsilon generators, linear combinations, contraction, and local relations. Gluing can create loops and cycles, so a tree-only set of diagrams is not closed. Kuperberg proves a generator-and-relation presentation and invariant-space bases; we rely on that theorem for full web completeness, not merely on our small tests. Our index identities should not be mistaken for a separately proved complete list of planar web relations.

By Schur-Weyl duality, the invariants in V^tensor(3m) carry the rectangular Specht module (m,m,m). The hook-length product gives

dim Inv(V^tensor(3m)) = 2(3m)!/[m!(m+1)!(m+2)!]
                         = 1,5,42,462,6006,... .

For nine ports there are 9!/[(3!)^3 3!]=280 unordered triple partitions, but only 42 invariant dimensions. This is a theorem-derived count, NOT a nine-port enumeration performed by this checker. A full rooted ternary tree has 2I+1 leaves; the six-port example is a two-baryon forest, not such a tree. Neither these dimensions nor the occurrence 42 establish equality with ternary-tree enumeration.

### Corollary: unbounded linear gluing memory

Let H_m=Inv(V^tensor(3m)), of dimension d_m=2(3m)!/[m!(m+1)!(m+2)!]. Suppose a linear summary F:H_m -> C^D allows recovery of every future dual gluing <w,v>, for every w in H_m and every input v, from F(v) (with the testing vector w supplied). Then F is injective, and D >= d_m.

Proof. If v lies in ker F, then F(v)=F(0). Every recoverable gluing must consequently satisfy <w,v>=<w,0>=0. Choosing w=v gives ||v||^2=0, hence v=0. Rank-nullity gives D >= d_m. Finally,

d_(m+1)/d_m = (3m+1)(3m+2)(3m+3)/[(m+1)(m+2)(m+3)] > 1

for m>=1, with ratio tending to 27, so d_m is unbounded. Thus no fixed-dimensional linear summary preserves all these unrestricted gluing tests across all boundary sizes.

The quantifiers matter: this is not an obstruction to arbitrary nonlinear encodings, to summaries for a restricted state or test family, or to dimensions increasing with m. It is not a claim about physical QCD entropy. Its content is exactly linear observational sufficiency under all allowed dual gluings.

## 6. Comparison with the audited tournament path-cover carrier

Prior repository work, audited THM3121, supplies the following result: if pc_T(d) counts unordered spanning d-path covers, W_Q=1^T X(I-A_Q X)^-1 1, then

pc_{Q[S_i]}(d)=sum_c ((product c_i!)/d!) [x^c]W_Q^d product pc_{S_i}(c_i).

The division by d! is justified after decoration by distinct actual factor paths, not on repeating bare port words. The C3 expression is (x+y+z+xy+yz+zx+3xyz)/(1-xyz); its three-C3 substitution has Hamiltonian-path count 3159. These are credited prior repository results, not new proofs in this note.

The justified shared principle is observational sufficiency under allowed composition: tournament substitution retains full path-cover data and named ownership; color composition retains typed boundary ports, intertwiner multiplicities, signs, and the exact Gram pairing. A scalar Hamiltonian count or a triality label generally loses information needed by later compositions. There is no claimed algebra isomorphism between these two carriers.

## Exact verification and primary references

check_color_fusion.py uses exact integers/rationals/radicals and explicit exceptions, never assertion gates. Normal and Python -O executions give identical substantive output. It checks all 729 six-color assignments, exact Lie-generator kernels, the pointwise relation, Gram identities, a K5-edge isomorphism, the two mixed channels, octet recoupling, and all local epsilon/delta index cases. Four deliberate false summaries are rejected.

- Greg Kuperberg, Spiders for rank 2 Lie algebras, Communications in Mathematical Physics 180 (1996),109-151: https://arxiv.org/abs/q-alg/9712003 . Primary theorem establishing rank-two spiders and invariant-space bases.
- Shayne Waldron, On the construction of equiangular frames from graphs, Linear Algebra and its Applications 431 (2009),2228-2242: https://www.math.auckland.ac.nz/~waldron/Preprints/Equiangular/equiangular.pdf . Primary author manuscript for the frame/Seidel connection, including Petersen examples.
- David Kabat, particle physics lecture notes, SU(3) invariant-tensor conventions: https://www.lehman.edu/faculty/kabat/particles.pdf . University primary lecture notes; distinguishes delta with one upper/one lower index from a noninvariant same-variance delta.

The invariant-dimension hook-length derivation and all finite identities are mathematical derivations in this report; the checker independently validates dimensions 1 and 5 only. No finite test establishes a new QCD dynamics theorem.
