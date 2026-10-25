# A real geometric reflection that cannot descend through a fixed Fano marking

Status: PROVED elementary connection, with finite controls. Source integration, not a historical novelty claim.

## 1. The two source objects

Incoming commit82691d391b4beef9682ef956526391bef5f91948, `05-knowledge/results/reflection_bridge_20261025/seven.md`, §§2–4, constructs the seven-set of Boolean subsets of{p,q,r} after identifying empty and full. Modified complement fixes the boundary class and interchanges each singleton with its complementary pair. It has cycle type1·2³. It is not an order-reversing map of the nonempty-subset poset and is not a Fano collineation.

The audited Session60 hourglass arrangement is E,E±x,E±z,E±(x+z). Its intrinsic central line E is fixed by the projective base transformation A:(E,x,z)↦(E,−x,−z), which exchanges the three opposite outer pairs. The normalized forms are literally permuted, so a pure permutation of the seven root coordinates lifts A to an involution of the hourglass minimal resolution Y.

Choose the following explicit identification (a marking): boundary class↦E; singletons p,q,r↦E+x,E+z,E+x+z; complementary pairs↦the respective opposite forms. Under this marking, modified Boolean complement is exactly the geometric reflection A on the branch labels. This preserves the fixed central object and the three pairs. It does not transport Boolean inclusion/divisibility or manufacture Fano incidence. The marking of the three pairs is necessary; their coefficients and linear dependencies remain additional geometric data.

## 2. No Fano marking is invariant under this reflection

**Proposition.** An involution of the seven-point Fano plane induced by a nonidentity element of GL(3,F2) fixes exactly three points. Consequently no permutation of cycle type1·2³ can preserve ANY Fano-plane structure on these seven labels.

Proof. In characteristic2 an involution T has T=I+N, with N²=0 and N≠0. Since im(N)⊂ker(N) in a three-dimensional space, rank(N)=1. The fixed space ker(N) has dimension2, hence contains exactly three nonzero vectors. The remaining four vectors form two pairs. Conjugating by a relabeling cannot change cycle type. QED.

In the audited thirty-kernel construction, each kernel defines a Fano incidence on the seven inertia labels. A label permutation stabilizes the kernel precisely when it preserves that incidence (equivalently, induces an invertible map of the rank-three quotient deck group). Therefore A fixes NONE of the thirty kernels: it pairs them into fifteen pairs, within the same geometric quotient type.

This does not obstruct A as an automorphism of the original hourglass. Rather, A gives an isomorphism Y/K→Y/A(K), not a descended automorphism of a fixed marked quotient when A(K)≠K. This is a concrete instance where an actual geometric symmetry exchanges valid quotient structures instead of preserving one. It agrees with the even D12 orbit sizes4,12,6,6,2, but the fixed-point-free argument is independent of that census.

## 3. The canonical representation distinguishes another kind of involution

The audited symmetry theorem gives det(A) times the permutation action on the seven canonical forms. Here det(−I_2)=+1. Thus A has four +1 eigenvectors and three −1 eigenvectors on H^0(Y,K_Y), by its one fixed coordinate and three transpositions. In particular the resolution of Y/⟨A⟩ has geometric genus4: finite-quotient cohomology H²(O) is the invariant part by the Reynolds projector, quotient singularities are rational, and resolution preserves H²(O); Serre duality gives the claimed geometric genus.

This is different from projective DECK involutions, whose invariant geometric genera are1,3,5 according to their even-representative weights6,4,2. A branch-permuting involution and a sign-changing deck involution must not be identified just because both have order two. No assertion about the canonical self-intersection, full automorphism group, or minimal birational type of Y/⟨A⟩ is made here.

## 4. Exact boundaries

This is an equivalence of explicitly marked seven-sets with involution, plus an actual geometric lift and a proved obstruction to Fano descent. It is not an equivalence of Boolean lattices, divisor arithmetic, Fano incidence, or the historical Keller3–1–3 inverse fibers. In particular incoming seven.md's failure of a second-level Keller reflection remains untouched.

No integer-square point or tournament-path assertion follows. The useful retained data are the marking, the linear base action, its determinant, and the chosen quotient kernel.
