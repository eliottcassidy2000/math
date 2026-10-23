# Split real equality: the minimal ambient algebra is E8(8)

Status: PROVED corollary of `UNGRADED_METRIC_COMPLETION.md`, using standard real-form and split-Cartan facts. This separately strengthens the earlier notes' explicit real-form nonclaim. It does not change their complex statements or assert uniqueness of marked embeddings.

Let

    g_R=e7(7)⋉H57(R)

use the real56-dimensional symplectic fundamental module of split e7. Its complexification is the irreducible complex56-dimensional E7 module; in particular it is absolutely irreducible. Suppose a finite-dimensional real Lie algebra L_R contains g_R and has a real nondegenerate invariant symmetric bilinear form.

**Corollary.** dim_R L_R≥248. Equality is attained by split e8(8), and every equality ambient L_R is isomorphic as a real Lie algebra to split e8(8). No grading is assumed.

**Proof.** Complexifying the embedding and metric gives a quadratic complex Lie algebra L_C containing e7⋉H57(C), with dim_C L_C=dim_R L_R. The ungraded theorem gives the lower bound. At equality L_C≅e8.

The real centralizer T_R=C_(L_R)(e7(7)) complexifies to the3-dimensional centralizer T_C≅sl2(C), because centralizers are kernels of linear commutation equations and commute with scalar extension. The real isotypic component is V_R⊗M_R with dim_R M_R=2. This follows either from complete reducibility over R or by complexifying the space of s-equivariant linear maps; absolute irreducibility gives End_s(V_R)=R. The nondegenerate invariant symmetric form restricts as the product of the alternating forms on V_R and M_R.

Consequently T_R acts on M_R through sp(M_R)≅sl2(R). This action is faithful: its complexification is the faithful action T_C≅sl2(C) on the2-dimensional multiplicity module proved in the ungraded theorem. Since both real algebras have dimension3, T_R≅sl2(R). Equivalently, the original nonzero z gives a nonzero nilpotent matrix annihilating the original real multiplicity line; a compact su2 centralizer could not do this.

Choose a7-dimensional split Cartan a in the embedded e7(7), and a nonzero split Cartan generator h in the commuting T_R≅sl2(R). Their sum a⊕Rh has dimension8 because e7 and its centralizer intersect trivially. On every finite-dimensional representation of a split real semisimple algebra, its split Cartan acts by commuting real-diagonalizable operators. Thus ad(a) and ad(h), which commute, are simultaneously diagonalizable over R on L_R. They give a split toral subalgebra of dimension8.

The complex rank of e8 is8, so L_R has full split rank8 and is the split real form e8(8). Conversely the split Chevalley realization of the E8 contact grading contains the specified split e7(7)⋉H57(R) and has nondegenerate Killing form. It attains248. QED.

Scope: the conclusion depends on the split E7 input and its specified real56-dimensional module. It is not a classification of embeddings with arbitrary real forms of E7, and it does not choose a unique embedding or metric scale.
