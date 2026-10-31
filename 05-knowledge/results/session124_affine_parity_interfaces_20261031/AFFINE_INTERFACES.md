# Affine parity interfaces form a standard dagger compact closed category

**Purpose and scope.** This gives a precise bridge from the proved parity-interface calculations to closed symmetric monoidal categories. It is a standard categorical model, not a new category-theoretic result, not a realizability topos, and not a claim that Hamiltonian construction or connectedness is affine.

The direct-sum compact structure for finite-dimensional **linear relations** is established in Baez and Erbele, [*Categories in Control*](https://math.ucr.edu/home/baez/control.pdf), especially Section2. Affine extensions of graphical linear algebra are studied by Bonchi, Piedeleu, Sobocinski and Zanasi, [*Graphical Affine Algebra*](https://discovery.ucl.ac.uk/id/eprint/10081075/). The elementary affine-relation construction needed here is proved directly below; no novelty is claimed.

## 1. Definition and closure

Use objects `V=F2^n`, n>=0. A morphism `R:V->W` is either the empty set or an affine subspace of `V direct-sum W`. Equivalently it is the solution set, possibly empty, of a finite system of affine equations in the input/output coordinates.

Composition is ordinary relational composition:

`S after R = {(v,u): there exists w with (v,w) in R and (w,u) in S}`.

The identity is the diagonal `{(v,v):v in V}`. Composition imposes the two systems simultaneously and eliminates the shared coordinates. An intersection of affine constraints is affine or empty, and a linear projection of an affine space is affine. Thus composition stays in the stated class. Associativity follows by reassociating existential quantifiers; diagonal relations give the identity laws.

**Empty relations are essential.** Preparing the state0 in F2 and then testing for state1 yields the empty scalar relation. The class of nonempty affine relations would not be closed under composition.

The tensor on objects is direct sum. On relations it places two independent systems side by side, with the coordinate permutation grouping the inputs and outputs. The unit is F2^0. Coordinate swaps give the symmetry. Relational converse `R dagger={(w,v):(v,w) in R}` is affine, involutive, reverses composition and preserves this tensor.

Here the zero-dimensional vector space is the monoidal unit; it is not a categorical zero object.

## 2. Explicit self-duality and snake equations

Every object is its own dual. Define the cup and cap by

`eta_V:0->V+V`, `eta_V={(0,(v,v)):v in V}`,

`epsilon_V:V+V->0`, `epsilon_V=eta_V dagger`.

They are linear, hence allowed affine relations. For the first snake, an input x passes through cup-tensor-identity to `(a,a,x)`. Applying identity-tensor-cap requires a=x and returns a=x. The composite is exactly the identity on V. The other snake sends x to `(x,a,a)`, requires x=a, and returns x. There are no additional outputs or missing inputs.

The cups and caps respect direct sums after the standard middle-coordinate swap, and cap is cup's dagger. Therefore the category is dagger compact closed, in particular symmetric monoidal closed.

## 3. Internal hom and currying, explicitly

Choose the internal-hom object `[V,W]=V+W`, using V's self-duality. For an affine relation `R:U+V->W`, its curry is

`curry(R):U->V+W`, with `(u,(v,w)) in curry(R)` exactly when `((u,v),w) in R`.

This is simply reassociation of the same coordinate set; other coordinate conventions express it as a permutation. It preserves affine spaces and the empty relation. The inverse is the opposite reassociation.

Evaluation `ev:(V+W)+V->W` consists of `(((v,w),v),w)` for v,w. Composing `curry(R) tensor identity_V` with evaluation reconstructs R, because the existentially chosen first V coordinate must equal the supplied V input. These assignments give the natural bijection

`AffRel(U+V,W) = AffRel(U,[V,W])`.

The internal hom is an object of the relation category, not the vector space of affine functions. Its states are affine relations, not merely point-valued maps.

## 4. Exact relation to the proved parity interfaces

Let A and B be piece indices and S their shared labels. With common reference input parity, compatibility requires

`u(a(s))+v(b(s))=0` over F2 for every s in S.

This defines a **linear relation** between the color vectors in F2^A and F2^B. Its solution space has one free bit per connected component of the piece-overlap graph, as established by the earlier pushout lemma.

Adding fixed colors gives affine equations, such as `u_i=1`. A fixed nonzero color vector is an affine point-state `{c}:0->F2^A`, but it is not a linear subspace when c!=0. Thus affine rather than solely linear relations are the appropriate exact class when these assignments are included.

**Notation clarification:** Baez–Erbele write `FinRel_k` for finite-dimensional linear relations. In that convention AffRel is the needed enlargement. If `FinRel` instead means finite sets and arbitrary relations, that broader category already contains all these affine relations; nonzero points do not require enlarging it.

Disjoint independent blocks correspond to tensor/direct sum. Gluing an internal color interface corresponds exactly to existential elimination, hence relational composition. It may produce an affine relation or detect inconsistency by producing the empty relation. Associativity and currying are therefore valid rewiring laws for this parity-constraint layer.

## 5. Hard boundary: affine compatibility is not Hamiltonian realization

A satisfying color assignment does not produce an ordering of atoms, enforce numerical square connectors, establish exact interval ownership, or prove that a traversal is connected. All those obligations remain separately certified. The affine relation is exact for the stated parity constraints; as a model of possible Hamiltonian words it is generally only a necessary-condition relaxation, not automatically the exact image of all realizations.

For a concrete warning, K4 has three undirected Hamiltonian cycle edge-sets. Their indicator vectors form a set of cardinality3, whereas every nonempty affine F2-subspace has cardinality a power of2. Thus Hamiltonian feasibility is not an affine relation in general. No category of all Hamiltonian constructions, or composition law preserving their existence, has been established here.

This AffRel category is also **not a topos**: it has no terminal object. For any candidate object W, there are at least two distinct morphisms from the monoidal unit to W—the empty relation and the zero point—so W cannot be terminal. Consequently the present model supplies no realizability-topos claim. Nor is its direct-sum tensor a categorical product in this relation category.

## 6. Exact small-dimensional checks

`test_affine_relations.py` uses only the Python standard library, enumerating affine subspaces by generating linear subspaces and then their cosets, with the empty set added. It checks all object dimensions0,1,2:

-143,332 composable relation pairs for affine closure and dagger reversal;
- all six snake equations (both directions for each object);
-36,042 affine relations for currying/evaluation and the inverse bijection; the reassociated ambient dimension can reach6;
-3,912 relation triples for exhaustive associativity with each object dimension at most1;
- the inconsistent point/effect example and the K4 nonaffinity example.

Normal and optimized runs are byte-identical, with **364,040 explicit checks**. Affine-set counts including empty in ambient dimensions0 through6 are2,4,12,52,308,2452,26388. These finite checks supplement the general proofs, rather than establishing category laws by sampling.

Reproduce with `python test_affine_relations.py` and `python -O test_affine_relations.py`. No frozen square-sum or mixed-radix publication artifact is changed by this standalone note.
