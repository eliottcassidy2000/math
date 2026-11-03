# Path-cover interfaces form a concrete symmetric monoidal closed category

**Status.** Exact constructions and proofs. This gives two distinct categorical products behind S164, an explicit positive internal-hom formula, and a genuine slice-category model of moving guards. It does not assert that every formal interface is the path-cover interface of a tournament. The category-theoretic machinery is standard; the formulas below specialize it to the path-cover calculus.

## 1. The interface carried by a finite directed graph

Let FI be a skeleton of the category of finite sets and injections. Write

\[
 h_k(S)=\operatorname{Inj}([k],S),\qquad
 \mathcal E=[\mathrm{FI},\mathbf{Set}].
\]

This covariant functor category is also the presheaf topos on FI-opposite. Let **P** be its full subcategory of finite coproducts of representables:

\[
 A\cong\coprod_{k\ge0}p_kh_k,\qquad p_k\in\mathbb N,
\]

with only finitely many nonzero coefficients. Define its cover-coordinate profile and slot count by

\[
 P_A(z)=\sum_kp_kz^k,\qquad |A([r])|=F_A(r)=\sum_kp_k(r)_k.
\]

The coefficients are uniquely determined: the falling factorials are linearly independent. Conversely the displayed profile specifies the object up to FI-natural isomorphism.

For a finite directed graph X, define C_X(S) to consist of a cover c of its vertices by disjoint nonempty directed paths, together with an injection from the set of paths of c into S. An injection of slot sets relabels occupied slots. Thus C_X is an FI-set, with

\[
 C_X\cong\coprod_{k\ge0}p_k(X)h_k.
\]

A graph isomorphism transports each path and its assigned slot; hence X↦C_X is a genuine functor on the groupoid of finite directed graphs and isomorphisms. Empty X gives h_0. For nonempty tournaments, p_0=0.

## 2. Two products, not one

### 2.1 Shared slots: cartesian product is the matching product

Given an injection [k]→S and an injection [l]→S, their images meet in t slots. The overlap determines a partial bijection between t selected input coordinates. There are

\[
 \binom{k}{t}\binom{l}{t}t!
\]

possible partial bijections, and their union occupies k+l−t slots. This decomposition is natural in S and proves

\[
 h_k\times h_l\cong
 \coprod_{t=0}^{\min(k,l)}
 \binom{k}{t}\binom{l}{t}t!\;h_{k+l-t}.
\]

Consequently

\[
 P_{A\times B}=P_A\star P_B
 =\sum_{t\ge0}\frac{z^t}{t!}P_A^{(t)}P_B^{(t)}.
\]

For ordinal join X⊕Y, all cross arcs point from X to Y. Each slot's X-path and Y-path concatenate uniquely, omitting an empty fragment. This is a natural bijection

\[
 C_{X\oplus Y}\cong C_X\times C_Y.
\]

Associativity and the unit follow from literal concatenation. Thus cover interfaces give a strong monoidal functor from ordinal-join digraphs to (P,×,h_0). The source is not symmetric; the symmetric image does not identify the input graphs.

### 2.2 Disjoint slots: Day tensor is ordinary multiplication

The disjoint-union monoidal structure on FI gives Day convolution ⊗ on E. On representables,

\[
 h_k\otimes h_l\cong h_{k+l}.
\]

Distributing over finite coproducts shows that P is closed under ⊗ and

\[
 P_{A\otimes B}=P_AP_B.
\]

Concretely, for the interfaces in P, A⊗B consists of pairs of interface elements whose occupied slot sets are disjoint. The evident inclusion A⊗B→A×B forgets this disjointness requirement.

For digraph disjoint union, a path cannot cross components. Hence

\[
 C_{X\sqcup Y}\cong C_X\otimes C_Y
\]

naturally. This is a strong symmetric monoidal functor. Disjoint union is normally not a tournament, so this second graph operation belongs to the larger digraph category.

Both tensor units are h_0. The Day structure is affine, since h_0 is terminal. It is not cartesian: there is no map h_1→h_2, because an FI-natural transformation h_1→h_2 would, by Yoneda, be an injection [2]→[1]. In contrast there is a cartesian diagonal h_1→h_1×h_1. Sharing and duplicating a slot are allowed by × but not by ⊗.

## 3. A finite positive internal hom

**Theorem.** The finite-interface category (P,⊗,h_0) is symmetric monoidal closed.

For A=∐p_kh_k and B in P, the internal hom is

\[
 [A,B]_{\otimes}(S)
 =\prod_{k\ge0} B(S\sqcup[k])^{p_k}.                 \tag{1}
\]

**Proof.** For a representable input, Day's adjunction gives

\[
 [h_k,B]_{\otimes}(S)=B(S\sqcup[k]).
\]

This identity can also be proved directly by testing against every h_m and using Yoneda:

\[
 \operatorname{Nat}(h_m,[h_k,B])
 =B([m+k])
 =\operatorname{Nat}(h_m\otimes h_k,B).
\]

An internal hom out of a finite coproduct is a finite cartesian product, which yields (1). To see that it remains in P, note the natural decomposition

\[
 h_j(S\sqcup[k])\cong
 \coprod_{t=0}^j\binom jt(k)_{j-t}h_t(S).             \tag{2}
\]

Choose which t argument coordinates land in S and inject the remaining coordinates into the k reserved slots. Finite cartesian products remain in P by the matching decomposition. This proves closure and the adjunction. ∎

Let D=d/dz. Shift by k reserved slots has profile

\[
 P_{B(-\sqcup[k])}=(1+D)^kP_B.
\]

Thus (1) becomes the explicit nonlinear positive formula

\[
 \boxed{\quad
 P_{[A,B]_{\otimes}}
 =\mathop{\star}_{k\ge0}
   \big((1+D)^kP_B\big)^{\star p_k}.
 \quad}                                                   \tag{3}
\]

An empty matching-product is 1. Formula (3) is a new available interface operation, not a claimed tournament construction.

Examples:

- [h_1,h_1]=h_0⊔h_1, profile 1+z. A unary slot program either returns the argument slot or one distinguished ambient slot.
- [h_2,h_1]=2h_0⊔h_1, profile 2+z.
- [h_1,h_2]=2h_1⊔h_2, profile 2z+z².
- If B is a triangle's cover interface, P_B=3z+3z²+z³, then [h_1,B] has profile 3+9z+6z²+z³.

The complete finite set of natural interface maps also has a closed count:

\[
 |\operatorname{Nat}(A,B)|
 =\prod_{k\ge0}F_B(k)^{p_k}
 =F_{[A,B]}(0).                                           \tag{4}
\]

This follows by applying Yoneda separately to every representable summand. Such maps act on cover interfaces, not necessarily on their underlying graphs.

If B has positive degree d and A is nonzero, then

\[
 \deg P_{[A,B]}=d\sum_kp_k=dP_A(1).
\]

This predicts the exact interface-degree cost of forming a higher-order function space before expanding it.

## 4. A shifted derivation and the multiplicative formal group

Because shift preserves cartesian products,

\[
 (1+D)(P\star Q)=((1+D)P)\star((1+D)Q).
\]

Equivalently,

\[
 D(P\star Q)=DP\star Q+P\star DQ+DP\star DQ.                \tag{5}
\]

The extra overlap term is forced by shared-slot semantics. Ordinary multiplication instead has the usual Leibniz rule. Thus the same cover coordinates carry both the ordinary derivative rule for disjoint resources and a multiplicative-formal-group derivative rule for shared resources. Iterated reserved-slot shifts form a monoid of ⋆-semiring endomorphisms.

The currying identity is equally concrete:

\[
 [A\otimes B,C]\cong[A,[B,C]].                             \tag{6}
\]

Formula (3) therefore has exact nested-composition identities, rather than merely resembling a closed structure.

## 5. Why cartesian closure requires larger objects

The ambient E is a cartesian-closed topos, but P is not cartesian closed. Let 2=h_0⊔h_0. The ambient cartesian exponential has

\[
 (2^{h_1})(S)=\operatorname{Nat}(h_{|S|}\times h_1,2).
\]

Since

\[
 h_n\times h_1\cong nh_n\sqcup h_{n+1},
\]

its cardinality is

\[
 |(2^{h_1})([n])|=2^{n+1}.                                \tag{7}
\]

No finite sum of falling factorials equals 2^{n+1} for every n. Any exponential internal to the full subcategory P would have the same values when tested against h_n, so such an exponential does not exist there.

There is an explicit description: an element at S assigns one of two values to every named slot of S and one further value to an unspecified fresh slot. Under S↪T, new slots receive that fresh/default value. The growth obstruction is therefore an exact interface-size statement, not an abstract nonexistence slogan.

This separates three notions:

1. finite positive profiles, which are closed for the disjoint-slot tensor;
2. arbitrary FI-set interfaces, which also support cartesian function spaces;
3. literal graph or tournament realizability, which is an additional condition in either setting.

## 6. Guarded covers live over the slot object

Let P_* be the full subcategory of E/h_1 whose underlying FI-set lies in P. A map e:A→h_1 marks an occupied slot. Every representable component has degree at least one and marks one of its coordinates; all choices of coordinate give isomorphic pointed representables.

Define

\[
 R(A,e)(S)=\{(a,t):a\in A(S),\ t\in S,\ t\ne e(a)\}.
\]

This is a functor on the slice: a pointer-preserving map keeps the inequality valid. Depending on the intended operation, mark either the old slot e(a) or the new chosen slot t. Its underlying profile is

\[
 P_{R(A,e)}=(N-1)P_A,\qquad N=z(1+D),                      \tag{8}
\]

because on a component h_k, a different slot is either one of k−1 other occupied slots or a fresh slot:

\[
 R(h_k,e)\cong(k-1)h_k\sqcup h_{k+1}.
\]

For a unique-in-neighbor guard q in X, S and U are the cover interfaces according as q starts or does not start its path, each pointed by the slot containing q. The defect interface d for X−q is pointed by the slot containing the distinguished predecessor rho. These are objects of P_*.

Let V=h_1 and 1=h_0. Then the S164 one-vertex guard operations admit the following actual coproduct/product constructions:

\[
 C(S,U,d)=
 (R(S)\sqcup(V\times U),\ S\sqcup U\sqcup d,\ S),           \tag{9}
\]

\[
 D(S,U,d)=
 ((V\times S)\sqcup d,\ (V\sqcup1)\times U,\ R(d)).         \tag{10}
\]

Pointer conventions:

- In the first coordinate of C, mark the chosen insertion slot t.
- In the second and third coordinates of C, retain the old marked slot of the relevant summand.
- In D, retain the old marked slot in every summand, including R(d).

The natural bijections are exactly the insertion cases in S164: insertion at an allowed slot start; insertion after q; and the exceptional segment q,a,rho. For D's V×S, the chosen old guard slot indexes insertion after q and the other slots index insertion at their beginnings. For (V⊔1)×U, V indexes all beginning insertions and the extra point indexes insertion after q.

Decategorifying (9)–(10) gives

\[
 C(S,U,d)=((N-1)S+NU,S+U+d,S),
\]

\[
 D(S,U,d)=(NS+d,(N+1)U,(N-1)d).
\]

The positivity of N−1 is now a genuine complement construction using a retained physical pointer. It would be unsafe to call an arbitrary subtraction of FI-sets a functor. Likewise the formal derivative D is not an endofunctor on all natural interface maps: a map may forget the coordinate whose distinguished-slot use defined differentiation. The slice construction avoids that mistake.

The right ordinal graft has an equally literal module-category realization:

\[
 (A,e)\odot B=(A\times B,e\circ\pi_A).
\]

Under the S164 hypotheses this acts componentwise on (S,U,d), and corresponds to (X,q)↦(X⊕Y,q). The module laws follow from cartesian associativity and the marked slot is retained. A nonempty block placed before X need not preserve the unique-in-neighbor condition, so the physical realization remains a right action.

## 7. An actual selector operad and an equality-branching gap

The endomorphism operad of a single occupied slot in the Day-closed category is completely explicit:

\[
 \operatorname{End}_{\otimes}(h_1)(n)
 =\operatorname{Nat}(h_1^{\otimes n},h_1)
 =\operatorname{Nat}(h_n,h_1)
 \cong[n].                                                 \tag{11}
\]

There are exactly n operations of arity n, each selecting one input slot. Under substitution, first select one input block, then select one input within that block. The symmetric-group action relabels the selected input. This is the set operad of pointed finite sets, obtained here as an actual endomorphism operad rather than an analogy. In particular no slot-natural operation on independent slots can inspect their names or manufacture a new slot.

More generally,

\[
 \operatorname{End}_{\otimes}(h_k)(n)
 \cong\operatorname{Inj}([k],[nk]),
\]

with (nk)_k operations. Substitution composes the corresponding blockwise injections.

Switch to shared-slot cartesian inputs. An n-tuple of slots has an equality partition \(\pi\) of [n]; its distinct entries occupy one slot per block. Consequently

\[
 h_1^{\times n}\cong
 \coprod_{\pi\in\operatorname{Part}([n])}h_{|\pi|}.
\]

A natural output chooses one of the occupied blocks separately for each equality pattern. Thus, for n≥1,

\[
 \left|\operatorname{Nat}(h_1^{\times n},h_1)\right|
 =\prod_{\pi\in\operatorname{Part}([n])}|\pi|
 =\prod_{k=1}^n k^{\left\{{n\atop k}\right\}}.               \tag{12}
\]

Here the braces are Stirling numbers of the second kind. At arity three this gives 24 shared-slot operations, compared with only 3 disjoint-slot operations. At arity four it gives 373248, compared with 4. The increase has a precise source: permitted branching on equality patterns, not hidden extra graph constructions.

## 8. Concrete consequences and boundaries

- **Higher-order compilation.** Formula (3) computes a finite positive profile for a family of slot-natural programs, and (4) counts those programs exactly.
- **Resource distinction.** Replacing Day tensor by cartesian product permits reuse of slots. The extra matching terms and the extra term in (5) are the exact cost of that reuse.
- **A rigorous guard type.** The marked occupied slot makes the forbidden-slot operation R functorial, supplying a categorical explanation for the closed positive three-state calculus.
- **A precise limit of finite profiles.** Day closure is finite; cartesian closure already produces 2^(r+1) on a one-slot input. Finite profile algebras and topos-level function spaces are different expressive levels.
- **No unsupported physical claim.** An interface with a nonzero constant term generally is not the cover interface of a nonempty graph. For example [h_1,h_1] has profile 1+z, which cannot be the cover profile of any single finite graph: the empty graph has profile 1, and a nonempty graph has zero constant term.

## Finite controls

Run `python3 check.py` and `python3 -O check.py`. The control checks 64 pairs of small nonnegative profiles, 512 currying and adjunction-cardinality instances, 96 literal injection-overlap cases, and 24 literal pointed-complement cases. These finite controls corroborate the displayed decompositions; the proofs establish the general statements.

## Sources

The tournament formulas being categorified are S164's [profile algebra](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/PROFILE_ALGEBRA.md) and [pointed profile algebra](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/POINTED_PROFILE_ALGEBRA.md).

The general convolution construction is Brian Day, [On closed categories of functors](https://people.math.rochester.edu/faculty/doug/otherpapers/DayReport.pdf), Reports of the Midwest Category Seminar IV, Lecture Notes in Mathematics 137 (1970), pp. 1–38. FI is the standard category used in Church–Ellenberg–Farb, [FI-modules and stability for representations of symmetric groups](https://arxiv.org/abs/1204.4533). Set-valued FI-functors are studied directly in Ramos–Speyer–White, [FI-sets with relations](https://arxiv.org/abs/1804.04238). No claim of literature priority is made.
