# Critical injectivity, incidence trees, and classification-free monodromy

2026-10-19. **PROVED; frozen after independent review.** A polynomial whose distinct
critical points have distinct images has alternating or symmetric
geometric monodromy, provided it has at least two critical points.
We give a proof using a tree of simultaneously chosen branch-cycle
supports and an elementary 3-cycle argument. This supplies an alternative
proof of the group-identification steps in packets 135 and 137 without
using the classification of finite simple groups. Those frozen proofs
are preserved unchanged. No historical novelty is claimed.

## 1. General theorem

Let K be any characteristic-zero field and let `f∈K[v]` have degree n.
Suppose f has at least two distinct critical points over `K̄`, and that
no two distinct critical points have the same image under f. Write the
multiplicities of the roots of `f'` as `d_1,…,d_q`, where `q≥2`.

**Theorem.** The geometric splitting-field group of `f(V)−Y` over
`K̄(Y)` is

    A_n  if every d_i is even,
    S_n  if at least one d_i is odd.                (1)

In the geometric S_n case the arithmetic group over `K(Y)` is S_n.
In the geometric A_n case n is odd; put

    δ=(−1)^((n−1)/2)n.

The arithmetic group is A_n exactly when δ is a square in K, and S_n
otherwise.

The assertion concerns the generic target Y. It does not by itself
claim preservation under arbitrary substitutions `Y=R(s)`. The actual
fixed-time transfer needed for the source flows is supplied separately
by packets 135 and 137.

## 2. Primitivity from critical injectivity

The required indecomposability argument is elementary. If `f=A∘B`,
with both degrees greater than one, let α be a critical point of A.
Every point of the fiber `B(v)=α` is critical for f and has the same
image `A(α)`. Hence critical injectivity forces

    B(v)−α=c(v−b)^(deg B).

There cannot be two distinct critical points of A: differentiating
these expressions would force their unique centers b to agree, and
then evaluating B there forces the critical points themselves to agree.
With only one critical point, A is a translated pure power. The displayed
formula for B then makes f a translated pure power, which has only one
critical point. This contradicts the hypothesis. Thus f is indecomposable.

The monodromy group G is transitive because `f(V)−Y` is irreducible.
A nontrivial block system would give a proper intermediate field of
`K̄(v)/K̄(f(v))`. Lüroth's theorem and the unique place above infinity
turn such an intermediate field into a nontrivial polynomial
composition. Therefore G is primitive. This standard bridge and the
critical-injectivity condition are discussed in [Kreso–Tichy,
*Diophantine equations and the monodromy groups*, Lemma 3.4 and
Proposition 3.9](https://arxiv.org/pdf/1601.07316); the needed decomposition
obstruction has been proved above rather than imported through a group
classification.

## 3. Simultaneous branch cycles over the complex numbers

First take the constant field to be the complex numbers. Let
`b_1,…,b_q` be the distinct finite critical values. Choose one regular
basepoint Y_0 and one labeling of its n preimages. Choose standard
based puncture loops around the b_i, with compatible noncrossing access
arcs. These loops freely generate the fundamental group of the complex
plane with the b_i removed. Consequently their monodromy permutations

    σ_1,…,σ_q

are a simultaneous generating set for G in that one labeling.
The cover is connected because its domain is the complex plane minus
finitely many points, so their action is transitive. A loop at infinity
is already a product of the finite loops; it is not an additional
independent generator needed for transitivity.

At b_i there is exactly one critical point, of local degree
`e_i=d_i+1`. Thus σ_i is one e_i-cycle, with all other roots fixed.
This is the standard local covering description, or follows by the
formal substitution `Y−b_i=η^(e_i)` and simple-root lifting. Importantly,
we have not chosen unrelated conjugates of local inertia elements:
the compatible based loops already give both the single-cycle property
and simultaneous generation.

For a primary description of this branch-cycle setup, including the
generating tuple, its product relation, and total ramification at
infinity for a polynomial, see [Peter Müller, *Primitive Monodromy
Groups of Polynomials*, section 2.1](https://ypfmde.github.io/Papers/mon.pdf).
Only this standard covering setup is used here, not that paper's
classification of primitive polynomial groups.

## 4. The support-incidence graph is a tree

Form a bipartite graph with:

- n vertices for the labeled sheets;
- q vertices for the finite branch cycles;
- an edge from sheet j to cycle i exactly when `j∈supp(σ_i)`.

This graph is connected: a union of its connected components on the
sheet side is invariant under every σ_i, and their generated group is
transitive. Its number of edges is

    Σ_i e_i = Σ_i(d_i+1) = (n−1)+q,

because the derivative has degree n−1. Its number of vertices is n+q.
Thus it is a connected graph with one fewer edge than vertices, hence
a tree.

Because q≥2, some two branch-cycle vertices are joined through a
common sheet vertex on a path in this tree. Their supports meet at
that sheet and nowhere else: two common sheets would create a
four-cycle in the incidence graph. Choose such cycles σ,τ and write

    supp(σ) ∩ supp(τ) = {c}.

Under the convention that the rightmost permutation acts first,

    στσ⁻¹τ⁻¹ = (c, σ(c), τ(c)).                     (2)

The three displayed points are distinct. Directly, the commutator
sends `c→σ(c)→τ(c)→c` and fixes every other point. Thus G contains an
actual single 3-cycle, regardless of the original critical
multiplicities.

The incidence tree depends on the chosen compatible paths and labeling.
Only its tree property and the resulting group element are needed;
we do not assert a canonical embedded tree or a spatial symmetry of
any source flow.

## 5. An elementary primitive-3-cycle lemma

We include the group argument to avoid any classification theorem or
small-degree exception.

**Lemma.** A finite primitive permutation group containing a 3-cycle
contains the full alternating group on its underlying set.

First, a subgroup generated by 3-cycles is the product of the full
alternating groups on the connected components of the hypergraph of
their three-element supports. To verify this, start with one support U,
on which its cycle generates A_U, and attach intersecting supports.
If a new 3-cycle meets U in two points and adds one point x, conjugating
by A_U and taking inverses gives all 3-cycles supported on x and a
pair of U. Together with A_U these generate `A_(U∪{x})`; A_U is
transitive on unordered pairs even when `|U|=3`.

If the new cycle is `(a x y)` and meets U only in a, choose
`b∈U\{a}` and conjugate by an element of A_U sending a to b. This
gives `(b x y)`, and

    (a x y)(b x y)⁻¹ = (a x b).

The preceding two-point-overlap argument first adjoins x and then y.
This proves the claim by induction along any connected support
hypergraph. It uses the standard elementary fact that all 3-cycles
generate the alternating group.

Now let N be the normal closure in G of the given 3-cycle. It is
generated by its conjugate 3-cycles. Since N is normal, its orbits
form a G-invariant block system. Primitivity and the nontrivial cycle
force N to be transitive. Its support hypergraph is connected, so the
preceding claim gives `N=A_n`. Therefore `A_n≤G`, as required.

Apply the lemma to (2). If every d_i is even, each σ_i has odd length
and is even. Since those cycles generate G, `G≤A_n`, and equality
follows. If some d_i is odd, the corresponding σ_i is odd, so G is
S_n. This proves (1) over the complex numbers in every possible degree.

## 6. Arbitrary algebraically closed characteristic-zero constants

The complex proof suffices for any algebraically closed field Ω of
characteristic zero; the reduction does not assume Ω itself embeds
in the complex numbers.

Choose a finitely generated subfield k_0 of Ω containing the finitely
many coefficients of f. Such a characteristic-zero field embeds in
the complex numbers. Let k be an algebraic closure of k_0 in Ω; it
is algebraically closed and admits a compatible embedding into the
complex numbers. Critical multiplicities and equality or inequality
of their values are unchanged by either embedding.

Let L be the splitting field of `f(V)−Y` over `k(Y)`. It is a field
and, as an algebra over the algebraically closed field k, is
geometrically integral. One standard criterion is [Stacks Project,
Lemma 10.49.3](https://stacks.math.columbia.edu/tag/0FWF): geometric
integrality is equivalent to remaining a domain after algebraic closure
of the ground field, which here is already k.

For any field extension Ω/k, the base change

    L ⊗_(k(Y)) Ω(Y)

is therefore a domain: it is a localization of `L⊗_k Ω`. It is also
finite dimensional over Ω(Y), so it is a field, of the same degree
as `L/k(Y)`. All automorphisms of L extend, and the roots still generate
the base-changed splitting field. Hence the geometric splitting-field
group, with its action on those roots, is unchanged by extension of
algebraically closed constants. Compare the extension to the complex
numbers and the extension to Ω to obtain (1) in general.

Taking Ω=K̄ proves the stated geometric theorem for arbitrary K.

## 7. The arithmetic parity character

If the geometric group is S_n, arithmetic equality is immediate.
Suppose every critical multiplicity is even. The derivative has even
degree n−1, so n is odd. There is a monic polynomial `h∈K[V]` with

    f'(V)=nℓ h(V)²,

where ℓ is the leading coefficient of f. Indeed the unique monic square
root of `f'/(nℓ)` is fixed by every algebraic-closure automorphism over K.
For the monic polynomial `F_Y=(f−Y)/ℓ`,

    Disc(F_Y)=(−1)^(n(n−1)/2)n^n Res(F_Y,h)².

Its square class in `K(Y)` is therefore
`δ=(−1)^((n−1)/2)n`. Since the arithmetic group lies between A_n and
S_n, the discriminant sign character proves the arithmetic statement.
A constant is square in K(Y) exactly when it is square in K.

This arithmetic calculation is generic-target only. It is compatible
with, but does not replace, the separately proved source-field
specialization and nonzero-square-factor checks in packets 135 and 137.

## 8. Consequences for the source-flow classification

The time primitives in packet 135 have derivative critical
multiplicities a and r, with pairwise distinct critical values. The
present theorem yields A_n exactly when a and r are both even, otherwise
S_n, without invoking Jones's theorem.

For packet 137 with `r≥3`, the nonzero critical multiplicities are r−2,
and the multiplicity at zero is b−2 when b≥3. All are even exactly
when b and r are both even. The same argument gives its full A_n/S_n
classification, including the cubic and quintic cases with no special
group-theoretic treatment. Its r=2 pure-power Kummer family has at most
one critical point, so it is correctly outside the present theorem.

The linear-disjointness transfers to fixed time, individual-coordinate
generation, exact arithmetic square classes, and radical conclusions
remain those already proved in packets 135 and 137. Together with this
alternative group argument, those conclusions have proofs independent
of the classification of finite simple groups.

The hypothesis of at least two critical points is essential: pure
powers have cyclic monodromy. Critical injectivity is also essential:
colliding critical values can produce several nontrivial cycles in one
inertia permutation, invalidate the incidence count used here, and
permit decomposable or other monodromy groups. No statement about all
polynomials, arbitrary rational covers, numerical specializations, or
planar JC follows.

## Verification scope

The [small exact checker](../../04-computation/experiments/polynomial_critical_injectivity_incidence_tree_20261019.py)
passed 13 specified permutation and graph witnesses, including the
commutator orientation and both 3-cycle attachment cases. Its
[normal output](polynomial_critical_injectivity_incidence_tree_20261019.json)
and [optimized output](polynomial_critical_injectivity_incidence_tree_20261019_optimized.json)
are byte-identical. These checks do not replace the universal branch-cycle,
field-extension, or group arguments.

The [independent review](polynomial_critical_injectivity_incidence_tree_20261019_independent_review.md)
accepted the complete proof without corrections and separately checked
the cited branch-cycle and geometric-integrality sources. Final changes
after its pinned draft concern only this review metadata and status.
