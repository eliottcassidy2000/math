# Relative evenness: a Fourier carrier for paired tournaments

**Status:** PROVED here from the stated finite-group/Fourier argument, with independent-audit status supplied separately by the session. FINITE-EXACT tests are scoped below. No external-priority claim is made. This is a decorated extension and application of the existing repository Walsh principle, not a new proof of the original equinumerosity theorem.

Source snapshot: `3eb49fa63d7ff4e58742b1ae84944f832e87a118`. Source paths below identify the dependencies.

## 1. Inheritance and the definition that matters

The published source is Royle, Praeger, Glasby, Freedman and Devillers, *Tournaments and Even Graphs are Equinumerous*, J. Algebraic Combin. 57 (2023), 515–524, [arXiv:2204.01947v2](https://arxiv.org/abs/2204.01947v2). Its abstract and publication record were checked on 2026-10-29.

For a simple graph X on a linearly labelled vertex set, define

    epsilon_X(g) = (-1)^(number of edges {i<j} in X with g(i)>g(j)).

When g preserves X, this is a character of Aut(X), independent of the reference orientation. X is Royle-even when this character is identically +1. This is **not** degree-evenness, and it is **not** the sign of the permutation induced on the edge set. The graph K2 has just one edge, hence trivial edge permutation, but its endpoint swap has Royle sign -1.

Closest inherited mechanisms:

- `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`: twisted fixed-set counting, Sylow refinement and Walsh basis. Its current correction explicitly refutes unrestricted odd-subgroup type refinement at C15 on ten vertices.
- `01-canon/theorems/THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md`: Sylow 3/5 Hall obstruction.
- `05-knowledge/results/session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md`: the constructive class bijection between aligned rooted consecutive pairs and odd-order self-converse tournaments; the retained input automorphism group is an anti-involution centralizer.

Canonical hostile: K2 with and without distinct vertex colors. Least-used sidecar: the stabilizer of the **joint** graph-and-decoration, rather than the graph alone.

## 2. The decorated Walsh theorem

Let G be any subgroup of S_n and D any finite G-set. Let T_n be the labelled tournaments. Consider the diagonal action on T_n × D. For a simple graph X and d in D, call (X,d) **relatively even** if

    epsilon_X(g)=+1 for every g in G fixing both X and d.

**Theorem A.** The number of G-orbits on T_n × D equals the number of G-orbits of relatively even pairs (X,d). More strongly, over Z[1/2], the invariant functions on T_n × D have a basis indexed by those relatively even orbits.

**Proof.** Put m=binom(n,2), and write each tournament as m orientation bits t. For every graph X let chi_X(t)=(-1)^(sum_{e in X} t_e). The functions

    chi_X(t) * 1_{d=d0}

form a basis over Z[1/2]: invert the ordinary Walsh matrix on the tournament coordinate, leaving D untouched. Relabelling permutes these basis vectors up to the cocycle sign epsilon_X(g). On one orbit of these vectors, invariant coefficients are determined by one coefficient. The coefficient is arbitrary precisely when every stabilizer sign is +1; if some sign is -1, invariance forces a=-a and hence a=0. This proves the basis statement. The orbit indicators on T_n × D give the other basis, so their cardinalities agree. QED.

Choose a base point (X,d) for a relatively even orbit O and transport its signs. The invariant signed orbit sum is

    Phi_O(t,d') = sum_{(Y,d') in O} epsilon(g,X) chi_Y(t),

where g sends (X,d) to (Y,d'). This is well defined: changing g by a stabilizer does not change its sign. A different base point can change the entire column's sign; the invariant line is canonical.

For source orbits A and target orbits O, exact orthogonality is

    sum_A |A| Phi_O(A) Phi_P(A) = 2^m |O| [O=P].

Only the tournament bits are Fourier transformed. The scaling is 2^m, **not** 2^(m+log2|D|) or the total labelled input cardinality.

This is an exact change of basis, not a positive set bijection. Negative entries are essential, and inverse coefficients may contain powers of 2 in their denominators. Over characteristic 2 the Walsh transform collapses and this basis argument cannot be imported.

## 3. Paired tournaments and odd self-converse tournaments

An aligned rooted consecutive pair is represented by (T,S,d) where T and S are tournaments on the same n vertices and d is the extra root's orientation mask. Set R=T XOR S. Since the two affine orientation cocycles cancel, R is an ordinary simple graph. Therefore

    (T,S,d) <-> (T,R,d)

is an S_n-equivariant bijection. Apply Theorem A with

    D = {simple graphs R on [n]} × {binary vertex colorings d}.

The result is the following exact counting and invariant-carrier statement:

    aligned rooted pair classes
        <-> basis lines indexed by relatively even triples (X,R,d)
        <-> invariant functions on odd self-converse classes of order 2n+1.

The first and last class sets are genuinely bijective by the inherited Session89 lift. The middle relation is a **linear-basis correspondence**; this theorem does not choose a canonical matching of its individual graph triples to self-converse tournaments.

For a permutation with all cycle lengths odd, let e be its number of unordered-pair orbits and c its number of vertex cycles. The common count is

    F_n(2,1) = sum_{lambda all odd} 2^(2e(lambda)+c(lambda))/z_lambda.

Through n=0,1,2,3,4 the count is 1,2,8,88,2752. All graph triples, without the relative-even filter, instead have 1,2,12,120,3400 orbits.

**A necessary correction to a tempting simplification.** One cannot demand that X alone be an even graph. For n=2, X=K2 is odd, but d coloring its endpoints differently destroys the swapping automorphism, so (X,R,d) is relatively even for either choice of R. Joint ports can expose Fourier coordinates invisible in the unmarked graph. Conversely, forgetting the colors can make those coordinates disappear again. An independent full-universe audit found that imposing the incorrect stronger condition “X itself is even” gives only 1,2,6,60,1092 classes through n=4, rather than the required 1,2,8,88,2752.

## 4. Decoration does not repair natural set bijectivity

**Theorem B.** Let n>=5 and D be a finite S_n-set containing at least one globally fixed decoration. There is no S_n-equivariant map

    T_n × D -> {relatively even (X,d) in Graph_n × D}

which is bijective on isomorphism classes.

**Proof.** First extend the inherited Sylow refinement to decorations. For an odd prime p and a Sylow p-subgroup P of S_n, the N(P)-orbits on P-fixed objects are exactly the full isomorphism classes admitting a conjugate of P. This is Sylow conjugacy inside the stabilizer of the whole decorated object. On graph objects, Frattini applied to the kernel of the Royle character inside that joint stabilizer says N(P)-relative evenness equals full relative evenness. The restricted twisted Burnside calculation of THM-4721 is unchanged: for every h in N(P), both the ordinary tournament fixed count and the signed graph fixed count are multiplied by the same integer |Fix_D(<P,h>)|. Thus the numbers of source and target classes admitting P agree.

For Sylow 3- and 5-subgroups P3,P5 of S_n, no odd-order subgroup can contain conjugates of both. The inherited THM-4722 proof uses overlapping 3- and 5-cycles; their group contains an alternating group on 5,6 or 7 vertices. Every decorated tournament stabilizer is odd. Hence the source families admitting P3 and P5 are disjoint. Their target families overlap: the empty graph with a globally fixed decoration is relatively even and has stabilizer S_n.

Any equivariant map sends an object's stabilizer into its image's stabilizer. The union of the two source families therefore has at least one fewer possible target class than source classes, by the two Sylow equalities. More generally the Hall deficit is at least |D^(S_n)|: each globally fixed decoration gives a different empty-graph target class in the intersection. Hall's necessary condition fails. QED.

This rules out a symmetry-respecting class bijection in the indicated direction. It does not rule out arbitrary bijections between unlabelled classes, signed linear equivalences, or maps retaining extra data. No major open graph-isomorphism problem is solved by this obstruction.

The globally fixed decoration is a sufficient hypothesis for the obstruction, not a necessary characterization of every obstructed D. Theorem A also permits D empty, when both spaces have dimension zero.

For D=Graph_n×binary colors at n=5, the independently enumerated Sylow-fixed universes give

- Sylow 3: 512 source and 512 relatively even target classes;
- Sylow 5: 8 source and 8 relatively even target classes.

The four globally fixed decorations are R empty or complete and d identically 0 or 1. Thus 520 source classes have at most 516 possible target classes. The independent session audit additionally computes the two target families and finds their intersection has exactly four classes at n=5, so this particular union bound is exact. This is not a claim that the full matching deficit is four. These are exact normalizer-orbit computations and independent full-stabilizer checks, not a full n=5 pair census.

## 5. General form: affine binary state spaces

The mechanism is not confined to graph syntax. Let a finite group G act affinely on a finite binary vector space V:

    g.v = L_g v + c_g,   c_(gh)=c_g+L_g c_h.

Partial Fourier transform over V sends its permutation representation into a signed permutation representation on the dual V*. The sign attached to a fixed dual character alpha under g is (-1)^(alpha(c_g)). With any finite decoration G-set D, invariant basis lines are precisely dual-decoration orbits whose joint stabilizers kill this sign character. The proof is the same basis argument as Theorem A. The tournament case is V=edge-bit space and c_g=edge-order inversion cocycle.

A globally fixed affine point a trivializes the cocycle, c_g=a+L_g a. Translation by a then identifies affine and linear G-sets. Conversely an equivariant bijection from the linear G-set V to the affine one must send 0 to a globally fixed affine point. Consequently:

    affine V is G-isomorphic to linear V
      iff the affine action has a fixed point
      iff c is a coboundary.

For odd |G|, a fixed point exists by summing the affine orbit of any point in characteristic two. For even |G| it may or may not exist. The special tournament action fails for every even-order subgroup because an involution would reverse its own transposed pair's arc.

This distinction is the useful unifying theory: a Fourier change of representation always exists after inverting 2, while a symmetry-respecting change of **objects** requires a vanishing cocycle. Decorations alter the stabilizers on which that obstruction is tested. This explains why new joint port information can make an old invisible sign observable without providing a canonical object-level matching.

## 6. Verification and scope

Run:

    python3 check_decorated_walsh.py
    python3 -O check_decorated_walsh.py

The standalone standard-library checker enumerates the full 2^(n²) labelled pair/triple universes through n=4; counts source orbits and relatively even target orbits; checks every source stabilizer is odd; verifies all pairings of the exact signed Fourier columns through n=3; independently constructs the n=5 Sylow normalizer fixed sets; and retains the K2 parity/color hostile examples. It uses explicit runtime checks, not Python assert statements.

Expected output: 11,222 checks; counts above; PASS. Exhaustive finite verification supports but does not replace the all-order proofs. No full n=5 universe, no Monster group computation, no external novelty, and no major open-problem closure is claimed.
