# Counting, explicit bijections, and the observer they must carry

## 1. Recover the actual two categories

Let T_n be tournament isomorphism types on n vertices. For an undirected
simple graph X and g in Aut(X), choose any reference orientation and let
epsilon_X(g) be (-1) to the number of edges whose directions g reverses.
This parity does not depend on the chosen reference orientation. Let E_n
be the graph types for which epsilon_X is identically+1.

**CITED:** |T_n|=|E_n| for every n. This is Royle, Praeger, Glasby, Freedman
and Devillers, [Theorem1.1](https://arxiv.org/html/2204.01947v2), checked
against the primary text on November4. It counts unlabelled types of the
same order. The paper presents a counting proof and asks for a natural
bijection; the constructions below are repository continuations with their
own stated scope.

Three parities must stay separate:

|Predicate|Meaning|Decisive example|
|---|---|---|
|Automorphism-sign-even|Every automorphism reverses an even number of oriented edges|P3 is sign-even although two degrees are odd|
|Degree-even / Eulerian|Every vertex has even degree|C3 is Eulerian but sign-odd|
|Even edge count|The total number of edges is even|K1,3 is sign-even but has three edges|

Automorphism sign here is neither vertex-permutation sign nor the sign of
the induced permutation of the edge set. The tournament/graph identity is
not an identity with all Eulerian graphs. At n=4 the respective unlabelled
counts are4 and3. Labelled cardinalities also differ: at n=5 there are1024
tournaments and396 sign-even graphs, despite12 types on both sides.

If “categories” means their isomorphism groupoids, there is an additional
elementary obstruction for every n>=2. A tournament automorphism group has
odd order: an involution would swap the ends of some edge. The empty graph
is sign-even and has automorphism group S_n, which has even order. An
equivalence of groupoids must preserve automorphism groups up to isomorphism,
so these groupoids are not equivalent. This does not prevent bijections of
their sets of isomorphism classes. It is a different assertion from the
more specific equivariant-map obstruction below.

## 2. What the prior bijections actually accomplish

|Map / source|Construction and inverse|Preserved information|Boundary|
|---|---|---|---|
|SCC/Lyndon type codec, October11|Read the ordered strong-component word, factor it into nonincreasing Lyndon words, match weighted atoms to connected sign-even graph types by a specified rank, then take component union; invert the ranks and concatenate factors|Order and the multiset of Lyndon-factor weights as graph component sizes|Depends on atom order; labelled names and local graph operations are not preserved|
|Primitive Fourier matching, October20|A normalized linear coalgebra map identifies primitive supports; a selected supported matching supplies connected atoms|Component/primitive organization, with a specified matching|The linear map is rank-free; the selected set map still uses global support and is not a local rule|
|Prime-local Walsh matching, October20|Choose a nonzero modular matching in the signed-copy matrix, with its inverse permutation|For ONE chosen odd prime p, the entire embedded Sylow-p action; converse/parity refinement is available|A different p can require a different map; simultaneous preservation of all primes fails|
|Prime-circulant derivative map, October20|Differentiate an antiperiodic selector, route its necklace along the specified cut derivative graph, and integrate along the inverse path|Prime-circulant type and affine covariance of the labelled forward map|Forward map is an equivariant SURJECTION inducing a type bijection, not a labelled bijection|
|Sign-aware rooted augmentation, October20|Track the sign character on the CURRENT root stabilizer and rank/unrank legal branches|Marked extension orbits and transporter signs|Relatively even children can arise from globally odd parents; restricting to globally even parents loses branches|

The exact source paths and audit qualifications are in SOURCE_LEDGER.md.
These constructions are real progress, not merely unsuccessful attempts.
They do not define a unique object called “the bijection.” Already at n=3,
the SCC codec sends TT3 to the empty graph and C3 to P3, whereas the Sylow-3
matching sends TT3 to P3 and C3 to the empty graph. Both are inverse pairs
on the two types; they preserve different predicates.

There are also sharp obstructions:

- **PROVED**, THM-4722 and the independently audited uniform obstruction:
  for every n>=5 no S_n-equivariant map from labelled tournaments to labelled
  graphs induces a bijection between all tournament types and sign-even
  graph types. This does not exclude arbitrary type bijections or every
  possible meaning of a natural unlabelled construction.
- **PROVED:** the empty graph admits both full Sylow3 and Sylow5 symmetries,
  whereas tournament stabilizers cannot jointly pay the required actions.
  Equal primewise marginal counts do not give a joint matching.
- **PROVED:** the mixed action C15 on ten names has32 tournament types
  but31 sign-even graph types. Normalizer-relative evenness loses a hidden
  odd part swap after false-twin fusion.
- The four sign-even graph types at n=4 have edge counts0,2,3,4. They cannot
  be four representatives of one affine edge-toggle square: a linear parity
  on a square has0,2,4 ones, not one. A set bijection does not transfer the
  tournament flip algebra.

The surviving target is a specified type map with a specified preserved
operation or observable. A direct local atom rule beyond the known circulant
sector is still an open improvement, not a consequence of equinumerosity.

## 3. The faithful common language is the signed Fourier basis

Fix old physical names and let V=F2^{E(K_n)}. Write t_e(T)=+1 for an upward
arc and -1 otherwise, and chi_X(T)=product_(e in X)t_e(T).
These functions, indexed by ALL labelled graph supports X, form the full
Walsh basis. For any scalar function f on labelled tournaments,

    f(T)=sum_X fhat(X) chi_X(T),
    fhat(X)=2^(-binom(n,2)) sum_T f(T) chi_X(T).          (1)

**INHERITED PROVED**, THM-4721 Theorem W and the October26 guarded-response
note: if f is invariant under a relabelling group Gamma, its supports combine
into signed Gamma-orbit sums. A support orbit contributes one invariant line
iff epsilon_X is trivial on Gamma intersect Aut(X). Otherwise its coefficient
is forced to equal its negative and vanishes. For Gamma=S_n the surviving
supports are precisely E_n. This explains equinumerosity as equality of the
dimensions of two bases of the same invariant-function space.

It is a linear isomorphism, not an assignment of one graph to each tournament.
The canonical positive coupling obtained by squaring normalized Walsh entries
is another valid bridge, but squaring discards signs and gives a probability
distribution over pairings, not a deterministic adjacency map.

**INHERITED trace identity.** Whole converse sends each chi_X to
(-1)^|X| chi_X. In the basis of tournament-type point indicators its trace
is the number s_n of self-converse types. In the graph basis the trace is
the signed edge-parity count. Thus

    s_n = #{X in E_n: |X| even} - #{X in E_n: |X| odd}.  (2)

At n=5:12=10+2 and8=10-2. The converse-even and converse-odd function spaces
have dimensions10 and2. Eight self-converse types are not eight specially
designated graph objects: selecting fixed graph images requires a matching.
The prior prime-local parity construction provides one such choice.

## 4. Specialize this bridge to the full named two-port response

Let F_T(a,c) count actual Hamilton paths after adjoining x->y, with
a_u=[u->x], c_u=[u->y]. Write its UNIQUE Boolean multilinear polynomial as

    F_T(a,c)=sum_(|A|<=2, |C|<=2) f_(A,C)(T) a_A c_C.   (3)

This retains every cubic and quartic and permits shared physical names
between A and C. Each port has at most two incidences in a Hamilton word,
which proves the bidegree bound. It is not an assumption that every partial
reversal contrast is cubic.

**PROVED specialization of the inherited basis.** The coefficient
f_(A,C) is invariant under the group fixing A and C separately as sets.
Its graph basis is therefore relative-even for that marked stabilizer.
Likewise a fixed actual query q=(a,c) uses the subgroup preserving its two
attachment masks. A fixed physical query need not be an S_n-invariant
function, so the unmarked E_n basis alone cannot represent it.

The small exact hostile is particularly direct: the named TT3 with rows
(6,4,0) and its converse (0,1,3) are isomorphic, but at (a,c)=(1,0) their
actual responses are3 and9. Every unmarked invariant agrees on them.
The marked Fourier expansion has a nonzero one-edge support coefficient
-12/8=-3/2. That support is globally sign-odd; the permitted stabilizer
fixes its marked endpoint and no longer contains the offending swap.

**Exact operation law, all orders.** If D is a fixed set of old pairs and
T^D reverses those pairs, then

    chi_X(T^D)=(-1)^|X intersect D| chi_X(T),
    f(T)-f(T^D)=2 sum_(|X intersect D| odd) fhat(X)chi_X(T). (4)

This is just multiplication of the reversed edge signs; it needs no module
hypothesis. Module legality and contextual invariance are separate claims.
For several D_i, mixed differences multiply by the product of the factors
1-(-1)^|X intersect D_i|. Thus the named reversal masks are exactly the
additional information the Fourier operation needs.

Whole converse has D=E(K_n) and selects odd-edge supports. The inherited
quartic cover formula shows every f_(A,C) of attachment degree4 is
converse-even, so ALL its odd-old-edge Fourier coefficients vanish.
Hence the whole-converse contrast is cubic. For a partial reversal an
even-edge support may meet D oddly, so genuine quartics survive. This
explains the earlier mixed quartic -2 without discarding that coordinate.
Attachment degree and old-edge support parity are different gradings.

For the disjoint-module decoder already proved in S1-S3,

    F_e(q_i)=A_i+e_i*d_i
            =(A_i+d_i/2)-(d_i/2)(-1)^e_i,
    d_i/2 is an odd integer.                            (5)

On that reversal cube the decoder has exactly one nonconstant bit-character.
The module theorem supplies this sparsity by exact exterior profile
invariance; generic Fourier expansion does not prove it. At unrestricted
queries, mixed characters and cubic/quartic interactions remain present.

This is the useful merger: use marked signed graph coordinates to transport
the complete observer and its operations, while keeping type codecs as
separate selected set maps. A graph-side response defined by decoding a
type is legitimate but inherits that codec's choices and lacks physical
port labels unless the marking/transporter is also carried.

## 5. What it does and does not say about the spectrum

H is a genuine unmarked invariant, so it has an exact sign-even graph-basis
expansion. Its full two-port refinement is a marked covariant object.
Equal dimensions, an invertible basis matrix, or a type bijection provides
neither positivity of graph-side coefficients nor a graph construction for
each prescribed H. Therefore none proves completeness outside7 and21.

The three strong seed states still supply the actual scalar interval
4*3^r+9*2^r-3+2j, j=0,1,2, r>=1. The ordered-cover skeleton explains its
two exponential channels. The present linear correspondence can encode the
same invariant but does not turn that thin family into all odd integers.
The exact next test for any proposed constructive graph rule is whether
it intertwines the required marked attachment operation and has the correct
integer realization domain, not whether its object count matches A000568.
