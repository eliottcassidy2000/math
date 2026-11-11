# Canonical recursion, spectral interfaces, and the prime reservoir

**Status:** PROVED elementary laws below, using the inherited canonical strong-module decomposition theorem. This is a complete normal form for finite tournaments, not a finite catalogue of prime tournaments, dynamical behaviors, or Hamiltonian-path values. Chan–Pak remains OPEN. No external novelty claim.

## Inheritance and scope

The canonical modular decomposition is an inherited input, as used by
[THM-1960](../../../../01-canon/theorems/THM-1960-tournaments-compose-from-regular-seeds-the-spectral-substitution-law.md)
and [Session268 modular §2](../../session268_composable_tournaments_20261111/modular/PROOF.md).
The latter's two historical October26 source paths are absent in this checkout;
we do not claim to have inspected them. The available theorem and Session268
state the inherited mechanism needed here. Regular-seed skew spectral splitting
is already THM-1960; all-word F2 splitting is already
[Session268 dynamics](../../session268_composable_tournaments_20261111/dynamics/PROOF.md).

Closest mechanisms: strong-module trees; full ordered-path-cover substitution;
centered-coronal switching in
[THM-3315](../../../../01-canon/theorems/THM-3315-tournament-cut-switching-centered-coronal-walk-compiler.md).
Hostile: a root change in the same TT2 becomes spectrally visible after a single
seam. Corrected near miss: scalar characteristic polynomial does not classify
F2 dynamics, and a bounded prime alphabet does not approach the global census.
Least-used sidecar: a matrix of resolvent entries at named ports. Anchor:
classify recursion precisely. Niche: irregular-child spectral composition.
Wildcard: finite-rank network scattering on the same tournament carrier.

## 1. A complete normal form, with an infinite alphabet

A module M is nonempty and each vertex outside M has one orientation to all
vertices of M. A strong module overlaps no module: if it meets N, then one
contains the other. A prime tournament here has order at least three and no
proper module of size between two and n−1. We keep TT2 in the linear constructor,
not in the prime alphabet.

**Inherited decomposition theorem.** For n≥2, maximal proper strong modules
partition the vertices. Their quotient is either transitive or prime. In the
transitive case its children are exactly the uniquely ordered SCCs. Recursing
on children yields a unique decorated tree up to isomorphism of prime quotients
with their child decorations. Linear nodes have at least two children and no
linear child; prime nodes carry their actual oriented quotient and at least
three children. Singleton leaves terminate the recursion.

**Normal-form consequence.** Equality of these decorated trees is equivalent to
isomorphism of the underlying finite tournaments. This is a completeness claim
about a representation, not about a finite list of allowed prime decorations.

Proof of consequence: graph isomorphisms preserve the module predicate, hence
strongness, maximality, and quotient arcs. They induce recursive decorated-tree
isomorphisms. Conversely recursively chosen child isomorphisms, permuted by a
quotient isomorphism, preserve every internal and cross-child arc. Induct on n.

Substitution is associative when a nested expression is flattened with its
port addresses retained. It is not commutative. Different expression trees
can describe the same graph, especially parenthesized linear sums; canonical
normalization removes that ambiguity. The tournament substitution operad is
a convenient syntax, whereas the decorated strong-module tree is its canonical
object normal form. A seam edit is not an operation of that operad.

### Automorphisms and converse

Let T=Q[T1,...,Tk] be a canonical prime node, and let
G=Stab_Aut(Q)([T1],...,[Tk]), stabilizing ACTUAL child isomorphism classes.
Then

    1 -> product_i Aut(Ti) -> Aut(T) -> G -> 1,
    |Aut(T)| = |G| product_i |Aut(Ti)|.

This extension splits after choosing consistent identifications of all copies
of each child isomorphism class. The resulting semidirect product is generally
not canonically split. At a linear root G is trivial, even when children are
isomorphic, because the quotient order cannot be permuted. Proof: restrict an
automorphism to the canonical partition; its kernel acts independently inside
children. Every stabilizing quotient permutation lifts using the chosen
identifications. No scalar-response equality may replace child isomorphism in
G. The observable stabilizer may be larger than the graph stabilizer.

Every tournament automorphism group has odd order. Indeed, an involution moves
some pair u,v into each other and would reverse their arc. Cauchy's elementary
group theorem therefore excludes even order.

Converse sends Q[Ti] to Q^op[Ti^op]. At linear roots it reverses the child order;
at prime roots self-converseness asks for a quotient isomorphism to Q^op
matching each child to the converse of the corresponding child. This gives a
recursive, exact chirality test. Neither q nor adjacency χ sees converse.

### Species interface

For any specified prime alphabet P, let A be nonempty tournaments and S strong
ones within its substitution closure. Then the unambiguous grammar is

    A = SEQ_nonempty(S),
    S = X + sum_(Q in P) Q[A,...,A]/Aut(Q).

For all finite tournaments take the infinite alphabet of every prime class.
At every fixed degree only finitely many terms contribute. Ordinary counting
uses cycle-index/Burnside action on actual child types. Reciprocal-automorphism
mass uses the corresponding labelled-species division by |Aut(Q)|. Session268
already gives the response-weighted refinement; this is its common carrier,
not a new scalar enumeration formula.

## 2. Finite alphabets are an exponentially thin part of the census

**Theorem S1.** Let P contain s fixed finite prime tournaments and permit any
number of substitution levels and arbitrary linear roots. Its number c_P(n)
of ordinary classes satisfies, for n≥1,

    c_P(n) <= [16(s+1)]^n.

Proof: a reduced tree with n leaves, no unary vertices, and internal arity at
least two has at most 2n−1 vertices. A plane rooted tree on v vertices has at
most 4^(v−1) encodings. Summing v≤2n−1 gives fewer than 16^n encodings. Each
internal vertex can be marked as linear or one of s prime types, contributing
at most (s+1)^(n−1). Fix a labelled representative of each prime and order its
children by those ports. Every canonical class has such a plane encoding; the
map need not be injective. Overcounting is safe and gives the bound.

Since a(n)≥2^binom(n,2)/n!, without ANY asymptotic theorem about A000568,

    c_P(n)/a(n) <= [16(s+1)]^n n! / 2^binom(n,2)
                = exp(−(log 2)n²/2 + O_P(n log n)).

The same upper bound controls the uniform labelled probability of this family.
A finite alphabet can produce rich symmetries, arbitrarily long computations,
and enormous numbers of objects while missing essentially the entire ambient
ensemble. Bounded prime order implies a finite alphabet, so has the same
obstruction. This does not prohibit a finite algorithm recognizing/generating
an infinite prime alphabet; it prohibits a fixed finite seed catalogue under
uniform substitution from being globally exhaustive.

**Theorem S2 (labelled prime reservoir).** A uniformly random labelled tournament
on n≥5 vertices satisfies

    Pr(nonprime) <= [4 binom(n,2)+4n+64] 2^(−n).

Proof: for a fixed S of size k, each of n−k outside vertices has 2 favorable
orientation patterns out of 2^k. Distinct outside vertices use disjoint edges,
so

    Pr(S is a module) = 2^(−(k−1)(n−k)).

Union-bound k=2 and k=n−1 separately: their contribution is
[binom(n,2)+n]2^(−n+2). For 3≤k≤n−2, (k−1)(n−k)≥2n−6. There are at most 2^n
subsets in this middle range, giving 64·2^(−n). This bound is deliberately
simple rather than sharp at small n. It proves the claim and existence of
prime tournaments of unbounded order.

For comparison the expected number of twin-pair modules is exactly
binom(n,2)2^(−n+2). A twin pair cannot be swapped by an automorphism: its internal
arc fixes the order. Thus the leading decomposability mechanism need not
create ANY symmetry. It is on the n²2^(−n) expectation scale, whereas the
leading automorphism correction in the forecast lane is on n³2^(−2n).
Do not identify these two rare-event populations. The expectation is not
silently declared an exact probability asymptotic here.

The parent session's ensemble-transport law converts this labelled result to
ordinary-class probability with an additional error at most δ/(1+δ), where
δ=a(n)/(2^binom(n,2)/n!)−1. See the session synthesis/forecast bridge for that
separate derivation. This section does not re-prove it.

### Sparse seams cannot repair the finite-alphabet entropy deficit

**Theorem S2b (edit-distance barrier).** Let F_P be the preceding fixed finite-
prime-alphabet family, including all transitive tournaments. Let d(T,F_P) be
the minimum number of arc reversals to an n-vertex family member, minimizing
over relabellings. For uniformly labelled T,

    d(T,F_P)/binom(n,2) -> 1/2 in probability.

Proof: put N=binom(n,2), K=16(s+1). There are at most n!K^n labelled family
members. A Hamming ball of radius floor(αN), 0≤α<1/2, has at most
2^(N h2(α)) elements, where h2 is binary entropy (the usual binomial entropy
bound). Therefore

    Pr(d(T,F_P) <= αN)
      <= n! K^n 2^(−N(1−h2(α))) -> 0.

For the opposite bound choose a uniformly random linear order of the vertices.
Its associated transitive tournament disagrees with T on N/2 arcs on average,
so some order disagrees on at most floor(N/2). That transitive tournament lies
in F_P. Thus d≤N/2 deterministically, proving convergence.

For a fixed budget t, the number of edited ordinary classes is at most
K^n sum_(j≤t) binom(N,j)≤K^n(t+1)n^(2t). The entropy proof additionally rules
out every o(n²) seam budget as a generic completion mechanism. An expression
that permits fresh arbitrary large primes is outside this bound. So is a
construction with Θ(n²) independent seam choices. Infinite depth, spectacular
local dynamics, and a few nonuniform ports do not erase this entropy barrier.

## 3. Spectral substitution without regularity: χ alone closes

Use row-to-column adjacency A (A_uv=1 for u→v), any field F, and an
indeterminate z. Work in F(z), not at a numerical eigenvalue. Put

    f_T(z)=det(zI−A_T),   r_T(z)=1^t(zI−A_T)^−1 1,
    g_T(z)=(-1)^n f_T(−z−1)−f_T(z).

**Theorem S3 (tournament coronal reciprocity).**

    r_T = g_T/f_T.

Proof: A+A^t=J−I, so (z+1)I+A=zI−A^t+J. The matrix determinant lemma gives

    (-1)^n f_T(−z−1)
      =det((z+1)I+A)
      =det(zI−A^t)(1+1^t(zI−A^t)^−1 1)
      =f_T(z)(1+r_T(z)).

The transpose leaves the scalar coronal unchanged. This proof works in
characteristic two as well; minus signs then reduce accordingly. The polynomial
g_T has degree at most n−1. No regularity assumption is used.

**Theorem S4 (polynomial quotient determinant).** For every tournament Q with
adjacency B and arbitrary nonempty children Ti,

    f_(Q[Ti])(z) = det(diag(f_Ti(z)) − B diag(g_Ti(z))).       (1)

Proof: let D=diag(A_Ti) and let U have child-indicator columns. Then
A_T=D+UBU^t and U^t(zI−D)^−1U=diag(r_Ti). The determinant lemma gives

    f_T = product_i f_Ti * det(I−B diag(r_Ti)).

Multiplying column i by f_Ti yields (1), a polynomial identity with no
spurious poles. The coronal of the parent, if wanted, follows again from S3;
alternatively it is 1^t R(I−BR)^−1 1 with R=diag(r_Ti).

Consequences:

- Adjacency-cospectral tournaments are interchangeable in EVERY uniform
  substitution context for the adjacency characteristic polynomial. The
  converse follows by the identity context. Thus χ is the exact contextual
  equivalence class for this observer, just as full q is for H in Session268.
- For a linear quotient, f_T=product f_Ti. For C3,
  f_T=f1 f2 f3−g1 g2 g3. The spectral interface is commutative in the latter
  three responses, despite the oriented quotient's automorphism group C3.
- Regular identical children recover familiar mode splitting, but irregular
  children now have an exact interface rather than a failed block-eigenvalue
  guess. The k-by-k determinant cost is separate from polynomial bit complexity
  and the cost of obtaining the prime quotient's data.
- Reducing (1) modulo p commutes with substitution. Characteristic polynomials
  do NOT supply invariant factors/Jordan blocks, hence do not classify F2
  transient heights or periods. Those require a stronger state, such as the
  direct-sum similarity theorem of Session268 in its regular-odd regime.

This extends the practical interface beyond THM-1960's regular-child splitting;
the determinant lemma itself and coronal methods are classical and inherited.

## 4. A seam changes the required observer: matrix scattering

Let A be any tournament adjacency, E a matrix whose columns select specified
ports, and K a matrix of changes supported on those ports, with A'=A+EKE^t
still a tournament. Set C=E^t(zI−A)^−1E. Then

    f_A' = f_A det(I−KC),
    C' = C + C K(I−C K)^−1 C = (I−C K)^−1 C.               (2)

Proof: matrix determinant lemma and the Woodbury resolvent identity. The two
forms in (2) agree by the geometric-series identity in a formal expansion at
infinity. They hold as rational functions, with denominator cancellation in
the final characteristic polynomial. For reversing an arc a→b, use ports
(a,b) and K=[[0,−1],[1,0]], a rank-at-most-two update.

A useful complete spectral wiring interface for a child is

    (f_T, G_T), G_T=U_T^t(zI−A_T)^−1 U_T,

where U_T contains the all-ones column and every named boundary vertex.
For disjoint children assemble G=diag(G_Ti). Every uniform quotient connection
is a coupling between all-ones columns; every finite seam is a coupling
between selected vertex columns. With their combined coupling K,

    f_total = product_i f_Ti det(I−KG),
    G_wired = G + G K(I−G K)^−1G.

Select or sum columns with a fixed matrix P to keep the next boundary:
G_out=P^t G_wired P. In particular the new all-ones column is the sum of child
all-ones columns. This is an exact compositional scattering system on the
actual adjacency operator. Duplicate/dependent observer columns cause no
problem. It is not an electrical/passive network claim.

**Sharp small hostile.** Start with TT2 on 0→1 followed by sink 2. Choosing
port 0 and reversing 0→2 gives C3 with f=z³−1. Choosing port 1 and reversing
1→2 gives TT3 with f=z³. The unrooted child, its χ, q, and automorphism group
are identical in the two inputs. Only the chosen boundary differs. Thus χ
alone ceases to be a congruence once arbitrary rooted seams are admitted.

No new graph-isomorphism recursion follows from (2): seams may destroy the
old canonical modules and create new automorphisms. Re-normalize the actual
graph before ordinary-class counting. Rooted constructions can have many
presentations of the same unmarked object.

## 5. Hamiltonian boundaries are different from spectral boundaries

The full unrooted q_r, the number of ordered r-path covers, closes for uniform
substitution. For marked seams a sufficient exact state is the following tensor:
for each r, count ordered path covers refined by the ordered list of start/end
labels of all r paths, with labels from {named ports,*}; * means an unmarked
vertex. Every named vertex occurs at most once in the cover, and endpoint
coincidences for singleton paths are retained.

To wire children, enumerate quotient run words and substitute these refined
covers. Test each cross-run adjacency against the uniform quotient orientation,
except where the two endpoint labels specify a seam override. The ordered
endpoint labels tell exactly which transitions were deleted or added. This
is an exact bijective compiler: cut a finished path cover into maximal child
runs to invert it. A port lying internally on a run cannot participate in a
cross-child seam and correctly contributes no endpoint label.

This state is sufficient, not claimed minimal. Its run index r is unbounded
with child size even with a fixed number of marked vertices. Fixed boundary
cardinality alone therefore does not justify a finite-state Hamiltonian
compiler. Session270 obtains a finite 17-dimensional rational representation
for a particular ordered-spine language by stronger structural restrictions;
that is not an all-tournament theorem.

For one reversed cross arc, this tensor contracts to Session268's exact
H(A⋈B)=H(A)H(B)−e_A(a)s_B(b)+L_A(a)R_B(b). The subtraction is an actual
forbidden transition, not an abstract signed numerical gadget.

## 6. Completeness boundary

These results completely specify (a) the finite-tournament normal form given
all prime decorations, (b) its automorphism recursion, (c) the χ observer under
uniform substitution and finite named seams, and (d) a sufficient H boundary
interface. They do not classify every prime, solve isomorphism from χ or q,
prove the Chan–Pak spectrum, give a fixed finite automaton for all contexts,
or classify infinite dynamical systems arising from arbitrary update rules.
That stronger phrase 'all emergent structures' has no finite exhaustive
meaning without fixing the object language, observer, equivalence, and scale.
The atlas fixes those coordinates and identifies which classifications close.
