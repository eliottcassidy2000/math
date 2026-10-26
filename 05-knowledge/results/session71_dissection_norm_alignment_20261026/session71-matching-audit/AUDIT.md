# Session71 independent matching / geometric-carrier audit

## Result and scope

**PASS, with the distinctions below retained.** The literal whole-side edge
matching diagram has an exact two-involution presentation. The two-involution
product alone is not faithful. The marked Penrose construction genuinely
realizes every Boolean conjunction by the same named translations in every
bit state, with exact nonoverlap. Euclidean seam holonomy determines compatible
frames but not a dissection. None of these statements independently proves
the complete Dudeney lower-bound case analysis or an infinite Penrose extension
theorem.

The independent checkers use active exceptions under both normal Python and
`-O`. `check.py` checks every ordered pair of involutions on at most seven
named sides, explicit endpoint-sign models, flat-side and parallel-edge
controls, and the five-square overlap witness. `check_penrose.py` uses a
separate arithmetic implementation over `Q(sqrt(5))`, not the producer's
`Q(phi)` code. It verifies shape angles, simplicity, metric, oriented area,
direct isometries, markings and 1,364 actual state/subset placements through
five gadgets. The all-size assertions also have proofs below; they are not
extrapolations from these finite universes.

Base repository: `c2039c32f9726b09ed0e09d5f5c126fc72975897`.
Only files in `session71-matching-audit` are owned by this audit. No shared
repository edit or commit was made. The read-only processor check was blocked
by the missing `.machine-id`; no identity was fabricated.

## 1. Exact statement: whole sides, not seam atoms

In the paper's Definition 3.1 (HTML Definition 2), let `S` be the disjoint
union of **intrinsic whole side occurrences** of the pieces. For each
assembly `X`, every side determines one node `q_X(s)` of ED, and every ED node
has one or two incident side occurrences. Thus `q_X` is a surjection with
fibers of size one or two. Define a partial fixed-point-free involution
`tau_X` by swapping the elements of each size-two fiber and leaving it
undefined on each singleton. Equivalently totalize it to `alpha_X` by fixing
the singletons.

The map

    ED node v on shore X -> (X, q_X^{-1}(v)),
    ED link s           -> s

is an isomorphism to the bipartite incidence multigraph whose shore-X nodes
are the `alpha_X` orbits and whose named edge `s` joins its two orbits. Its
inverse sends a singleton to a fixed point and a two-element node incidence
to a transposition. The proof is a direct bijection on nodes and links;
parallel edges are not collapsed. The abstract converse holds for every
bipartite multigraph without isolated nodes and with all node degrees one or
two. It does **not** assert that every such graph is geometrically realizable.

Boundary/flat/internal node types, the owner of each side, whole-side lengths,
piece cyclic order, corner incidence, and assembly vertex identifications are
additional data. Fixed points are not automatically outer boundary nodes.

### Flat hostile and correction lineage

My initial startup concern suggested common atomization. Reading the actual
definition corrects that suggestion: **common atomization is not the literal
paper ED.** At a T-junction a long side opposite two shorter sides is represented
by three whole-side links whose nodes at this assembly are singletons; the
long side gets a virtual node. An atomic seam refinement instead splits the
long side and pairs each segment with a short-side atom. These objects have
different link sets and degrees. The paper also removes certain residual
cut edges between oppositely flat endpoints because they are no whole piece
side. A refinement can be useful, but needs a map back and extra grouping
information; it cannot silently replace the object under an isomorphism claim.

### Product-only hostile

On two named sides, compare

    alpha_A = alpha_B = identity,
    alpha_A = alpha_B = (01).

Both have product `alpha_B alpha_A = identity`. The first ED consists of two
one-link paths; the second has two parallel links making one two-link cycle.
So product-only state loses even connected components. A two-step product
can be used as an operation, but the faithful state retains both involutions
(or one plus their product), fixed-point stopping semantics and node labels.
In partial-involution language undefined moves must remain failure, not be
mistaken for actual seam contacts.

There is a partial-product hostile as well: take `tau_A=(01)` and `tau_B`
nowhere defined, versus both nowhere defined. Both `tau_B composed tau_A`
are nowhere defined, but the first ED is one two-link path and the second
is two one-link paths. Product loss is therefore not merely an artifact of
the fixed-point totalization convention.

The producer's double-coset statement is correct for the stated conjugation
action: independently fixed conjugacy classes have relative alignment
`C_G(a) \\ G / C_G(b)`. If assembly decorations vary as well, stabilizers of
the complete decorated states replace plain centralizers. This is a
combinatorial classification, not a geometric sufficiency theorem.

## 2. Endpoint parity, corner incidence and flips

Order a side's two endpoints by the orientation of its piece boundary. Two
pieces meeting on a full seam occupy opposite local halfplanes. If both
placements preserve their chosen boundary orientation, seam crossing swaps
the endpoint slots. Therefore `r` seam-contact steps give slot permutation
`t -> t xor (r mod 2)`. This is the elementary mechanism used by the
paper's no-flip endpoint lemma, with the step count stated explicitly instead
of relying on an ambiguous edge-versus-node path convention.

For general orientation signs `f_X(P)` in assembly X, one contact from side
`s` to side `s'` changes endpoint bit by

    1 xor f_X(owner(s)) xor f_X(owner(s')).

Compose these signed permutations along the path. Under no-flip common
orientation conventions this reduces to the preceding parity rule. If a
piece changes orientation between assemblies, its two occurrence signs are
different and must **not** be canceled as one fixed sign. Already two
contacts, with sign pairs `(0,0)` and `(1,0)`, have net reversal despite an
even contact count. A one-contact geometric control is the unit square at
the origin and a second reference square placed by `(x,y)->(2-x,y)`: their
shared east-side endpoint indices agree, rather than swap.

The split-chain reflection `x -> L-x` has a different parity: it counts
divide-in-two events. No universal equality with seam-contact parity was
proved or used.

Neither the involutions alone nor their product reconstructs VD. The source
must retain the corner set, the two endpoint-to-corner maps for each side,
the piece cyclic boundary order, and the assembly maps from corners to cut
vertices. Full-seam contacts generate some vertex identifications, not all
flat/higher-valence incidence data. This prevents the loss of precisely the
vertex components and angle sums used in the Dudeney argument.

## 3. Holonomy gives frames, not nonoverlap

For a connected declared seam graph, let `h_ij` be a specified isometry from
piece j's frame to piece i's, with `h_ji=h_ij^{-1}`. Placements satisfying
`g_j=g_i h_ij` exist iff every closed-walk product is identity. Necessity is
telescoping; sufficiency follows by multiplying along root paths, with cycle
consistency proving path independence. The placements are unique up to
common left composition. A spanning tree plus its fundamental cycles is a
finite equivalent test. With no flips, the group is `SE(2)`, not all `E(2)`.

This theorem fixes frames only. It does not prove injective development,
interior disjointness, prescribed target coverage, target boundary, topology,
or absence of unlisted contacts.

The exact hostile is five unit squares at origins

    (0,0), (1,0), (1,1), (0,1), (0,0),

with consecutive seams E,N,W,S. At each seam the interiors lie on opposite
sides. Every interior piece uses two distinct local ports. The declared
seam graph is a tree, hence all closed-walk checks hold. Yet the fifth square
coincides with the first, with positive overlap area one. Total piece area is
five and union area four. The graph is the **declared** seam graph, not the
full incidental-contact graph. The E,N,W,W control has five distinct unit
grid cells and no positive-area overlap. This cleanly isolates the missing
global geometric acceptance condition.

## 4. Independent Penrose geometric audit

Primary source: Penrose, [*Pentaplexity*, printed pp.33–35, Fig.9](https://web.ma.utexas.edu/users/radin/Pentaplexity.pdf).
The actual Figure 9 pixels were inspected. The adjacent text explicitly
requires matching vertex colors. The displayed dart has black 72-degree
and reflex 216-degree vertices, with white 36-degree tips; the kite has
white opposite 72- and 144-degree vertices, with black 72-degree flanks.
Only the corresponding finite same-color matching rule is imported here.

Let `phi=(1+sqrt(5))/2`, `q=2-phi`, and let U,V be the unit directions at
plus/minus36 degrees. In oblique coordinates their Gram form is

    N(x,y)=x^2+y^2+(phi-1)xy.

The independent checker reconstructs

    D = [0, U, (q,q), V],          colors [black,white,black,white],
    K = [(1,1), U, (q,q), V],      colors [white,black,white,black].

It verifies exact side lengths `1,phi^{-1},phi^{-1},1`, simple boundaries,
dart angles `(72,36,216,36)`, kite angles `(72,72,144,72)`, and area ratio
phi. Angles are checked by exact Gram cosines together with convex/reflex
turn signs, so a 144-degree and a 216-degree angle cannot be confused by
their equal cosine. Areas in oblique coordinates are multiplied by the same
positive absolute basis determinant, which cancels in the ratio.

For the two candidate poses the checker independently reduces the vertices to

    K_good = [(1,0),(1,-1),(q,1-phi),(0,0)],
    K_bad  = [(0,0),(1,0),(1,1-phi),(phi-1,-1)].

The bad pose arises from matrix `[[1-phi,-1],[1,0]]` and translation
`(phi,-1)`. Its full Gram matrix is preserved and determinant is +1;
this is a rotation, not a reflection. Both kites lie in `-1<=y<=0`, the dart
in `0<=y<=1`, and their common boundary is exactly the full side `[0,U]`.
Both are interior-disjoint from the dart. The good pose matches both seam
endpoint colors; the bad pose mismatches both.

Thus counts, bare adjacency `K2`, shared length and nonoverlap do not decide
this marked local guard. The bad final configuration is not a valid marked
patch; that is the observable failure, not a counterexample to a tiling rule.

### All-size separated conjunction construction

Put gadget j's dart at `3jU`, and stage the bit-selected candidate kite at
`3jU-10V+K_{b_j}`. Initially all 2k pieces are disconnected and locally legal.
Every gadget lies in the affine x-strip `[3j,3j+1]`; these strips are pairwise
disjoint. Darts lie in y-band `[0,1]` and staged candidates in `[-11,-10]`.
Selecting a named subset J means moving precisely its kites by the same
translation `+10V`, regardless of their bits. Their final band is `[-1,0]`.
No inter-gadget collision or accidental contact is possible. Within a
gadget the only final contact is the displayed seam. Even the continuous
translation is collision-free until that seam is reached.

Consequently its validity guard is exactly `g_J(b)=product_{j in J} b_j`.
This is a genuinely geometric realization with real marked kites/darts,
not an analogy inferred from a shared numeral or a bare Boolean model.

The full subset-by-support matrix is the Boolean zeta matrix, triangular
with diagonal one in bitmask order. It has rank `2^k` over every field.
Selected rows of order at most t have rank `sum_{r<=t} binomial(k,r)` for
`0<=t<=k`. Singleton guards separate the `2^k` individual states, giving k
bits in this fixed gauge. The linear law dimension is `2^k` for arbitrary
signed laws; probability laws have a fixed total-mass coordinate and affine
dimension `2^k-1`.

Uniform even/odd parity laws agree on all proper-subset guards, while the
all-coordinate guard differs by `2^{1-k}`. The transfer to Session70 is thus
exact at the level of the stated native guard family. No conclusion is made
about extension of these disconnected placements to infinite Penrose tilings,
minimal carriers for arbitrary tilings, or Dudeney's four-piece lower bound.

## 5. Cross-review of the Dudeney producer's deductions

The component angle-charge deduction is correct: non-target-corner weights
are integral multiples of pi, so regular n-/m-gon corner counts a,b in one
VD component obey `2a/n-2b/m in Z`. If n,m are coprime this is equivalent to
`n|2a` and `m|2b`. The odd-n versus square specialization is therefore valid.
It is only a necessary component constraint, not a piece-count theorem.

The abstract positive split chain has exact state either empty or
`F(x)=epsilon*x+c` on a nonempty open interval I. Composing with
`G(y)=delta*y+d` on J gives

    (delta*epsilon, delta*c+d, I intersect F^{-1}(J)).

The producer's domain inequalities and composition code were reviewed and
agree with this formula, including infinite upper endpoints and empty
intersections. Equality transport is identity. The final positivity of every
intermediate value is an essential invariant; sign and offset alone cannot
replace the domain.

The minimality proof can be strengthened beyond numerical outputs: if named
positive inputs and **all** future positive split totals L are allowed,
Boolean validity alone separates the triples. Unequal domains are separated
immediately. On a common nonempty open domain distinct affine functions
differ at some x; their positive outputs u<v are separated by appending a
split with u<L<v. One succeeds and the other fails. This is not a minimum
theorem for a fixed finite alphabet such as `{1,2,3,4}` without a separate
context-separation proof. The producer was advised of and accepted this
scope condition.

No claim was made that every abstract split word is realized by a polygonal
dissection. The paper's full optimality result remains CITED. Its requested
v3 artifact/version uncertainty and v4 correction notice remain in the
producer report; this audit does not resolve them.

The later closed-chain strengthening also passes: Q-independent allowed totals
force every type to occur equally often in odd/even positions for an even
formal closed chain. This is not positivity. In the paper's normalized lengths
`sigma=3^(1/4), tau=2`, the word `(sigma,tau)^3(tau,sigma)^3` is formally
identity and perfectly coefficient-balanced. Its fifth output would be
`3sigma-4-x<0` for all x>0, since `3sigma<4` follows from `243<256`.
Its legal domain is empty. This is an abstract-chain hostile, not a purported
geometric dissection or a dependency attributed to the paper.

### Equal-area golden norm extension

The Penrose producer's later P5 also passes. Suppose equal-area similar copies
of the native kite and dart both have vertices in a finite abelian field F.
In the abelian compositum `L=F Q(zeta_5)`, their similarity scalars a,b belong
to L by ratios of corresponding edge differences. Reflected templates still
have coordinates in the same cyclotomic field. Equal area forces
`(b/a) conjugate(b/a)=phi`. Complex conjugation commutes with every automorphism
of L, so every Galois conjugate of that norm is nonnegative. The automorphism
of `Q(sqrt(5))` taking phi to `1-phi<0` extends to L, a contradiction.

The proposed escape `K` and `sqrt(phi) D` has degree eight over Q:
`[Q(zeta_5):Q]=4`, and phi is not a square in that field. Any complex square
root of positive phi would be real, hence would itself have the forbidden
complex norm phi. Adjoining it is therefore quadratic. This theorem concerns
**equal-area rescalings** of unlike shapes, not the ordinary unequal-area
native kite/dart pair. It proves no minimum-piece theorem.

## 6. Reproduction and evidence

From this directory:

```sh
python3 -B check.py > normal.json
python3 -B -O check.py > optimized.json
cmp normal.json optimized.json
python3 -B check_penrose.py > penrose.normal.json
python3 -B -O check_penrose.py > penrose.optimized.json
cmp penrose.normal.json penrose.optimized.json
```

Both checkers accept `--mutate` controls that must terminate nonzero in both
modes. These are explicit semantic hostile assertions (and a malformed
involution input), not a claim of general file-corruption detection. The
corresponding `*.corruption.json` files record the exact rejection reasons.
`verify_evidence.py` separately tests actual one-byte corruptions of copies
of archived source/output files and rejects them by frozen SHA256. This is
an integrity check, not a mathematical proof checker. Source/output hashes
and cross-reviewed producer file identities are recorded in `EVIDENCE.json`;
checks are not assertions that disappear under `-O`.

Inheritance used: the current AGENTS/research protocol and active guardrails;
the Session70 joint named-guard synthesis; the fivefold-boundary packet's
separation of local counts, conjugate/window data and genuine geometry.
The board was whole-side incidence, signed endpoint transport, split-domain
memory, relative frame, marked tile guards and global nonoverlap. The decisive
method is to identify the exact observed operation before claiming a quotient
or connection. No new universal meta-pattern or priority claim is asserted.
