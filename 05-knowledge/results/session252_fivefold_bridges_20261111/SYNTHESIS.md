# Five cores, golden geometry, and the algebra behind the shared numbers

**Session252, 11 November 2026. COMPLETE research deliverable.** New
proved-in-note structural identities, recovered older mechanisms, and small
exact construction controls. No independent canon audit or priority claim.
**Chan–Pak and five-core all-height coverage remain OPEN.**

The useful conclusion is richer than five things being in bijection. The core
bank sits simultaneously in a path-interface algebra, a finite/affine triangle
overlap diagram, a golden cycle-packing family, and an explicit binary arithmetic
compiler. Those are actual maps. The historical atlas also recovers many fives
that are valuable comparisons but do not yet have a map to the bank.

## 1. The most useful new theorem: nine coordinates become seven

The inherited core bank is

    K1, C3, S4, S5a, S5b,
    with H = 1,3,5,9,9.

Its last two cores cannot be identified merely because their counts agree.
Session251 supplies the exact compiler with arbitrary binary flank addresses
A,B and a core cut. This session proves the universal identity

    K_C(−1)=−H(C), for every nonempty tournament C.

Here K counts directed paths weighted by H of their unused core. The proof is
an odd-cycle antipode: in squarefree variables, the generating function of all
induced-subtournament H values is an exponential of odd-degree monomials, so
h_C(x)h_C(−x)=1. Extracting the full vertex coefficient proves the identity.

Consequently the apparent interaction of the two address parity bits vanishes
universally. Only seven coordinates are needed:

    (d,e,u,u',v,v',q)
      =(H,H+K(1),O(1),O(−1),I(1),I(−1),J(1)+J(−1)).

The child count is exactly

    d(g(A)+g(B)−1)+eAB+vA−v'(A mod2)+uB−u'(B mod2)+q.

The seven coordinates also have an explicit associative ordinal-gluing law.
They suffice after every future gluing operation in this language, not just
for one final query. Seven actual core/cut responses have determinant−48,
and the seven query functions are linearly independent, so this is a minimal
linear carrier over Q for the displayed response family.

This is a practical reduction toward the open cover: fewer response directions
must survive a target-reducing rewrite. It is not a finite-state or all-height
cover theorem. [Proof, algebra, and exact basis](algebra/PROOF.md).

## 2. A genuine five-dimensional algebra and a genuine pentagon

Before evaluating at ±1, every coloured core has the matrix

    [[H, O(z), J(z)], [0, H+K(z), I(z)], [0, 0, H]].

Ordinal gluing is exactly multiplication. This is a five-dimensional algebra:
the centralizer of the square-zero matrix E13, equivalently the endomorphism
algebra of the dual-number module k[ε]/ε² ⊕ k. Thus there really is a natural
five-coordinate mathematical object behind the compiler. It exists for every
core and is not the assertion that its five coordinates correspond to the
five chosen templates.

Four ordered blocks have five Catalan bracketings. Their lawful reassociations
form the A2 associahedral pentagon, and all five evaluate to the same full
marked interface. The regular pentagon and its golden triangles give a geometric
representation of this actual flip graph. The coefficients and cut marks must
remain attached: a bare pentagon forgets the counts it is meant to compose.

This is the cleanest link to other five-element collections: a finite alphabet
of templates can have a composition algebra whose coherence is itself pentagonal.
[Path classification and pentagon](golden/PROOF.md#1-polynomial-transfer-matrix-and-its-pentagon).

## 3. The two H9 cores lie on a finite/affine boundary

Use an intrinsic graph, rather than inventing an ordering of the five cores.
Its vertices are the core's directed triangles, and adjacency means sharing a
directed arc. For the five literal cores the resulting graphs are

    empty, A1, A2, A3, affine A2.

In S5a the triangles are 012,123,234: a chain. In S5b they are 012,123,124:
all three pass through arc1→2. Their adjacency spectral radii are respectively
sqrt2 and2. Hence the first Cartan matrix 2I−A is positive definite, and the
second has a one-dimensional kernel. This is a genuine finite/affine ADE
observable on the two equal-H cores.

Both full odd-cycle conflict graphs are K4. Their H values and their entire
univariate packing polynomials agree. The difference only becomes visible
when cycle overlap remembers arcs, not merely whether cycles meet. That is
exactly the kind of information loss the compiler's coloured endpoints avoid.
We have not proved that this overlap graph determines the full cut language.

A recovered Platonic map sharpens the picture. Two directed triangles on five
labels meeting at one vertex determine a double transposition in A5, hence an
icosahedral half-turn axis. S5a has exactly one such marker, (01)(34), fixing2;
S5b has none. This is a direct five-label-to-icosahedron incidence map, not
one Platonic solid assigned to each core. [Core/ADE and Platonic details](golden/PROOF.md).

## 4. A single tournament carries both the golden ratio and Fermat prime17

Extend the first four cores along the canonical tournament T_n whose consecutive
edges point forward and all longer edges point backward. Its directed cycles
are exactly consecutive intervals. This recovers an old tower mechanism,
including THM-4952 and THM-4865, rather than rediscovering it under a new name.

The triangle arc-overlap graph is the path A_(n−2). At n=6 it is A4, with
Perron vector (1,φ,φ,1) and eigenvalueφ. Folding by reversal gives
λa=b, λb=a+b, hence λ²=λ+1. The corresponding side lengths (φ,φ,1) form the
golden triangle. The very same six-vertex tournament has H=17.

There is an even more informative common carrier. Its odd-cycle packing
polynomial Z_n(x) has generating function

    sum Z_n(x)z^n=(1−z²)/(1−z−z²−(x−1)z³).

At fugacity1 it gives Fibonacci numbers; at fugacity2 it gives the actual
tribonacci Hamiltonian counts 1,1,3,5,9,17,31,57,105,.... Golden and tournament
growth are two observations of the same interval-packing object. We do not
pretend H itself is Fibonacci, or that every later count is a Fermat prime.

The resulting useful picture is a finite-type stem plus an affine branch.
The all-height question is whether those different incidence patterns supply
all needed arithmetic responses. The inherited holes21901 and86641 show that
neither H9 branch can simply be discarded from the current bank.

## 5. Fermat arithmetic supplies an exact theorem, and exposes the missing lift

A binary word has a four-coordinate linear state (g(n),n,n mod2,1), or five
coordinates when its length is retained. The two-flank tensor then calculates
both H and graph order. The known Fermat factors give clocks2,4,8,16,32 and a
CRT decomposition modulo

    3·5·17·257·65537=2^32−1.

This is a structural arithmetic link, but an explicit new construction reveals
its limit. Let B_k be the binary word (100) repeated 2^s k times and put
M_s=2^(2^s)−1. Then

    B_k≡0 mod M_s,
    g(B_k)≡1+(2^(s+2)/7)k mod M_s.

The slope is a unit, so **K1 alone already covers every residue modulo every
such Fermat-product modulus**. In particular, five-prime CRT fullness cannot
explain why the two H9 cores repair genuine integer holes. It even accepts
residues7 and21.

The positive lesson is a precise next obligation: obtain a small-output lift.
For M_s>N, a legal witness with H≡N mod M_s and 0≤H<M_s must have H=N.
Our modular witnesses are enormous; they do not provide that bound. The five
cores should therefore be studied as ways to lower the integer cost within
residue fibres, rather than as ways to obtain residues already available.

The recovered Lucas/Gilbreath identity also gives actual flank addresses:
W(t)=product of Fermat numbers indexed by the1-bits of t. Its two bit
substitutions feed directly into the graph constructor, independent of
whether later Fermat numbers are prime. [Arithmetic proofs and constructions](arithmetic/PROOF.md).

## 6. What the older archive adds, including the Ogg correction

The recovery scanned329 thematic candidate files and followed35 explicitly
recorded source paths, with additional proof routes in the three notes. It
reaches March–November, including old notebooks rather than only recent
session syntheses. Two sparsely materialized packets were recovered directly
from Git. A prior atlas's much larger search counts are clearly attributed
to that older atlas, not claimed as new reading here.

The most useful recovered collections are:

- Five Platonic solids, with one self-dual solid and two dual pairs, governed
  by curvature and three rotation-group families. Their classification is
  not a demonstrated symmetry of the core bank.
- Five known Fermat primes, with binary constructibility and exact run clocks.
  Their exhaustiveness remains open.
- The genuine modular/Monster pentad {2,3,5,7,13}: ordinary genus-zero X0(p)
  levels, B12 denominator support, and primes whose cubes divide |Monster|.
  The full Ogg set has15 primes and uses the Fricke quotient X0(p)+ instead.
- A5's five conjugacy classes, including two distinct size12 classes; Petersen's
  five maximum independent stars; five exceptional Lie types; five Bravais
  lattice types; five partitions of4; five historical fugacity evaluations;
  and several operation-defined or finite-census fives.
- The golden A4 Fourier embedding, whose visible contraction hides expansion
  in its conjugate coordinate. It supplies exactly the warning needed for a
  proposed decreasing inverse: prove proper height in the full state.

[Full historical atlas and source trails](HISTORICAL_ATLAS.md).

## 7. The next research programme

1. Use the seven-state gluing law to build target-preserving rewrites between
   the two H9 response languages. Keep actual cuts and row masks.
2. Measure a proper height that includes address length, integer output and
   hidden carry. Test whether the affine branch can lower this height where
   the finite-type stem cannot.
3. Use the costed binary tensor and Fermat moduli as a search organizer for
   bounded lifts, not a substitute for them. Lucas addresses give a specific
   structured sublanguage to test.
4. Retain cycle-by-arc incidence and the icosahedral marker when comparing
   cores. A collision of full response languages is the decisive test, not
   a matching H value or a common five-element count.

These are concrete construction and proof-reduction tasks. No universal
rewrite, bounded lift, minimal-five theorem, A5 action on the bank, or completed
Chan–Pak proof is claimed. The inherited five-core finite range remains262143.

The three small reproducible engines exercise new identities only; they do
not rerun inherited large coverage censuses. Both ordinary and optimized
Python runs retain explicit checks and agree. Documentation routing passes.
Methods used: recover the statement before deriving it again; retain local
incidence through counting; separate local support from bounded-height coverage.
