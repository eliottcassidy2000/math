# Unit-distance reduction keeps edges but loses the first arithmetic carry

Date: 2026-10-19. Status: **PROVED** for the elementary results below;
**FINITE-EXACT** for the accompanying checks; **CITED / source-checked** for
external results. No new value of the chromatic number of the plane is claimed.

## The useful outcome

Two plausible shortcuts in the current field programme can be separated cleanly.

1. Combining the eight conjugation-fixed places above 29 does not improve on
   one place, even when different finite precision depths are used. The exact
   joint quotient is a categorical graph product, and its chromatic number is
   exactly that of the deepest single-place quotient. This is an elementary
   consequence of simultaneous norm-one lifting and compatible diagonals.
2. Passing from residue precision one to precision two is genuinely different.
   At every odd residue characteristic, the residue unit-distance graph has no
   graph-homomorphic section into its unramified precision-two graph. A single
   straight p-cycle is already an obstruction. Its first carry has trace 2,
   while all allowed first-order corrections are tangent and have trace zero.

The second point does **not** say that chromatic number decreases, that the
residue graph has no unrelated embedding upstairs, or that a prescribed finite
Euclidean graph cannot lift. It says that retaining every residue vertex and
every residue edge simultaneously is impossible.

A third result identifies exactly what all fixed-place finite observations
miss: at any prescribed split prime pair, norm-one directions in the kernel of
those observations have arbitrarily large pole order. Denominator complexity
cannot be inferred from finitely many such residue digits.

A final barrier rules out a tempting repair: stress-free Hensel lifting can
only produce graphs with at most four colors. Five- and six-chromatic critical
frameworks necessarily have nontrivial, quantitatively large stress spaces.
The four-color ceiling is sharp, as an exactly checked Moser spindle shows.

## 1. Inheritance and boundary

Closest proved mechanisms:

- THM-4558, `residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md`:
  norm-one directions are integral at conjugation-fixed places; residue
  colorings extend by components.
- THM-4628, `residue-colourings-of-the-unit-distance-graph-of-q-zeta55-are-kloosterman-forced-at-every-conjugation-fixed-prime.md`,
  and `tree_resonance_composite_pell_20261013b.md`, sections 2–3:
  prime-power stationary phase, simultaneous norm-one lifting, and the CRT
  categorical-product description are already present. We do not claim these
  mechanisms as new.
- THM-4700, `the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md`:
  the zero component is the stated localization; the first excluded rational
  denominator prime is 29.
- THM-4701, `kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md`:
  explicit finite-field configurations and the golden Petersen reduction.
  Nothing here repackages its Petersen construction.

The current `merge_18_19_11_22_collatz_hn_20261018.md`, section 2.2, explicitly
leaves products of the eight 29-places unanalyzed. Section 2 below answers that
question, including unequal finite depths. The older composite-ideal note
establishes the tensor graph and its Hoffman bound; it does not state this
exact homogeneous-local-family chromatic equality.

Canonical hostile: a homomorphism onto a high-chromatic quotient does not give
a lower bound for the original graph. A second hostile is confusing
simultaneous coordinate surjectivity with simultaneous norm-one surjectivity.
Both issues are addressed explicitly below.

Least-used sidecar: the first-order carry/stress compatibility equation. The
external note `decalion89/chromatic-number-of-the-plane/notes/local_colourings.md`,
section 7, reports inconsistent section systems at q=5,11,17. The universal
odd-characteristic p-cycle proof below supplies their common mechanism.

Concept board: norm-one torus; component ring; tensor quotient; first carry and
stress; split-prime pole divisor; cyclotomic direction classes. Anchor: exact
residue-coloring limits. Niche: mixed-characteristic lifting. Wildcard:
CM-integrality and rational-angle coloring.

## 2. Equal local fields: many places are no stronger than the deepest one

Let L be a number field in C, stable under complex conjugation c, with
[L:L+] = 2. Let P_1,...,P_s be distinct c-fixed primes. Assume every completion
L_(P_i) is, as a field with involution, the unramified quadratic extension of
the same unramified extension of Q_p with residue field F_q. This includes
all eight fixed places above 29 in the current degree-48 field, by its inherited
local decomposition. Write S for the common quadratic local integer ring.

For n>=1, let

    A_n = S/p^n S,
    U_n = {u in A_n : u c(u)=1},
    G_(q,n) = Cay(A_n,U_n).

Let U(L)={u in L:u c(u)=1}, and let Gamma be its additive span. For positive
depths n_i set N=max n_i. Every unit step is integral at these places:
since c(P_i)=P_i, one has 2 v_(P_i)(u)=v_(P_i)(u c(u))=0. Hence Gamma is
also integral there. Reduction maps Gamma to the product of the A_(n_i).

**Theorem 1.** This reduction has full vertex image and full unit-step image.
Consequently the graph controlling all colorings of Gamma that factor through
these residues is exactly

    G_(q,n_1) x ... x G_(q,n_s),

where x is the categorical product: an edge changes every coordinate by a
unit step. Its chromatic number is exactly chi(G_(q,N)).

### Proof: full unit-step image

Fix arbitrary norm-one residues u_i. In the local ring A_(n_i), choose t_i so
that w_i=t_i+u_i c(t_i) is a unit. Such a choice exists already modulo p:
the F_q-linear map t -> t+u_i c(t) cannot vanish identically. If it did,
evaluating at 1 and then at arbitrary t would make the nontrivial quadratic
residue involution equal to the identity. Moreover

    u_i c(w_i) = u_i c(t_i)+u_i c(u_i)t_i = w_i.

Thus u_i=w_i/c(w_i). CRT gives one global integral element w with all the
specified w_i residues modulo P_i^(n_i). Then w/c(w) is an **exact** global
norm-one element with the desired residues. It is a unit at every P_i.
This is the norm-one lifting step, not merely CRT on vertices.

### Proof: full vertex image even on Gamma

Modulo p, the norm-one group U_1 in F_(q^2) additively spans the whole field.
Indeed, its F_p-span is a subring (products of generators stay in U_1), hence a
subfield; a proper subfield of F_(p^(2f)) has at most p^f=q elements and cannot
contain the q+1 distinct norm-one elements.

Its differences span the same field: because sum_(u in U_1) u=0 and
q+1=1 in characteristic p,

    sum_(u in U_1) (u-u_0) = -u_0.

The full unit-step image permits two tuples differing at only one coordinate.
Their difference isolates any chosen element of U_1-U_1 in that coordinate.
Therefore the image H of Gamma surjects onto the whole product modulo p.
The target A=product A_(n_i) is a finite abelian p-group, and H+pA=A implies
H=A: repeatedly substitute to obtain A=H+p^N A=H. This proves full vertices.
Every target edge now lifts to an actual edge x -> x+u in Gamma.

### Proof: exact chromatic number

Choose an index j with n_j=N. Projection to factor j is a graph homomorphism.
Conversely the synchronized map

    a in A_N  ->  (a mod p^(n_1), ..., a mod p^(n_s))

is a graph homomorphism into the product. Use the fixed isomorphisms of the
local rings with involution. The two maps give both chromatic inequalities.
There is no appeal to a general product chromatic conjecture.

**Application at 29.** At depth one, all eight places together cost exactly
kappa(29), the same as one place. At arbitrary finite depths they cost exactly
the deepest one-place value. Thus the only possible benefit within this
homogeneous 29-adic family comes from greater precision, not from combining
more of the eight places. The value at greater precision is not determined here.

The local arithmetic in this application can also be checked directly.
For K=Q(zeta15,sqrt(-11)), Frobenius at 29 is complex conjugation, so K
has eight unramified degree-two primes. The real cubic factor is

    X^3+2X^2-4X-4 = (X-3)(X^2+5X+11) modulo 29.

The quadratic discriminant is 10, a nonsquare in F_29, while the cubic
discriminant is 592, prime to 29. Over F_(29^2) all three roots are distinct.
Thus F=K(w) has 24 unramified degree-two places above 29: eight associated
with w=3 are c-fixed; the other sixteen are exchanged in eight conjugate
pairs. The eight fixed completions have the involution, local norm, and
precision quotients used here. This also supplies same-rational-prime examples
for Theorem 3.

**Limits.** This theorem does not compare different residue field sizes or
nonisomorphic local fields. It does not establish chi(Gamma)=chi(G_(q,N)).
Nor does it concern arbitrary colorings of a small selected finite configuration.
The phrase “factoring through the residues” must be retained.

## 3. A straight p-cycle cannot lift through the first carry

Let R be the integers of a finite unramified extension of Q_p, with p odd,
and let S/R be its unramified quadratic extension. Put k=R/pR and
ell=S/pS. The following proof works for every q=|k|, not only prime q.

**Theorem 2.** The reduction map G_(q,2) -> G_(q,1) has no graph-homomorphic
section. More specifically, no closed walk of p unit steps modulo p^2 can
have every directed step congruent to the same u in ell with u c(u)=1.

### Proof

Choose one norm-one lift U of u modulo p^2. Such lifts exist, for example by
the local norm-one construction in section 2. Write the hypothetical steps

    d_i = U + p h_i,  0<=i<p.

Their norm conditions give

    Tr(c(u) h_i) = 0  in k.

Here and below h_i is reduced modulo p. Since the steps close,

    0 = sum_i d_i = p U + p sum_i h_i  modulo p^2,
    u + sum_i h_i = 0  in ell.

Apply Tr(c(u) times -). The left side has trace 2+0=2, a contradiction.
The residue vertices 0,u,2u,...,(p-1)u form precisely such a directed
p-cycle. Any section preserving every edge would lift this cycle. QED.

### Every depth, and a quantitative loss

The same obstruction holds for a section from depth r to depth r+1, for
every r>=1. Use the straight cycle of length p^r and write its hypothetical
steps U+p^r h_i. The norm expansion remains linear modulo p^(r+1), and
closure again gives u+sum h_i=0 modulo p with trace 2. Thus no successive
reduction in this unramified tower has a graph-homomorphic section.

At depth one, the straight p-cycles partition the undirected edges: for each
pair of unit directions {u,-u}, take all affine F_p-lines parallel to u.
There are q^2/p such cycles in each direction, and (q+1)/2 directions.
Every section must fail at least one edge on each cycle. Therefore

    failed unit edges >= q^2(q+1)/(2p) = |E(G_(q,1))|/p.

For q=p=29, at least **435 of the 12,615 residue edges** must fail. This
counts failed adjacency equations for a chosen section, not missing edges
of an independently chosen Euclidean configuration. The bound is attained
at q=3 (6 failures among 18 edges) and q=5 (15 among 75), by the exact
section witnesses checked by `verify_section_witnesses.py`.

At general depth r the analogous cycle partition gives at least
|E(G_(q,r))|/p^r failures for every section to the next depth.

A decisive hostile to an overstrong interpretation is q=3:
**chi(G_(3,r))=3 for every r**. Reduction gives the three-color upper bound,
and the straight odd cycle of length 3^r gives the lower bound. Hence
chromatic number can remain unchanged at every depth while no successive
section exists.

### The literal carry certificate

Take u=1 and the usual representatives 0,1,...,p-1. The first p-1 edge
increments are 1, while the closing increment is 1-p. Their squared-norm
errors divided by p are respectively 0,...,0,-2 modulo p. The first-order
edge rows telescope to zero, but their right-hand sides sum to 2. Thus the
inconsistent system has an explicit one-row dependency certificate, for
all odd p. At p=29 it has rank 28 and augmented rank 29.

### General homogeneous version

The same proof is not specific to the norm form. Let F be a homogeneous
polynomial of degree d over R, with p not dividing d, and let F(u)=1 in the
residue field. Consider the directed Cayley rule F(y-x)=1. There is no
closed p-step walk upstairs all of whose steps reduce to u.

Choose any lift U, write F(U)=1+p a, and write the steps U+p h_i. The norm
linearization is replaced by

    dF_u(h_i) = -a.

Closure gives sum h_i=-u. Summing and using Euler's identity yields

    0 = -p a = dF_u(sum h_i) = -dF_u(u) = -d,

a contradiction. No smoothness or chosen exact lift is needed for this version.

### Sharp boundaries

- If p divides d, the argument can fail completely. For F(x)=x^p and p odd,
  the p vertices 0,...,p-1 lift the directed residue cycle: both 1 and 1-p
  have p-th power 1 modulo p^2. This is checked exactly.
- In equal characteristic, F_q[epsilon]/epsilon^2 admits the constant section.
  The obstruction is the mixed-characteristic carry pU, not just “one more
  infinitesimal digit.”
- Ramified mixed-characteristic rings are outside the stated theorem.
- A graph-homomorphic section is a prescribed lift of every residue vertex.
  Its nonexistence does not rule out other homomorphisms, other realizations,
  or preservation of chromatic number.

## 4. The exact stress test for an individual finite residue framework

This is the standard first-order rigidity equation with its arithmetic carry
retained. Let x_v in ell be a finite residue framework with unit edges, choose
arbitrary lifts X_v in S/p^2, and orient its edges. For e=(v,w), put

    a_e = x_w-x_v,
    c_e = (N(X_w-X_v)-1)/p modulo p.

A lift X_v+p h_v preserves every unit edge **if and only if**

    Tr(c(a_e)(h_w-h_v)) = -c_e  for every edge e.       (L)

All terms omitted from this expansion are divisible by p^2, so this is an
exact criterion at the stated precision, not an asymptotic approximation.
Equivalently, every left-kernel vector omega of the rigidity matrix must satisfy

    sum_e omega_e c_e = 0.

Changing the arbitrary representatives X_v changes c by an element of the
matrix image, so the obstruction class in its cokernel is independent of them.
The constant stress on the straight p-cycle is the certificate in section 3.

This also makes the minimum section defect an exact finite linear-code problem.
For a fixed residue framework, let C be the image of its rigidity matrix in
k^E. Every section has edge-error vector c+R h, and every h is allowed.
Consequently the smallest possible number of failed unit edges is exactly

    min_(z in c+C) #{e : z_e != 0}.                   (C)

The stress space is C's dual parity-check space; its pairing with c is the
fixed carry syndrome. Thus this connection preserves the actual edge-loss
objective, rather than merely comparing two uses of linear algebra. The
edge-disjoint cycle stresses in the full graph G_(q,1) give the proved lower
bound |E|/p. This bound is not asserted for arbitrary subframeworks, such as
forests whose defect is zero. Equality is certified at q=3,5. A concrete finite open target is whether equality holds
at q=29, or whether an additional stress certificate forces more than 435
failed edges. No value beyond the proved lower bound is asserted there.

At q=3, the full rigidity matrix has rank 12, so its six-dimensional stress
space is exactly the span of the six disjoint triangle-cycle stresses. The
code therefore splits into six length-three single-parity-check blocks; one
error per block is attainable. At q=5, rank 47 leaves a 28-dimensional stress
space, with thirteen independent constraints beyond the fifteen cycle checks.
Attaining the cycle bound is then genuinely a compatibility problem, solved
by the supplied witness. These exact ranks are checked by the section
verifier. The code here is defined by its rigidity matrix; no identification
with a particular named error-correcting code is implied.

If the rigidity rows are linearly independent, every carry is solvable and the
framework lifts to this precision. More strongly, a finite framework with
full row rank lifts to an exact local unit-distance framework by multivariable
Hensel; distinct residue vertices remain distinct. Forests give elementary
positive controls by choosing and lifting their edge steps successively.

### A necessary correction: stress-free lifting stops at four colors

The full-row-rank sufficient condition is attractive, but cannot supply a
five-chromatic witness. The following elementary barrier explains why.

**Theorem 4 (stress-free ceiling).** At odd residue characteristic, a finite
unit-distance framework with distinct vertices and linearly independent
rigidity rows has a 3-degenerate, hence 4-colorable, underlying graph.

For any subframework on v>=2 vertices, translations give two independent
vectors in the rigidity kernel over k. Choose a nonzero delta in ell with
c(delta)=-delta. The infinitesimal rotation h_v=delta x_v is a third:

    Tr(c(a_e) delta a_e) = N(a_e) Tr(delta) = 0.

It is independent of translations because at least two positions are distinct.
Thus every such subframework has rank at most 2v-3. Independence of all rows
implies independence of each subset of rows, so every subgraph has m<=2v-3
edges. Every nonempty subgraph therefore has a vertex of degree at most 3.
Successively removing such vertices and coloring in reverse order proves
the assertion. QED.

The ceiling is sharp. In ell=F_17[theta]/(theta^2-3), the seven vertices

    (0,0), (1,0), (9,2), (10,2), (15,1), (5,5), (3,6)

form a Moser spindle: two unit diamonds sharing vertex 0, with their tips
joined by a unit edge. There are exactly eleven unit edges, and their
rigidity matrix has rank eleven. In a three-coloring, each diamond forces
its tip to share the central vertex's color, contradicting the tip edge.
An explicit four-coloring is (0,1,2,0,1,2,3). The standard-library verifier
checks all edges, all 3^7 labeled three-color assignments, the exact rank,
and actual carry lifts modulo 17^2 and 17^3.

There is also a useful quantitative version. For k>=2, if a k-chromatic residue graph
contains a k-critical subgraph on n vertices, that subgraph has minimum degree
at least k-1. Its rigidity stress space therefore has dimension at least

    m-(2n-3) >= ceil((k-5)n/2)+3.                    (S)

In particular, every five-chromatic critical framework has at least three
independent stresses. A six-chromatic critical framework has at least
ceil(n/2)+3. These are dependencies that the arithmetic carry must satisfy;
the bound does not by itself prove any particular carry incompatible.

**Precise next problem.** At q=29, test chromatically critical subgraphs by
computing their forced stress spaces and the first carry class in each
cokernel. For those with zero first obstruction, test the subsequent carry
classes and identify a mechanism that solves them at every depth. Full-row-rank
Hensel lifting alone cannot answer this question for five or more colors.
A successful compatible lift would still be a local statement: realization
in the prescribed global number field, the chosen real embedding, and any
global chromatic consequence require separate arguments. No such
high-chromatic lifting witness is supplied here.

## 5. Arbitrarily deep poles are invisible to fixed-place finite precision

Let I=product P_i^(n_i) involve finitely many c-fixed primes of L. Let Q be
any prime with c(Q) != Q. Fix M>=1. Congruence modulo I here means
congruence in the semilocal ring, equivalently v_(P_i)(u-1)>=n_i for
every i; u may have poles elsewhere and need not belong to O_L.

**Theorem 3.** There exists u in U(L) such that

    u = 1 mod I,
    v_Q(u)=M,   v_(c(Q))(u)=-M.

### Proof

Use CRT to choose an algebraic integer t satisfying t=1 mod I,
t in Q^M but not Q^(M+1), and t=1 mod c(Q). The nonzero class in
Q^M/Q^(M+1) exists because this is a nonzero residue-field vector space.
The three sets of local conditions are disjoint. Since c(I)=I,

    u=t/c(t)

has norm one, is 1 modulo I, and has the stated valuations. Other split-prime
poles may occur; no claim of support at only this pair is made. QED.

Thus every finite collection of fixed-place residue digits identifies
infinitely many unit directions with unbounded split-prime poles. This
kernel is the denominator coordinate that the quotient loses. In the
current degree-48 field, Q may be one of the sixteen **conjugation-swapped
places above 29 itself** described in section 2. Even knowing arbitrarily
many but finitely many digits at all eight fixed 29-places does not bound
the pole order at those other 29-places.

An explicit control is L=Q(i), alpha=(3+4i)/5=(2+i)/(2-i). For r>=1,
alpha^(4*3^(r-1))=1 modulo 3^r, while its two valuations above 5 have
opposite magnitudes 4*3^(r-1). Powers multiply that magnitude without
changing the observed residue. The general proof does not depend on this
example or the exact minimal order of alpha modulo 3^r.

## 6. CM-integral directions: a clean denominator baseline

Assume now that L is a CM field, so all embeddings send complex conjugation
to ordinary complex conjugation. If u in O_L has u c(u)=1, every conjugate
of u has absolute value 1. Kronecker's theorem gives

    O_L intersect U(L) = mu(L),

the finite group of roots of unity. Accordingly the unit-distance graph on
O_L has components that are translates of Z[mu(L)]. Its chromatic number
is 2 when |mu(L)| is a power of 2 and 3 otherwise. The cyclotomic coloring
classification is CITED to Davidson, arXiv:2610.09474, Theorem 3.3, whose
prime-power zero-sum-label proof was source-checked here.

In particular, any finite graph realized by integral CM coordinates is
3-colorable. This is an arithmetic corollary, not a new proof of Davidson's
rational-angle theorem. It fails as an integrality-to-root-of-unity argument
outside CM fields: the nonreal unit-circle root of

    X^4-X^3-X^2-X+1

is an algebraic-integer unit that is not a root of unity. Its other conjugates
include real reciprocal roots off the unit circle.

More invariantly, for unit directions u,v in a CM field,

    v_P(u)=v_P(v) at every finite prime P
        iff u/v is a root of unity.

Thus the finite-prime divisor records exactly the direction modulo rational
angle. A common unit rotation by an element of U(L) translates every divisor
signature by the same amount. The widths

    max_edges v_P(u_e) - min_edges v_P(u_e)

are rotation-invariant. If every width is zero, the graph is 3-colorable.
Nonzero width is necessary for a 4-chromatic graph, but emphatically not
sufficient: all of Q(i)'s field plane is bipartite, although its directions
have arbitrarily large pole widths above 5.

For a finite denominator budget, let D be a nonzero integral ideal and
consider norm-one u in D^(-1). At a split pair {Q,c(Q)}, write
b_Q=v_Q(D) and b_cQ=v_cQ(D). The exponent v_Q(u) lies in
[-b_Q,b_cQ], and every c-fixed valuation is zero. A fixed complete valuation
signature admits at most |mu(L)| directions, because its ratios are roots
of unity. Therefore

    #{u in D^(-1):u c(u)=1}
      <= |mu(L)| product_(split pairs) (b_Q+b_cQ+1).   (D)

Only pairs meeting D need appear. This is an upper bound: ideal-class and
unit-norm conditions can prevent signatures from being attained. It is
sharp for L=Q(i), D=(5^m): the directions are exactly
mu_4 times ((2+i)/(2-i))^j, -m<=j<=m, giving 4(2m+1).

## 7. What changed, and what remains open

The two relevant operations are now sharply distinguished:

- More copies of the same local observation: no chromatic improvement beyond
  the deepest factor, by Theorem 1.
- More precision in one local observation: no section preserving all residue
  geometry, by Theorem 2; the effect on chromatic number remains open.
- Bounded local observations: no bound on split-prime pole depth, by Theorem 3.
- Stress-free local lifting: at most four colors, by Theorem 4. Higher
  chromatic targets necessarily bring many carry-compatibility equations.

The connection contract is explicit. Source: exact norm-one steps and their
finite frameworks. Target: finite local graphs. Map: simultaneous reduction.
Preserved predicate: unit adjacency. Destroyed data: the mixed-characteristic
carry and split-prime divisor. Sidecars: stress compatibility and valuations.
Cheapest decisive tests: the straight p-cycle, full simultaneous norm-one
lifting, and the t/c(t) kernel construction. No LRC, Collatz, or Jacobian
claim is inferred from this dictionary.

The local 29-adic chromatic sequence is nonincreasing in depth and takes
integer values, so it eventually stabilizes. This tautology supplies neither
a stabilization depth nor its terminal value. It does not turn finite testing
into an all-depth decision procedure.

## 8. Reproduction and external-status notes

Run `python3 reproduce.py`, `python3 verify_section_witnesses.py`, and
`python3 verify_rigidity_barrier.py`. All three use only the Python standard
library.
`exact-results.json` and `reproduce.out` are matching saved outputs.

- Exact rank/stress certificates: p=3,5,7,11,13,17,19,23,29,31,43,53,101.
- Independent complete step-set dynamic programming: p=3,5,7,11. After p
  steps reducing to 1, the first coordinate is p modulo p^2, never zero.
- Degree-divisible-by-p hostile: F(x)=x^p at p=3,5,7,11; the cycle does lift.
- Finite tensor diagonal/projection checks: q=3,5,7. The proof, not these
  small examples, establishes Theorem 1 at q=29 and all depths.
- Sharp edge-loss witnesses at q=3,5 are checked directly modulo p^2,
  independently of the heuristic search and rigidity-matrix implementation.
  Their lower bounds come from the proved cycle partition, not a solver.
- The sharp stress-free ceiling is checked by `verify_rigidity_barrier.py`: a
  seven-vertex, eleven-edge Moser spindle over F_(17^2) has rank eleven, needs
  four colors, and lifts through two successive carries.
- The sharp section witnesses were discovered by finite search, but their
  existence is certified by the displayed coordinates and direct integer norms.
  Their optimality follows from the proved cycle lower bound. Discovery-search
  code and incomplete searches are not needed in the publication package.

Current global status was checked against primary sources on 2026-10-19.
OpenAI's September 23, 2026 preprint **reports** that arbitrary five-colorings
fail, with a Lean formalization; together with the seven-coloring it states
6<=chi(R^2)<=7. The current catalogue and source implementation were inspected.
The source file `OAI/Geometry/PlaneColoring/Five.lean` contains
`no_proper_five_coloring`; its comparator configuration permits only the three
standard axioms propext, Quot.sound, and Classical.choice. This session did
**not** rebuild the imported dependency tree or independently audit that proof.
The comparator challenge files are theorem specifications with intentional
`sorry` placeholders, not the implementation modules. None of this note's
new elementary statements depends on the claimed global improvement.

Primary links:

- OpenAI preprint and date: https://github.com/openai/math/tree/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026
- Formalization scope: https://github.com/openai/math/blob/main/lean/docs/158.md
- Actual lower-bound implementation: https://github.com/openai/math/blob/main/lean/OAI/Geometry/PlaneColoring/Five.lean
- Comparator/axiom configuration: https://github.com/openai/math/blob/main/lean/ComparatorChallenges/EuclideanFiveColor.json
- Current corrections page: https://github.com/openai/math/blob/main/history.md (no family-158 correction listed when checked)
- Davidson, *Coloring Monomial Cayley Graphs: The One-Parameter Case*, submitted October 7, 2026: https://arxiv.org/abs/2610.09474
- Davidson's theorem/proof: https://arxiv.org/html/2610.09474v1#S3.SS2
- Madore, *The Hadwiger–Nelson problem over certain fields*: https://arxiv.org/abs/1509.07023
- External local-coloring experiments: https://github.com/decalion89/chromatic-number-of-the-plane/blob/main/notes/local_colourings.md

Novelty boundary: the exact equal-local-family consequence, the universal
carry obstruction with its homogeneous extension, and the denominator-kernel
application were developed in this session. Their ingredients are elementary;
no exhaustive novelty or publication-priority claim is made. The CRT lifting,
finite-field reduction technique, general rigidity linearization, Kronecker
argument, and rational-angle coloring have prior literature or repository
antecedents as identified above.
