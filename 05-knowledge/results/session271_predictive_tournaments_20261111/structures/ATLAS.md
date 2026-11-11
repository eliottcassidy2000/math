# Structural and interaction atlas

**Status legend:** PROVED = theorem or explicitly inherited theorem in scope;
FINITE-EXACT = bounded computation; OPEN = classification frontier;
DESIGN = specified implementable object, not a claimed complexity theorem.

## 1. The five coordinates of a classification claim

Every row below declares a **carrier**, **allowed constructors**, **equivalence**,
**observer**, and **scale**. Leaving any one implicit changes the question.
The normal form classifies finite graphs under ordinary isomorphism. The same
graphs support much smaller equivalence classes for H-contexts, spectra, or
linear dynamics. No one observer deserves the name 'the complete structure'.

### A. Ordered strong-component chains

- Carrier: the unique sequence of SCCs, singletons allowed.
- Constructors: ordinal sum, with adjacent linear nodes flattened.
- Exact ordinary object: ordered sequence of actual strong isomorphism types.
- Laws: H multiplies, adjacency χ multiplies, Aut is a direct product. The full
  q product is the ordered-cover contraction from Session268 modular (6).
- Lost by these commutative observers: SCC chronology, often chirality.
  K1 followed by C3 and its reverse have equal H, q, χ, and |Aut| but different
  graph types. Retain the ordered tree, or use a targeted oriented sidecar.
- Concrete niche bridge: the skew deletion-response current of
  [THM-3369](../../../../01-canon/theorems/THM-3369-skew-deletion-response-and-ordered-join-orientation-current.md)
  contracts an intrinsic skew adjacency with deletion responses and detects
  the lost order bit. It is not a full chronology reconstruction theorem.
- Regime: multiplicative target constructors and symbolic words. Exact finite
  alphabets can have rational scalar generating functions in this restricted
  chain language. This does not imply unrestricted A000568 is rational.

### B. Finite-prime-alphabet substitution trees

- Carrier: canonical linear/prime decorated tree, fixed finite prime alphabet.
- Constructors: arbitrary-depth substitution and ordered SCC assembly.
- Exact object classification: child isomorphism types modulo the actual prime
  quotient automorphism group. Ordinary, self-converse-identified, and labelled
  counts remain different.
- Laws: H consumes full q; χ consumes only χ via the polynomial quotient
  determinant; Aut consumes full colored quotient action and child groups.
- Examples: C3-only family of Session268; uniform C3 towers; mixed C3/Paley7
  factor words. The C3-only family is much larger than uniform C3 towers.
- Regime: exponential class growth at most K^n; rare globally, rich locally.
  Iteration depth can be huge without approaching global prime diversity.
- Interaction: a prime-node decoration is a multi-input operation, not an
  independent product. Observable symmetry can exceed object symmetry.

### C. Infinite-prime-alphabet trees, including all finite tournaments

- Carrier: the same normal form, but allow every finite prime quotient.
- Completeness: PROVED conditional only on the inherited decomposition theorem.
  The prime labels themselves carry unrestricted orientation data.
- Obstruction to a finite catalogue: labelled random tournaments are prime with
  probability 1−O(n²2^−n), by the explicit bound in PROOF.md. The depth-one
  prime layer already carries essentially all labelled entropy.
- Consequence: small decomposition depth is NOT a generic tractability test.
  The important complexity parameters are prime quotient size, port width,
  observer-state dimension, and integer precision as well as tree depth.
- Forecast bridge: the global class census is close to labelled mass divided
  by n!, but the rare correction is driven by genuine automorphisms, not twin
  modules. The session forecast/ensemble transport separates these scales.

### D. Regular-odd homogeneous towers

- Carrier: uniform products Q[B] with regular odd B; mixed regular factor words.
- Exact F2 state: global degree parity plus internal elementary divisors.
- Inherited law: A_(Q[B]) is similar over F2 to
  (A_Q+d_B I) ⊕ (A_B|1-perp) repeated |Q| times.
- Regime: the C3/Paley7 alphabet yields only fixed points and eventual
  3-cycles, with exact multiplicities for every word. Exponentially large
  adjacency matrices admit short word-level dynamics descriptions.
- Interactions: child degree parity shifts the ancestor modes; order of factors
  changes multiplicities. This is not a commutative product of invariants.
- Boundary: regularity and odd child order supply an invariant direct split.
  Arbitrary irregular children still admit our χ compiler, but not this
  similarity theorem. χ alone misses Jordan data and transient depths.
- Source: [Session268 dynamics](../../session268_composable_tournaments_20261111/dynamics/PROOF.md).

### E. Finite-port seams and boundary grammars

- Carrier: tournaments with named distinct ports; graph plus port tuple is the
  object BEFORE marks are forgotten.
- Constructors: uniform substitution plus explicit cross-arc overrides.
- Exact spectral state: χ and a matrix of resolvent entries on all-ones/port
  columns; Woodbury/determinant wiring closes.
- Exact Hamiltonian state: ordered path-cover tensor with start/end port labels.
  One seam contracts to the endpoint and two-path-cover law of Session268.
- Exact graph state: retain the actual graph or recompute canonical modules.
  Old child groups cannot automatically be multiplied after a seam edit.
- Regime: local hidden state becomes visible, and mixed interactions appear.
  No fixed number of ports implies a fixed finite Hamiltonian state without
  additional run restrictions. No bounded number of seams repairs the generic
  entropy deficit of a finite prime alphabet.
- Sharp boundary: a typical tournament has normalized arc-edit distance 1/2
  from every fixed finite-prime-alphabet family, even under relabelling.

### F. Ordered-spine and finite-core word languages

- Carrier: fixed k-vertex marked core plus an ordered transitive spine;
  letters choose each spine vertex's orientations to the core.
- Equivalence: words/marked graphs first; different words can be the same
  unmarked tournament. Retain marking weights when changing the census.
- Exact observer: a finite rational linear representation for H for fixed k;
  Session270's particular two-core series has minimal rational dimension 17.
- Regime: fast target compiler, regular-language restrictions, rational
  generating functions after summing words of each length. Mean/mass does not
  determine the set of attained integer values.
- Interactions: matrix products depend on letter order; a scalar output does
  not iterate by itself. Reachable/observable reduction preserves the full
  series, not positivity or integer-state realizability at minimal dimension.
- OPEN: a decreasing target algorithm or returning cone covering all eligible
  odd integers. The million-target atlas is finite evidence, not the proof.
- Source: [Session270](../../session270_chanpak_enumeration_20261111/SYNTHESIS.md).

### G. Interacting cell words and imposed dynamics

- Carrier: intrinsically decoded SCC words of seam-coupled cells, with cyclic
  neighbor relation on their ordered positions.
- Exact static enumeration: Session268 has 2^N rigid ordinary classes but only
  floor(N²/4)+2 H values; transfer sum trace(W^N) is C-finite in spatial size.
- Dynamics: a specified cellular rule updates the decoded word; a local bit
  change alters two actual arcs. This is a genuine implemented action.
- Boundary: the update rule is additional structure, not an automatic property
  of an arc-flip graph. Finite rings are finite; observed Rule110 behavior is
  not a proof of computational universality or an all-time H recurrence.
- Source: [Session268 automata](../../session268_composable_tournaments_20261111/automata/PROOF.md).

## 2. Operator interaction laws and what each one destroys

**Uniform substitution ↔ normalization.** Associativity holds with ports; normal
form quotients redundant syntax. Prime child isomorphism orbits must be retained
for class counts. Replacing actual child classes by equal q or equal χ is valid
only for the corresponding observer, not graph classification.

**Converse ↔ substitution.** It reverses the quotient and every child, and
reverses a linear sequence. H, q, χ and |Aut| are invariant. The canonical tree
with orientation and child chirality restores the lost distinction.

**Seam ↔ substitution.** A seam is a marked finite-rank adjacency update and a
transition override in the path-run compiler. It can destroy the old modules.
The mark tuple is indispensable, as the TT2 hostile proves. Applying an edit
before/after substitution requires explicitly lifting its endpoint addresses;
there is no unqualified commutation law.

**Cut switching ↔ spectrum.** For centered C=2A−J and signs D, switching gives
C'=DCD; centered spectrum survives. Recovering ordinary adjacency introduces
the changed signed all-ones observer. THM-3315 supplies the sidecar. H is not
preserved merely because centered eigenvalues are. Our ordinary χ law and
that centered switching law concern different spectral quotients.

**Reduction modulo p ↔ substitution/seams.** Polynomial determinant laws commute
with reduction. Rational formulas must be interpreted symbolically over F_p(z),
not by dividing by an integer that becomes zero. Characteristic polynomial is
not rational canonical form. Store invariant factors for full finite-field
linear dynamics.

**Averaging ↔ class quotient.** Sum over labelled realizations, ordinary classes,
marked words, or response buckets only with the correct weights. Aut stabilizers
and reciprocal mass are the bridge. H-value support cannot be read from moments
alone. The parent session's ensemble transport quantifies when labelled and
ordinary expectations are close; it does not erase markings for free.

## 3. Concrete bridges, not analogies

1. **Modular tree -> polynomial signal network.** Map each child to (χ,G) and
   each quotient/seam to K. Preserve global χ and boundary resolvent by exact
   determinant/Woodbury laws. Lose actual graph and Hamiltonian covers. Test:
   irregular-child substitutions and the two-root seam in receipt.json.
2. **Finite grammar -> coding-theoretic edit ball.** Map n-vertex orientations
   to N-bit words. A family with K^n unlabelled bases has at most n!K^n labelled
   centers. Binary entropy bounds the union of Hamming balls. Preserve actual
   number of arc reversals, lose graph labeling after minimization. The result
   is a proved half-density edit barrier, not just a metaphor for complexity.
3. **Canonical automorphisms -> group action/cycle indices.** Map repeated child
   isotypes to a colored quotient stabilizer. Preserve the automorphism order
   and exact orbit census. Lose nothing for this task if child group actions
   are retained. Scalar response colors alone give a strictly coarser action.
4. **Boundary tensor -> weighted language.** Map maximal module runs to words
   with endpoint labels, then contract ordered path-cover weights. Preserve H
   exactly. Lose graph isomorphism. Finite-state collapse needs a bounded-run or
   fixed-core theorem; it does not follow from the word notation.

## 4. Live concept board after these results

- Canonical trees: complete object carrier; bottleneck is the unbounded prime
  alphabet, not merely recursion depth.
- Hamiltonian response: contextual completeness under uniform modules is
  inherited; root refinement opens seams but target inversion stays open.
- Spectral response: χ alone unexpectedly suffices for uniform modules;
  matrix resolvents repair the first named-seam failure.
- Symmetry: actual colored quotient action is retained; observer symmetries
  and twin modules must not be mistaken for graph automorphisms.
- Typicality: prime entropy and edit-distance barrier explain why the richest
  controllable grammars need not dominate A000568.
- Dynamics: similarity invariants and supplied update laws classify specified
  actions; there is no operation-independent finite behavior catalogue.

Next concrete question: combine the ordered-spine 17-state H representation
with the polynomial spectral observer and port-resolvent interface, then seek
returning cones whose moves both shrink a target and remain in a realizable
boundary-state domain. This is a precise OPEN program, not a solved recursive
classification of the Chan–Pak spectrum.
