# A finite forbidden basis exists and can be synthesized exactly

**PROVED existence and terminating algorithm; NOT a computed list or an explicit numerical order bound.** Depends on `EXTENSION4.md` and the hereditary-slack observation in `EQUALITY_SKELETON.md`. This classifies the *type of decision problem* left by equality; it does not prove Chan–Pak.

## 1. The theorem

Let E be the class of finite rooted tournaments (T,z) satisfying

    H(T)−H(T−z)=d(G_z)−1.

Then E has a **finite forbidden rooted-induced-subtournament characterization**. Moreover, a finite complete forbidden family, and then its inclusion-minimal rooted-isomorphism types, can be computed by a terminating exact algorithm.

This theorem supplies neither the cardinality of that family nor a maximum forbidden order. In particular it does not claim that patterns through order seven suffice. The algorithm below has not been run to completion in this session.

## 2. Reduction to one fixed triangle and a transitive spine

Write S(T,z)=H(T)−H(T−z)−d(G_z)+1. It is nonnegative and nondecreasing under rooted induced extension.

If G_z is empty, S=0, as explained in the companion note. Otherwise choose any G_z edge ab, with a∈A,b∈B. The three vertices form the fixed directed triangle

    z→a→b→z.

If U−{a,b} is not transitive, it contains a directed triangle on three additional vertices. Together with z,a,b these give a rooted six-vertex induced tournament with S≥4, since the ab fibre has a cyclic complement. There are only finitely many possible rooted six-vertex types of this form; let F₆ be the family of all of them. Each is a genuine positive-slack witness.

If U−{a,b} is transitive, order its vertices in their unique transitive order. Encode each by its three orientations relative to z,a,b, giving an alphabet Σ of size eight. A word w∈Σ* defines one literal rooted tournament T(w): the core is the fixed triangle, the spine is transitive in word order, and the letters give every core-spine edge.

Deleting letters gives a rooted induced subtournament retaining the marked core. Consequently the language

    L₊={w:S(T(w),z)>0}

is upward closed under the subsequence order. This observation is where rooted heredity becomes a finite-language tool.

## 3. Finite-alphabet subsequences have no infinite bad sequence

We include the needed combinatorial lemma rather than relying on an unexamined external classification.

**Lemma.** For a finite alphabet, every infinite sequence of finite words has an earlier word that is a subsequence of a later word.

**Proof.** Suppose there is an infinite bad sequence. Construct one by taking, at each index, a word of minimum length among those permitting an infinite bad continuation after the already fixed prefix. No word is empty. Infinitely many selected words have the same final letter, since the alphabet is finite. List their indices as i₀<i₁<⋯ and write w_{i_j}=v_j a.

The sequence

    w₀,…,w_{i₀−1},v₀,v₁,v₂,…

is also bad. A prefix word embedding in some v_j would embed in the original w_{i_j}. An embedding v_j into v_k would extend by the common last letter to an embedding w_{i_j} into w_{i_k}. Both are forbidden by the original bad sequence. But v₀ is shorter than w_{i₀}, contradicting the minimal choice at index i₀. This proves the lemma.

It follows that any upward-closed word language has finitely many minimal members. Apply this to L₊. Their rooted tournaments, together with F₆, give a finite forbidden family for E: any positive-slack tournament has either the six-vertex witness, or a fixed-triangle/spine encoding containing a minimal positive word. Conversely every listed pattern has positive slack and forces positive slack in any rooted induced extension.

This already proves finite-basis existence. The rest establishes effectivity rather than pretending the lemma gives a useful size bound.

## 4. Exact weighted automata for the three counts

A finite weighted automaton is used here only as an explicit integer linear representation: every input letter acts by an integer matrix, and the output is an initial row times the matrix product times a final column.

### Hamilton paths on a fixed core plus a transitive spine

For a core of r vertices, fix one of its r! linear orders as the order in which a Hamilton path visits core vertices. The spine vertices are assigned to r+1 bins: before the first core vertex, between successive core vertices, and after the last. Within each nonempty bin they necessarily occur in increasing spine order.

Scan the input word. Assign the current vertex to any bin, subject to the following finite checks:

- on the first insertion into a bin with a preceding core vertex, require that core vertex to point to the new vertex;
- record whether the latest vertex in that bin points to its following core vertex, if any;
- an empty bin between two core vertices is legal at acceptance exactly when the corresponding core edge points forward.

Each bin needs at most three states: empty, occupied with legal final endpoint, occupied with illegal final endpoint. At acceptance every occupied bin must have legal final endpoint; exterior empty bins are allowed. This counts every Hamilton path once, since the path uniquely determines its core permutation and bin assignment.

Thus H(T(w)) has a positive integer automaton with at most 3!·3⁴=486 states. H(T(w)−z), using core a→b and the same spine word, has at most 2!·3³=54 states. Unused letter coordinates can simply be ignored by this second automaton.

### Dominating sets of the marked bipartite graph

Choose first the selected subset C⊆{a,b}, one of four possibilities. Process spine vertices in order, branching on selected versus unselected. Retain four Boolean flags:

1. whether an earlier spine A vertex was selected;
2. whether some earlier unselected spine A still needs a later selected spine B;
3. whether core vertex a is selected or dominated;
4. whether core vertex b is selected or dominated.

Core flags are initialized from C, including domination across the edge ab. A selected spine A sets flag1 and marks b dominated if it is adjacent to b. A selected spine B clears flag2 and marks a dominated if adjacent to a.

For an unselected spine A, a selected adjacent core b immediately dominates it; otherwise set flag2. For an unselected spine B, require either flag1 or a selected adjacent core a; otherwise reject this branch.

These rules use the precise spine incidence: a spine A is adjacent to a spine B exactly when the A occurs earlier. Thus all pending A vertices are simultaneously dominated by the next selected B, and every unselected B must be dominated at arrival. At the end require flag2 false and both core flags true.

The automaton has at most 4·2⁴=64 states. It counts every dominating set once. The four flags are genuine boundary data; dropping pending domination or core coverage would invalidate the calculation.

### The slack representation

Take the direct sum of the two Hamilton automata, the domination automaton and a one-state constant-one automaton, with final signs +1,−1,−1,+1. This yields an exact integer linear representation of S(T(w),z), of dimension at most

    486+54+64+1=605.                                    (3)

The signed representation is not asserted to be entrywise nonnegative. Its output is nonnegative by the universal theorem. The distinction matters for the next decision step.

## 5. Deciding whether an unlisted obstruction exists

Suppose finitely many positive words B have already been found. The language of words that avoid every member of B as a subsequence has an explicit finite deterministic automaton: for each b∈B, track the longest prefix of b greedily embedded so far, rejecting when all of b has embedded. The product of these trackers suffices.

Take the product of this avoidance automaton with the 605-dimensional slack representation. Its integer output is S(T(w),z) for avoiding words and zero for rejected words.

Whether this rational series is identically zero is decidable by exact finite-dimensional linear algebra. Start with the initial row, repeatedly multiply retained rows by each letter matrix, and keep a basis of the reachable row span. Test the final column on the basis rows. If every pairing is zero, the series vanishes on every word. If some pairing is nonzero, the row's recorded word is a concrete witness. At most the representation dimension many independent rows are retained, so the process terminates. Exact rational elimination is sufficient; no numerical tolerance is used.

The resulting witness has positive slack because slack is globally nonnegative. Greedily delete letters whenever slack stays positive, until no deletion is possible. The result is minimal under subsequence: any smaller positive subsequence would, by upward closure, allow at least one such single deletion.

Add the minimized word to B and repeat the equivalence test. This whole synthesis loop terminates. Otherwise the successive witnesses would form an infinite bad sequence: each newly found word avoids every earlier one as a subsequence. Section3 rules this out.

At termination, B is a finite complete minimal-word basis for L₊. Convert it to rooted tournaments, add F₆, and remove rooted-isomorphic duplicates and any pattern containing another listed one. These are finite induced-subgraph checks. The result is a finite complete minimal rooted forbidden family.

## 6. Scope, failure boundaries, and concrete next step

This is stronger than a bounded search proposal: it proves that an exact obstruction-synthesis computation must eventually terminate. It is weaker than an implemented classification: the complete basis and its maximal order are not given.

The fixed triangle is essential. Arbitrary tournaments are not well-quasi-ordered under induced subtournaments, so one cannot apply the word lemma before proving the edge-feedback/transitive-spine reduction. Likewise, the weighted representation alone does not make its nonzero support regular; subsequence monotonicity supplies the special structure.

A practical next step is to implement the 605-state-or-smaller representation and first validate its three count outputs on short words by direct DP. Then begin synthesis with the explicit rooted order-five patterns and the six-vertex feedback obstructions, reducing the avoidance automaton before each exact linear test. That computation is not silently represented as completed here.

## 7. Focused automaton controls

`equality_automata_probe.py` implements the finite-state transitions directly as integer state-mass maps. It compares all three outputs against literal Hamilton subset DP and direct domination enumeration for all eight-letter words of lengths0..3, plus100 seeded words of lengths4..7: **685 cases**, all passed. It also checks the explicit necessary-skeleton counterexample has slack2. The receipt is `equality_automata_receipt.json`.

These controls validate the new finite-state interfaces only. The full avoidance-product linear algebra and forbidden-basis synthesis have **not** been run. Run the checker from the repository root with `python3 05-knowledge/results/session261_creative_chanpak_20261111/bridges/equality_automata_probe.py`.
