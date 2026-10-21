# Independent audit: canonical symmetry-slack compiler

Date: 2026-10-21.

Status: **PASS, scoped mathematical audit.** I found no proof defect in the proposed all-fixed-symmetry-slack rationality theorem, its self-converse/converse-merged refinements, or its stated denominator bounds. This is an independent-agent proof review and independent focused exact computation, not a formal proof or an all-tournament census at unbounded order. No repository files were changed and no Git operations were performed.

## Sources and scope

- Candidate: `session23-connections/SYMMETRY_SLACK_RATIONALITY.md`, §§1–7.
- Inherited dependency: `math/05-knowledge/results/continuation_width_symmetry_20261020/README.md`, Theorems 3–6 and Corollary 5C. The load-bearing inherited theorem for the new proof is Theorem 5's invariant set B, with |B| ≤ 2δ, whose remaining nonfixed vertices are homogeneous directed triangles that are actual full-group orbits.
- I checked the relevant inherited proof, rather than treating the bounded set or the orbit condition as a black-box assertion. The unrelated analytic/Collatz results in the inherited packet are outside this audit.
- The targeted current mistakes-ledger search did not locate a correction to the inherited symmetry-slack statements. This audit makes no historical-priority claim.

## 1. Proof obligations

### Canonical partition: valid

For overlapping, nonnested modules A and B, the three nonempty differences X=A\B, Y=A∩B, Z=B\A are forced into one of the two orders X→Y→Z or Z→Y→X. Uniformity across each pair follows by using the module condition first for a vertex of X outside B, then for a vertex of Z outside A, and varying both witnesses. Thus this is not an unjustified assumption that arbitrary modules are strong modules.

The induced subtournaments on the three pieces belong to W by heredity. Their ordered join therefore belongs to W. The union is itself a module by the common intersection witness. Every vertex has a singleton W-module; maximal W-modules consequently exist, cover all vertices, and are disjoint.

The exact maximality test also works for non-strong proposed modules. If a larger W-module intersects one proposed block, repeatedly taking its union with every intersected proposed block produces a W-module that is a union of whole blocks. Every such union is a module exactly when its quotient indices form a module. This proves both directions of the test, not merely a necessary condition.

The flag test is sufficient. A nontrivial strongly connected quotient remains strongly connected after nonempty substitution; its inflation belongs to W only when it is exactly a triangle whose three inflated blocks each have size one. Singleton SCCs allow arbitrary W substitutions. Hence internal order/content beyond the flag “exactly the one-vertex word” cannot change maximality.

### Feedback straightening: valid for arbitrary homogeneous modules

For a fixed vertex order, uniformly and independently selecting one representative from each module gives a random quotient order. For each quotient arc, the weighted backward indicator has expectation exactly the number of backward cross-arcs in the original order. Summing and taking a realization at most the mean proves the required lower bound; concatenating internally optimal module orders proves the reverse bound.

This proof does not require the quotient or modules to be transitive, and does not silently assume an optimal order was block-contiguous before the averaging step. It applies to the cyclic hostiles used in the computation.

### Actual support and the inherited bound: valid

The canonical partition makes every automorphism induce an automorphism of the fully decorated quotient, and every decorated quotient automorphism lifts. The kernel is the independent product of internal triangle rotations. A moved quotient block moves every one of its vertices. For a fixed quotient block, only its internal triangle vertices can move. Therefore the actual support is

s = 3 Σ t(w_i) + Σ_(i moved) u(w_i).

The subtraction must be moved singleton mass U, not whole-block mass V. This distinction is substantive and was checked by two separate hostiles below.

A moved maximal block cannot meet a fixed vertex outside B. Nor can it meet a quiet triangle orbit outside B: overlap closure puts the entire quiet orbit inside that maximal block, whereas moving the block sends its vertex into another disjoint maximal block, contradicting that the quiet orbit was contained in the first block. Hence every moved block lies in B and V ≤ 2δ. This proves C=(δ+U)/3 ≤ δ without replacing actual orbits by score classes or independent apparent triangle factors.

The inherited theorem's proof supports the needed orbit condition: after marking all nontriangle nontrivial orbits and every orbit incident with an external backward arc, each unmarked triangle has no external backward arc. Its uniform count of out-neighbors in every other orbit forces contiguity, hence homogeneity. The size calculation is |B| ≤ 5d/4+6e ≤ 2δ. No assumption of independently rotating the marked triangle orbits is made.

### Bounded quotient and independent long slots: valid

An optimal weighted quotient order of cost C ≤ δ has at most δ backward arcs, so its endpoint set S has size at most 2δ. Each gap outside S is a homogeneous transitive module of the quotient. Two vertices in a gap would merge their W blocks, contradicting canonical maximality. Thus |Q| ≤ 4δ+1.

A block of size >δ cannot meet a backward arc in an order of cost at most δ, since all other block sizes are positive. The long slots therefore occur outside S and induce a transitive tournament. There are at most 2δ+1 of them.

For any fixed short-word assignment, every inadmissible order costs >δ for every permitted long filling, while an admissible minimizer has constant cost C0 ≤ δ. Thus C=C0 for every filling. An automorphism preserving the long-slot set must fix it pointwise, because a finite transitive tournament has no nonidentity automorphism. This proves equality of the configuration stabilizer and the actual stabilizer for every filling. Equality of some chosen long words introduces no extra automorphisms.

### Finite compiler and class counting: valid

Every actual tournament of positive slack has a canonical quotient satisfying the finite bound and induces exactly one configuration orbit. Conversely, maximality, C0, and the actual short-word stabilizer are filling-independent and certify the exact slack of every permitted filling.

When transporting an orbit to its chosen representative, two possible transports differ by a configuration stabilizer, which fixes every long slot. Hence the filling is uniquely defined at that representative. Distinct long fillings cannot become isomorphic through a hidden quotient permutation. This is precisely the injectivity missing from the inherited exceptional grammar.

Enumeration of tournaments up to order 4δ+1 is finite and decidable, even though impractical for general δ. The theorem does not disguise unknown numerical tournament counts as inputs. It proves an exact finite compiler without asserting a tractable runtime.

### Anti-isomorphisms and palindromes: valid

An anti-isomorphism preserves the canonical partition and reverses each W-word. For a self-converse configuration orbit, all realizing anti-isomorphisms differ by a configuration automorphism, so they have the same action on long slots. That action reverses their unique transitive order, giving paired slots and at most one middle slot.

Paired fillings are w and reverse(w); a middle filling is a palindrome. The cutoff remains n(w)>δ for each member of a pair, exactly as encoded by Lδ(z²,y²). The middle contributes P>δ(z,y). No involutivity of the anti-isomorphism on short blocks is needed to produce an isomorphism to the converse.

The polynomial short contribution z^a y^b is counted once per self-converse configuration orbit, not once per anti-isomorphism. Hence there is no centralizer or anti-isomorphism multiplicity missing from the formula. Burnside for the two-element converse action then gives (Fδ+Cδ)/2.

### Denominator and asymptotic claims: valid

The exact positive configuration sum gives the asserted powers of 1−z−yz³ and 1−z²−y²z⁶. At y=1, multiplying a top-long-slot term by the maximal denominator power and evaluating at the positive root leaves a strictly positive quantity. Thus top poles cannot cancel. The positive root is the unique smallest-modulus root, and a nonempty layer admits arbitrarily long source/sink W words, ensuring at least one long slot occurs. The stated eventual recurrence and exact leading polynomial-exponential asymptotic follow.

## 2. Independent finite controls

Reproduction:

    python3 session23-slack-review/check_actual_tournaments.py

The source imports no candidate/inherited verification implementation. Output is written to `actual-tournament-receipt.json`; a captured stdout copy is `actual-tournament-run.txt`.

Exact universe:

- Every labelled tournament of orders 1–4.
- Twelve fixed-seed arbitrary tournaments at each order 5–8.
- Named positive and hostile actual tournaments, including long-word/converse channels through order 15.
- 548 proposed W-block substitutions: every labelled quotient of order 1–3, each block chosen from 1, 11, 111, 3. This includes both maximal and deliberately nonmaximal proposed partitions.
- All 27 cyclic quotient weight profiles in {1,2,3}³, plus nontransitive-module substitutions over cyclic/strong quotient controls.

Methods:

- Feedback by subset dynamic programming over the last actual vertex, independent of representative sampling.
- Full actual vertex automorphism and anti-isomorphism sets by adjacency-consistent extension with degree filtering.
- Modules enumerated as actual vertex subsets and checked against every outside vertex.
- W membership by actual strongly connected components.
- Exact comparison with the quotient's weighted DP, stabilizers, short/long signature, proposed maximality flags, and palindrome test.

Result: **139 actual tournaments; 548 proposed partitions; 4,576 assertions; all passed.** These are exact bounded controls, not a claim that the finite compiler has been executed for an entire positive slack layer.

## 3. Particularly informative hostiles

### Diagonal C3 synchronization at δ=6

Use two identically directed triangles a_i and b_i, with a_i→b_j iff i=j. Direct results:

- f=4, s=6, δ=6, |Aut|=3, self-converse.
- Neither original triangle is homogeneous.
- The canonical W partition is **three transitive pairs**, not six singletons: {a1,b0}, {a2,b1}, {a0,b2} in the chosen indexing.
- The canonical decorated quotient is C3[11,11,11].
- C_Q=4, U=V=6, hence δ=12−6=6.

This is a useful explicit certificate for the new compiler. An initial *test expectation* that the canonical quotient would have six singleton blocks failed and was corrected after direct module inspection; the mathematical candidate did not make that mistaken claim.

### Moved whole-block mass is not moved singleton mass

For C3[C3,C3,C3]:

- n=9, f=12, s=9, δ=27, |Aut|=81.
- Canonical words are (3,3,3), weighted quotient cost C=9.
- U=0 but V=9. Correct: 3C−U=27. Incorrect whole-mass subtraction would give 18.

For C3[13,13,13]:

- n=12, f=19, s=12, δ=45, |Aut|=81.
- C=16, U=3, V=12.
- It is not self-converse; replacing the actual short-word decoration by mere lengths would incorrectly lose this distinction.

### Converse pairing versus arbitrary equal long words

For the δ=3 strong core C3[1,11,1], adding long exterior words 13 on the left and 31 on the right gives a self-converse tournament. Adding 13 on both sides does not. The direct actual anti-isomorphism test agrees with reversal pairing, not equality pairing.

For C3[1,131,1] with equal transitive length-four words on both ends, the three long slots satisfy the paired-tail/palindromic-middle conditions; the resulting order-15 tournament is self-converse. One-sided long tails are correctly rejected.

## 4. Optional exposition improvements

No repair is required for correctness. Useful additions would be:

1. Spell out the synchronized hostile's canonical C3[11,11,11] form, rather than leaving the canonical quotient implicit.
2. Include one explicit U≠V hostile to make the singleton-mass correction impossible to miss.
3. In the injectivity paragraph, state explicitly that two transports to the chosen configuration representative differ by a stabilizer fixing each long slot.

Scope boundaries remain essential: rationality is for exact symmetry slack, not fixed fundamental black degree; the general compiler is not claimed efficient; no complete δ=6 series is supplied; no marked Hamilton-path/port or metagraph adjacency conclusion follows from scalar class counting.


## 5. Final publication binding

Final review performed 2026-10-21 after inspecting `/tmp/session23_slack_final_exposition.diff`, the final proof, and the frozen publication allowlist. The status promotion, explicit U≠V hostile, and representative-transport explanation introduce no mathematical defect. The added low-slack bivariate/converse formulas in §9 also follow from the inherited normal forms as stated. The additional order-six hostile with six singleton canonical W blocks was independently reproduced using oppositely directed input triangles with matching cross arcs; it has f=4, s=6, δ=6, and automorphism order 3.

This PASS is bound to the exact final proof SHA-256:

    22cb6aca899d9a485b2432d37a90d13040dfc6353ee41529a243d8489aa974a1

Producer verifier SHA-256, checked unchanged:

    14f993131caf5903b83aa468f77886aeaa714664d2427dc4b4b2052a9aab9fb9

Independent audit verifier SHA-256:

    4aaa93ed94fb8bfac615659ec10ecd1c247427641e49141d98aedb74349fb4bd

The accompanying machine-readable binding is `final-review-receipt.json`. No sharper pole-exponent follow-on is included in this review or required for this packet.
