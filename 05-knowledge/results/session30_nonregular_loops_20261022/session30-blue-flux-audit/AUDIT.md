# Independent audit: completion boundary flux

**Verdict: PASS for even n >= 2.** The algebraic transport identity is valid.
The complete supplied order-8 bank gives the claimed distribution, independently
recomputed without importing either root checker. A minor domain clarification
is needed in the flux statement: its endpoint formula does not apply to n=0.

## Structural proof audit

Let S be the skew sign matrix, d=S1 in {−1,+1}, and sigma an anti-involution.
Anti-symmetry gives d_(sigma(v))=−d_v. Thus extending sigma to fix the new
vertex preserves the anti-automorphism relation for the completion matrix
[S,−d; d^t,0]. Completion is intrinsic: each old vertex has its new arc forced
by the demand that its completed row sum be zero.

Restriction and extension identify Aut(T) with the completion vertex stabilizer
in Aut(R). The extended sigma has exactly one fixed vertex. A commuting
permutation carries its fixed set to itself, hence fixes that vertex. It follows
that the two centralizers are equal as groups, not merely in cardinality.

A sigma-fixed Hamilton word on T must have second half equal to the reversed
sigma-image of its first half. Its first half contains exactly one vertex per
sigma-pair. Conversely any directed such half-path extends to a fixed Hamilton
word iff its endpoint v satisfies S_(v,sigma(v))=1. All second-half arcs follow
from anti-symmetry. The two constructions are mutually inverse.

For R the central vertex of every sigma-fixed word is the unique fixed vertex.
Deleting it leaves the same half-path and reflected half-path. The two arcs
through the central vertex are equivalent by anti-symmetry, and their acceptance
condition is S_(v,*)=−d_v=1. Again extraction and assembly are inverse.

The automorphism group of a tournament has odd order: an involution would
reverse the arc inside a swapped pair. In the extended ordinary/anti group,
Sylow 2-subgroups consequently have order two. Their involutions are conjugate
under ordinary automorphisms (an anti conjugator can be multiplied by the
source involution). Each symmetric Hamilton word determines precisely one
anti-involution by reversing its positions. Aut acts freely on words. For a
chosen sigma with centralizer c, counting the |Aut|/c conjugates and dividing
by |Aut| proves b=f_sigma/c. This argument applies to both T and R, with the
same c established above. Subtracting the two acceptance vectors proves

    b(T)−b(R) = sum_v h_v [S_(v,sigma(v))+d_v] / (2c).

No equality of total Hamiltonian counts or full metagraph wiring is implied.
The repaired n=0 exception in SELFCONVERSE_POINTING_IDENTITY.md is correct:
that class-count identity has an elementary empty-object check, whereas the
present endpoint formula must simply exclude n=0 or define a separate case.

## Independent finite verification

The Python implementation enumerates all fixed-free perfect matchings, rather
than the root script's low-score/high-score pairing permutations. Testing the
anti condition on all 85 supplied order-8 representatives finds exactly the
same 15 self-converse entries. Completeness of the 85-class input bank remains
inherited from its separate degree-sequence/mass audit, not reproved here.

For each selected involution a pair-mask dynamic program counts directed
half-paths by their last vertex. Its transitions append an adjacent vertex
from an unused sigma-pair. This differs from the root checker, which explicitly
enumerates choices of one vertex per pair and all permutations. All 15 endpoint
vectors match exactly. Every additional anti-involution found by the matching
enumeration gives the same two fixed-word gate counts.

Separately, full_words.cpp iterates every permutation of each T and R. It
counts automorphisms, anti-automorphisms, centralizers, Hamilton words, words
whose positional reversal is an anti-automorphism, and words fixed under the
selected sigma. No half-path gate is used in this second count. Its results
verify the two gates, the centralizer equality, and b=f/c directly for all 15
order-8 classes and both order-6 controls. In total this checks 6,059,520 full
permutations, including all 15 completions of order 9.

Order 8: flux −2 occurs 4 times; flux 0 occurs 7 times; flux +2 occurs 4 times.
The two order-6 controls have blue counts 9 on both sides and completion
Hamilton counts 175 and 171, respectively.

The denominator is genuinely tested, not assumed always trivial. For bits
0100110010011010010100010010 the centralizer has order 3, Aut(T) has order 9,
Aut(R) has order 81, and blue counts are 13 and 11. The other nontrivial
centralizer example has blue counts 11 and 13.

## Literal gate witnesses

For bits 0001011000101000100001000000 and sigma(v)=7−v:

- Half-path (0,4,1,5) yields the directed T word (0,4,1,5,2,6,3,7), but
  inserting vertex 8 at the center does not yield a directed R word.
- Half-path (0,4,2,1) fails the T gate but gives the directed R word
  (0,4,2,1,8,6,5,3,7).

explicit-paths.json records these and further literal paths, checking all four
terminal gate patterns in this negative-flux example and a positive-flux
example with c=3. It directly checks all arcs and vertex distinctness.

## Reproduce and scope

Run `python3 session30-blue-flux-audit/audit.py` from the parent workspace.
Requires Python 3 and g++; it compiles the local full-word checker.
Inputs are the sibling root JSON receipts and the order-8 census bank.
Outputs: receipt.json, explicit-paths.json; run.log records the executed cases.
Source and receipts are hashed in SHA256SUMS. The generated executable is not
a publication artifact. No Git operations or repository edits were made.

This audit establishes the stated identity and finite census, not literature
novelty, an all-order flux distribution, a uniform bound, or black-core wiring.
