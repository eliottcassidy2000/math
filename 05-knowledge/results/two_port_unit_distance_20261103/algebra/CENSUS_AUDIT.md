# Independent review of the exact two-port census through order nine

**Verdict: FINITE-EXACT conclusion accepted, with the inherited-bank
dependency stated explicitly.** I reviewed
[`../census/RESULT.md`](../census/RESULT.md), `check.py`, `exact.cpp`,
`run_cpp.py`, `audit.py`, and their stored receipts. I independently checked
the current hashes and receipt agreements and replayed the smaller
permutation audit. I did not rerun the complete order-nine census or
re-enumerate the inherited tournament bank.

The conclusion is precisely that every labelled tournament of order at
most nine is determined by H,S,E,C, hence by its at-most-two-port response.
It does not assert an unbounded reconstruction theorem or an equivalence
between H,S,E,C and the whole two-port response.

## 1. Correct observable and exact arithmetic

Session27's response proof, sections7--8, recovers the ordered C matrix,
both labelled endpoint vectors S,E, and H from actual queries with at most
two external vertices. The census retains the directed entries C_uv, not
only their symmetric sums or current signs. H need not be a separate key
entry because H=sum S=sum E at each fixed order.

The subset recurrence for `ends` counts directed Hamilton paths on each
subset with a specified last vertex. Applying it to the converse gives
the first-vertex counts. The C formula sums

    end_A(u) start_(V\A)(v)

over every nonempty proper subset A. Each ordered spanning two-path cover
has exactly one such first-path vertex set, so the formula is bijective.
There are no diagonal contributions because the two paths have disjoint
vertex sets. Empty-path cases are correctly excluded from C and do not
alter the endpoint vectors.

All values fit well within uint64. Each Hamilton path is a vertex
permutation, and each ordered two-path cover is a vertex permutation with
one of n-1 cuts. Thus the **total** number of ordered two-path covers,
and hence every C entry, is at most `(n-1)n! <= 8*9! = 2,903,040` in this
universe. Individual summands in the nonnegative subset recurrence cannot
exceed the same bound. Endpoint counts are at most9!. Canonical labels
and counters are also safely bounded: each input has at most9! leaves,
and even191,536*9! is below2^64. Adjacency codes have at most36 bits;
the unsigned vertex masks require only9 bits. No numerical comparison
uses floating point or a probabilistic fingerprint.

## 2. Why the refined canonicalization is complete

Initial colours are exact pairs(S_u,E_u), assigned by their sorted values.
At every refinement round, a vertex retains its old colour and receives
the sorted multiset of triples

    (colour(v), C_uv, C_vu),       v != u.

This construction is equivariant under every sidecar isomorphism. Because
the old colour is the first part of each key, colours never merge. If no
class splits, their order remains unchanged and the procedure terminates;
otherwise the number of classes increases, so there are at most n-1
strict splits. In particular, no heuristic signature is being treated as
an isomorphism certificate by itself.

After refinement, both implementations enumerate **every permutation
inside every remaining colour cell**. Cell order is fixed by its invariant
colour. Their Cartesian product is an invariant family of labellings.
Every sidecar automorphism preserves every cell and transports that whole
family to itself. Taking the minimum full sidecar key over this family is
therefore an exact canonical form. It need not equal the lexicographic
minimum over all n! unrestricted labellings; that stronger property is
unnecessary and is not assumed in the proof.

Two isomorphic sidecars have identical candidate-key sets and hence the
same minimum. Conversely, any equal full key directly gives a sidecar
isomorphism through the two displayed labellings. Thus comparing these
keys tests sidecar isomorphism exactly.

## 3. Why all minimizing leaves establish a labelled result

Fix one minimizing labelling. Every sidecar automorphism composed with it
is another permitted minimizing labelling. Conversely, any two minimizing
labellings differ by a sidecar automorphism because their complete labelled
sidecar keys agree. Hence the set of minimizers is a torsor for the entire
sidecar automorphism group; it is not merely a sample of that group.

The Python implementation keeps the set of tournament codes at all these
leaves. The C++ implementation keeps their minimum and maximum. These are
equivalent tests: there is exactly one resulting tournament code if and
only if minimum equals maximum. This is equivalent to every sidecar
automorphism preserving T. The reverse containment, Aut(T) contained in
Aut(sidecar), holds because the sidecar is defined by intrinsic path
counts. The separate equality with the inherited automorphism-group order
is therefore a useful control, not a premise needed for this implication.

Now suppose two labelled tournaments have equal sidecars. The no-class-
collision check first puts them in the same tournament isomorphism class.
Choose an isomorphism to a common representative. Their two namings then
differ by an automorphism of its sidecar, and the all-minimizers check
shows that this permutation also preserves the tournament itself.
Consequently the original labelled tournaments are equal. This validates
the passage from the class census to labelled injectivity.

## 4. Receipts, replays, and the inherited universe

The stored C++ receipt reports191,536 distinct sidecar classes at order9,
zero collisions of either kind, zero disagreement with the inherited
automorphism orders,1,362,537 canonical leaves, and a worst input of9!
leaves. The earlier rows report1,1,2,4,12,56,456,6,880 classes at orders
1 through8, all separated with the same symmetry condition.

I checked the following against the current files:

- All three normal/optimized receipt pairs agree byte for byte.
- The complete C++ order1--8 rows equal the Python rows.
- The C++ and wrapper source hashes match `cpp9.json`.
- The Python and audit source hashes match their receipts.
- Both bank hashes match the bytes actually present in the repository.
- `audit.py` replays successfully, with387 active checks. It independently
  counts paths and two-path covers by literal permutations for all76
  classes through order6, and canonicalizes them by all n! labellings
  without refinement. It also independently counts six explicitly named
  larger classes, two each at orders7,8,9.

The final `cpp9.json` SHA256 checked in this review is

    abc2dfb5e42ec43dfccd3c624bee34e0d9b65d4da51d0e0f8279c814dbdee95d

and the order-nine bank SHA256 is

    b4e559ad6cb31980852034f9c7f61e366a7cb479fee4d9a641ad3356745200d3.

The C++ and Python producers share their algorithmic design; their
agreement is not represented as a wholly independent mathematical
algorithm. The literal-permutation audit supplies that separate control
within its stated smaller universe. The wrapper's normal and optimized
Python modes both compile C++ with-O3; they do not test two different C++
optimization levels.

Completeness of the representative bank remains inherited from the
October26 enumeration packet's proved construction and exact audited
outputs. Merely observing the expected count would not prove completeness
on its own. The new census assumes that established universe and computes
a new faithful observable on it; neither this review nor the result note
claims to re-prove its enumeration here.

## 5. Consequence and stopping boundary

Since equality of full two-port responses implies equality of H,S,E,C,
the finite injection result excludes every two-port counterexample through
order9. A counterexample to universal two-port faithfulness must therefore
have order at least10. This is a valid lower bound on a possible
counterexample's size, not evidence that one exists.

Conversely, a future collision of H,S,E,C does not by itself refute full
two-port faithfulness: the full response retains additional ordered
three-path-cover information. The result note explicitly preserves that
distinction. No overstatement or blocking correctness issue was found.
