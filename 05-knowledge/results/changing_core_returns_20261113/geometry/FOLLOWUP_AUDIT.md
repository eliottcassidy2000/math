# Final independent reviews: minimal 881 patches and fixed balanced coordinates

**Audit conclusion: mathematical statements pass within their declared
operation and census scopes.** One wording clarification was requested in the
881 note: preserving first birth H1689 is distinct from the later balanced
entry (b,t)=(5207,1689). No correction to its count or construction was needed.

## Minimal endpoint patches at 881

[NEUTRAL_881.md](../preparation/NEUTRAL_881.md) separates three operations:
single-arc neutral paths, simultaneous endpoint-neutral arc patches, and
subsequent neutral tail deletions. The minimum-width conclusions refer to a
fixed named source/root and a fixed original vertex set before deletion.
They do not give a minimum in an arbitrary mixed birth/deletion/edit grammar.

The filter root-indegree>=2 is necessary for the chosen return's native guard.
A source root has (h,s,q,F)=(h,h,h,h), so Delta=0. Any indegree-one rooted
tournament is literally the newborn of its deletion parent at its unique
predecessor; the N formula gives Delta=0 again. Then b_k=3t_k excludes the
guard. Thus discarding all lower-indegree width-four patches cannot discard
a successful endpoint, whether strong or weak.

The C++ enumeration was reviewed directly. Its increasing arc-index recursion
visits each of the binom(55,4)=341055 patches exactly once. The H counter is a
valid subset-last recurrence, and 32-bit cells are safe at order11 since even
11! is below2^32. F uses 64-bit arithmetic. The script computes q and F before
overwriting its table for the root-start count. All16 retained states are
recounted independently in the producer, and the filtered width-three universe
is cross-checked against a complete separate Python scan.

Therefore the finite exclusions at widths0..3 and the successful weak width4
endpoint establish the stated minimum4. Filtering at width4 leaves only one
strong endpoint, with Delta60, so that endpoint is blocked. Together with the
explicit strong width5 role exchange, this establishes the strong fixed-order
minimum5. These conclusions rely on the declared exhaustive enumeration;
this audit reviews its loop and independently recounts endpoints without
rerunning the large patch scan.

`preparation_audit.py` independently recounts all44 supplied H881 endpoint
records: five width-two records, nineteen width-three records, sixteen
filtered width-four records, and the four original/full selector-cube
endpoints. Every recorded first birth misses1691. This last statement is
about those records, not an unfiltered exclusion of every possible width-four
patch whose root has indegree less than2.

The four-arc chosen B endpoint has SCCs on old names
(0,1,2,3,4,5,6,7,10), then(8), then(9). Removing9 then8 preserves the whole
state at each step and yields exactly

```
rows=(354,228,425,179,231,192,140,1,250),
root index8 (old name10), state=(881,169,201,977).
```

Its independently counted RNR output is(5195,1689,1689,5195). All16 filtered
width-four states, their strongness, and their stated wheel predicates through
k400 are independently checked. Normal and optimized runs agree.

## Complete support at fixed balanced deletion/start count

[FIXED_BALANCED.md](../returns/FIXED_BALANCED.md) has the correct quantifiers:
strong marked sources with full profile(H,t,t,H), for a **fixed** odd t.
The finite bank is not a bank for every zero-readiness-defect source at once.

The inverse balanced root has indegree1. Its deletion parent U has
H(U)=t, zero start count at the predecessor, and F_U=H-t. Strongness after
birth is equivalent to that predecessor lying in the final SCC of U. The
necessity follows because every old vertex must reach that sole entry to
the newborn. For sufficiency, the final component reaches the predecessor,
all earlier components reach it, and the newborn beats all other old vertices.
There is an old vertex other than the predecessor when t>=3.

The nontrivial strong SCC list is finite: each factor H>=3 divides t and the
published strong-order bound bounds every component. If the marked final SCC
is strong, the amplitude is2^a times the product of earlier two-bin counts
times its marked end-cover count. With only one strong component, a>=1 is
necessary for zero start; with two or more, a>=0 is sufficient. If the marked
final SCC is a singleton, all other singleton factors contribute2 and the
amplitude is2^(L-1) times the product of strong-component two-bin counts.
Every placement of the singleton gaps is retained by the templates.

Both converse constructions are valid and yield the exact balanced profile,
so the two displayed unions are equalities, not just necessary residue lists.
Every vertex of a strong tournament starts a Hamilton path: take a longest
path starting there. Any vertex that can neither be appended nor inserted
must beat every vertex of the path, contradicting strongness if an outside
vertex exists. This justifies the zero-start boundary used in the proof.

The t1 case is also correct: its parent is transitive, its marked predecessor
must be last, and the parent has at least two vertices. Thus H=1+2^a for a>=1.
The prime members are exactly the prime Fermat numbers by the stated odd-factor
factorization. No claim about their number is implied.

Odd-part normalization with the minimum exponent gives disjoint full tails
of dyadic rays. For fixed t their count is O_t(log X), and bounded ranges of t
remain finite unions. Excluding coverage of all primes or7p additionally uses
a standard prime-abundance fact; a minimal route is that each finite dyadic
union has a convergent reciprocal sum, whereas the prime reciprocal sum
diverges. This is unrelated to proving prime infinitude within a selected ray.

## Independent t71 support and positive construction checks

The extended `birth_budget_audit.py` recomputes the odd-part/minimum-exponent
compression from all independently recounted81 literal core-root records and
their P(2) values. It reproduces all62 rays, excludes odd part405, and gives
odd part101 minimum exponent3. It then independently counts all four actual
neighboring-ray graphs, checking every N/R adjacency update and strongness:

```
(879,71,71,879), (1687,71,71,1687),
(3303,71,71,3303), (6535,71,71,6535).
```

The weak parent's deletion count at its selected old root is not assumed71;
only its actual H71, start0, and F808 are needed. The newborn then has deletion
and start coordinates71 by the N law. These are strong balanced scalar
families with zero readiness defect, not unit-ready returns.

Normal and optimized audit results agree. The core-census completeness and
the general fixed-t finite-bank theorem remain distinct from these finite
literal checks.
