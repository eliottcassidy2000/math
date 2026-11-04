# Independent audit of the two apex suffix obstruction

**PASS after one harmless wording guard.** The bridge agent reviewed W1–W3,
the interior-insertion repair, and their boundary conditions independently
of the producer. This is a mathematical and executable audit, not Lean
formalization or a proof of the full spectrum conjecture.

## Independent implementation and exact checks

`independent_audit.py` imports no producer or repository computation code.
It implements the gap decomposition with36 states, keeping the empty first
gap separate from a nonempty first gap ending correctly. Its reachable
physical span has dimension28; the producer's merged24-state implementation
has physical dimension20. Both retain the same17-dimensional count observer.
These different intermediate dimensions are consistent, not an identification
of the physical state spaces.

The independent implementation verifies:

* The literal two-apex graph count for every word of length0 through4:
  341 cases, using a separate remaining-vertex Hamilton recurrence.
* All289 entries in the17-by17 strong-source profile matrix by literal
  graph counts; determinant150994944; the constant function representation;
  all68 query row-closure identities on the independently closed36-state
  reachable span.
* All two-apex strong source classes at orders3–6, with totals1,1,6,34.
  Complete independent Hamilton-frame generation at orders5–6 gives6 and35
  strong classes, so exactly one strong order6 class has no such presentation.
* All four strong H15 classes, reconstructed from1088 Hamilton frames.
  Every admissible ordered apex pair is recovered from its actual relation,
  checking transitive deletion directly. This gives17 distinct marked words.
  Direct counts on all549 possible suffix graphs through order8 show that
  every source marking misses47. Both order5 classes also miss49, while
  exactly one order6 class supplies49.
* The H9 source represented by112 misses31 under every admissible old
  marking and all suffixes through the complete strong cutoff. The actual
  extension10102 has H31 and preserves every old arc under the explicit
  vertex embedding `(0,1,2,4,6)`.
* Eighteen induced weak-source witnesses for the unboundedness argument
  (`m=1,...,6`, each non3 letter) and the positive controls29,31,47,49.

The independent normal and optimized receipts agree byte for byte. Their
hash fields identify the exact producer script and receipt audited.

## Proof review

**W1 passes.** Closure of the actual reachable vector span, rather than a
finite word-depth sample, certifies that all suffix counts lie in the query
span. The displayed strong profile minor has full rank, and the constant
function is in that span. An affine identity vanishing on all strong words
therefore vanishes on all actual words. The independent36-state construction
reproduces each necessary algebraic interface.

**W2 passes.** Applying that identity to the weak words3^m is legitimate only
after W1. A non3 suffix vertex and a suitably chosen apex induce a cyclic
pattern around the old transitive block, forcing unbounded H by induced
monotonicity. Use `m>=1` in this induced count: at m=0 there are only two
vertices, and the displayed `2^m+1` cycle formula does not apply. This wording
boundary does not affect the unbounded argument. An all3 suffix is exactly
sink padding, and the distinct strong counts3 and5 force the unique affine
parameters a=1,b=0.

**The finite suffix cutoff and W3 pass.** Remove trailing3s, which are
universal sinks. If anything remains, its final non3 vertex beats an apex
and receives an edge from every old chain vertex. Hence it joins the old
strong source, and each intervening appended vertex is reachable from the
old source and reaches that final vertex. The trimmed graph is strong, so
the cited Moon–Busch cutoff applies. This argument is essential: an arbitrary
weak extension would not inherit the same order bound.

For H15, the strong floor25 at order7 bounds every possible strong neutral
reset source by order6. The four-class enumeration is complete, every
permitted old marking is included, and the finite suffix menus are
independently recounted. Thus the corollary covers arbitrarily many strong
H15-neutral resets followed by one suffix phase. It does not cover native
moves interleaved with growth, weak neutral sources, or moving marks during
the suffix phase.

**The interior repair passes.** Old chain positions0,2,4 of10102 have word112
and retain the old relation. The new vertices are interior insertions, so
the suffix-only no-go does not apply. The stated family recurrence has the
correct initial values31,53 and recurrence roots1,2. Its comparison with the
proposed ternary+4 expression is equal only at m=3; the isolated success is
not a scalable induction.

## Audit implementation correction

An initial version of the independent marking extractor computed sorting
keys from the same list while Python's in-place sort temporarily exposed it
as empty. This incorrectly discarded known legal markings. The hostile
requirement that exactly one independently recovered H15 source reaches49
caught the bug. The final implementation freezes the deletion-vertex tuple
before sorting, and directly enumerates all its marked relation words. No
producer theorem or accepted receipt depended on the faulty audit version.
