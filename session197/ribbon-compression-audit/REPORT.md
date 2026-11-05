# Independent compression and sparsity audit

2026-11-05. **PASS**, for corrected `PROOF.md` SHA256
`7cdde933d32db6568e412402673fc7ebf780de7cec27a84b183f05d1be0859a2`.

## Conclusions

The logarithmic witness theorem is valid as stated for fixed or variable
integer `w>=1`, `h>1`: every realizable value h has a span-w realization on
at most `(4w+1) floor(log_3 h)` vertices. The algorithmic NP and numerical-h
polynomial-time consequences are correctly restricted to fixed width.
The width-1 bound `A(X)<=exp(pi sqrt(2 log_3 X))` is also proved correctly.
Neither theorem supplies an unproved global spectrum classification.

## Proof checks

1. The head-comparison merge preserves each input path as a subsequence and
   gives a directed path. At a switch, the prior head beat the still-current
   opposite head. Its deterministic output recovers both inputs by restriction,
   proving injectivity for arbitrary, possibly interleaved vertex parts.
   Iteration proves the product lower bound for any finite partition.
2. Tournament SCCs form an ordered join. A Hamiltonian path cannot leave
   and reenter a component or traverse components against the join order.
   Counts therefore multiply. Removing singleton SCCs leaves all nontrivial
   SCCs intact and preserves the product. Reindexing can only shorten each
   backward edge, so the span does not increase. No SCC interval assumption
   in the supplied vertex order is needed.
3. In a nontrivial strong component every vertex v has both in- and
   outneighbors. If no outneighbor points to an inneighbor, the outneighbor
   set cannot reach v. This proves that v lies on a directed triangle.
4. The two triangle orientation patterns are exhaustive. One has a backward
   extreme edge, giving diameter at most w; the other has both adjacent
   backward edges, giving diameter at most 2w. Centers spaced `4w+1` have
   disjoint integer neighborhoods of radius `2w`, including boundary-truncated
   neighborhoods. Thus `ceil(N/(4w+1))` vertex-disjoint triangles exist.
5. The product inequality gives `h>=3^ceil(N/(4w+1))`. Since the ceiling is
   an integer, it is at most `floor(log_3 h)`, proving the claimed constant.
6. For fixed w the frontier dimension and transition coefficients are
   constants. A witness uses `O_w(log h)` symbols, each of w bits. At cap
   h+1 every coordinate uses `O(log h)` bits. Each transition uses a constant
   number of fixed-coefficient additions/multiplications and comparisons,
   so verification uses `O_w((log h)^2)` elementary bit operations. This also
   covers computation of the integer logarithm by repeated multiplication
   by 3. Exhausting at most `2^{w(4w+1)floor(log_3 h)}` words, up to a
   constant factor for shorter lengths, is polynomial in numerical h,
   not in binary input length. No variable-width NP assertion is made.
7. For width-1 sparsity, nondecreasing tribonacci values imply
   `h_n>=3h_(n-3)`; the three base cases give the displayed exponential
   lower bound. Mapping a product representation to its exponent vector
   only overcounts spectral values. The partition-product estimate for
   total weighted degree at most `3 log_3 X` is valid, including the empty
   product. Expanding the positive logarithmic series and optimizing t
   gives exactly the stated constant `pi sqrt(2 log_3 X)`.

## Edge cases and correction

The initial draft said unconditionally that deleting singleton SCCs
preserves H, although all vertices disappear for transitive tournaments and
the companion implementation assigns H(empty)=0. This was reported and
fixed before signoff: the deletion theorem now assumes `H(T)>1`, and h=1
has a separate one-vertex witness. No empty-path convention is needed.

Width 0 is outside the stated main theorem and is harmless: every supplied
order is transitive, the only positive spectral value is 1, and one vertex
is a witness. The n=1 and h=1 cases are explicit. Nonstrong tournaments may
contain arbitrarily many singleton SCCs at a fixed count: tested examples
have H=3 at orders 3, 30, 100, and 300, all reducing to a three-cycle.
This shows why the theorem must bound the reduced witness, not the
unreduced input. Noninterval SCCs cause no problem.

## Independent computation

`audit.py` uses independent graph construction, repeated-search SCCs,
predecessor-removal Hamiltonian-path counting, and an index-cursor merge.
It checks:

- Every ordered tournament through order 6: 33,867 graphs, including all
  orientations, their actual backward span, and every SCC configuration.
- Every nonempty unordered two-part vertex partition through order 5:
  15,834 partitions and 31,322 input path pairs. These include noncontiguous
  and interleaved parts, not just prefix/suffix splits.
- Sixteen random bounded-span examples at widths 1–4 and orders up to 200,
  including deletion/reencoding and exact count preservation, with up to
  20 simultaneously selected disjoint triangles.
- Explicit width-0/h=1 cases, singleton padding, the generator growth
  inequality through n=500, and exact width-1 monoid counts through 100,000.

The source verifier also passes from the audit-local frozen copy. These
finite checks support implementation details; the preceding proof checks
establish arbitrary-order validity.

Reproduce with `python audit.py`. Results are in `AUDIT.json` and
`audit.log`. Original source copies are under `source/` with hashes in
`SOURCE_SHA256.txt`. No author source files were edited by the auditor.
