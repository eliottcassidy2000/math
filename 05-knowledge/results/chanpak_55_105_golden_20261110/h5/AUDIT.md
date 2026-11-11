# Independent audit of the fixed-core compiler

**Verdict: mathematical proof and critical exact consequences checked.**
Audit by the55/105 source-class lane, independently of the compiler author.
This is a local packet audit, not a claim of an external refereeing process
or of a proof of Chan–Pak.

The proof was checked against the inherited squarefree reciprocity identities.
In particular the sign identities for the endpoint sums exchange starts
and ends, and the `epsilon*delta` contribution cancels between paths
crossing both flanks without visiting K and paths visiting K. The OCF
classification includes paths wholly in either flank, crossing both without
K, and meeting K; the residual weight is H of the deleted core path's
complement. No contiguity assumption is used. The quotient correction
`(q−2h) A(B+1)` is consistent with this enumeration.

The finite inverse bounds use only nonnegative terms. The uniform cases
require the separate A=0 or B=0 normalization: stripping terminal1 bits
preserves g, and a nonzero even representative2r has g(2r)≥r+1. Thus
the inverse proves an all-order numerical absence, rather than stopping
at an arbitrary graph order. The H5 classification proof correctly uses
the two-triangle budget and the R→L arc forced by strongness. The claim
about all H1 sources was clarified to account separately for the empty
source, whose sole output1 is already present.

Additional checks performed during this audit:

* 160 general compiler cases on four cores beyond the author's main three:
  TT2, the order5 code571, the order6 code10023, and C3→C3. Every core
  used four named cuts including both uniform cuts and ten pairs(A,B),
  including asymmetric pairs and values of binary width3. Graphs were
  assembled independently with the helpers in `../targets/check.py` and
  counted by its Hamilton-word DFS, not the compiler's subset DP.
* A rectangular inverse check through173 examined all26 named cuts of
  singleton/C3/Q4 and every0≤A,B≤346:3,130,634 arithmetic triples. This
  rectangle contains every normalized solution at that target bound,
  and does not use the compiler's row-dependent enumeration limits. Its
  three supports exactly match the terminating inverse. Their union
  misses only7,21,173; the H5 support excludes105.
* The173 witness and its entire deletion deck were independently recounted
  by the other lane's Hamilton-word DFS:173 and
  `(27,71,51,27,51,71,77,9)`.
* The eight original source hashes in `SOURCE_PINS.json` matched at the
  pinned source HEAD. After incoming changes, compare them with that
  immutable Git version rather than mutable working-tree documents.
* Following recovery of Session197 `ONE_CORE_BINARY_INVERSE.md`, the
  compiler author audited its precise restriction O=V(K) minus b and
  the physical prefix/suffix bit conventions. The coefficient substitution
  `(h,q,e,s,u,v,t)=(H,C,L,C-H-Q,S,H-E,H-S+L)` reproduces its equation8,
  including both parity corrections. The checker now verifies all seven
  named C3/Q4 ports and448 pairs0<=A,B<=7. This post-review expansion
  credits the physical binary compiler and bounded inverse as inherited;
  the arbitrary-cut extension and complete H5-parent image remain the
  present claims. It does not change any previously audited support.
* The refreshed normal and optimized receipts after those restriction
  checks are byte-identical, with SHA256
  `52bbce1a6bf91396e6f67434ac671b9b9bd59ae92c3e119084314dafcbd6a568`,
  and their embedded script hash matches the current checker. This is
  the post-restriction receipt, replacing the original review's hash.

The additional arithmetic audit reused the inspected formula and its
coefficient routine; it is independent of the target-bound enumeration,
not a second derivation of those coefficients. Those coefficients were
audited mathematically and through the independent literal graph checks.
The author's full target1001 run remains the receipt for the further
misses335,663. The audit independently establishes the critical first
new gap173 and its sharp deletion floor9. Neither this floor nor the
fixed-core compiler proves closure under repeated ternary returns.
