# Independent narrow audit of incoming W1–W2

**Verdict: PASS for the W1 fixed-port argument and W2 common-support theorem, with the native K4 witness independently reproduced.** This is an audit of concurrent incoming work, not a claim that Session104 newly discovered these theorems. No publication was performed. The earlier `audit_descent` artifacts remain frozen.

## Provenance and exact scope

Read using the GitHub connector, pinned to commit `ba6481180a0c988e8ed8a11e08dd42155eb0dbf9`:

[`05-knowledge/results/wedge_descent_168_20261029/THEORY.md`](https://github.com/eliottcassidy2000/math/blob/ba6481180a0c988e8ed8a11e08dd42155eb0dbf9/05-knowledge/results/wedge_descent_168_20261029/THEORY.md), file blob `eeec19ecbd00756b5a52d446a3b6d7f2061d2a07`.

Only W1–W2 were audited. This report does **not** certify W3–W7, the whole incoming packet, its other receipts, or the claimed n=6,k=3 sharpness computation. The all-k potential lemma inherited from Session103 was also read locally in the preceding audit. No raw or lazy git-blob fetching was used.

## W1: correct finite rank and nilpotence proof

Canonical cycle cuts provide the required bijection Phi between ordered start tuples and ordered end tuples. This is essential: B is the column-matched table M, not an arbitrary repeated use of M with an unspecified seam. Telescoping gives weight difference w(b)-w(a). It is positive for every noncanonical cover because every noncycle arc has positive integer cost. Thus B=I+N with N nonnegative and every nonzero N entry strictly increasing tuple weight.

For distinct k-tuples, the weight range is k(n-k), including k=n where all paths are singletons and N=0. A chain of h+1 nonzero N entries is impossible, proving N^(h+1)=0 without an asymptotic argument. The finite binomial power, inverse, and resolvent formulas then follow algebraically. Nilpotence index equals longest support-path length plus one because nonnegative path contributions cannot cancel. The statement that actual degree may be lower is necessary and borne out by the K4 test (index 3 rather than the general upper bound 5).

The q-weight similarity is valid for invertible q, or formally over Laurent/rational functions. At q=0 its displayed D(q)^-1 expression is not literal matrix inversion; the cover-polynomial itself still has its natural specialization. This is a standard domain convention, not an obstruction to the proof. The weighting only records boundary potential difference, and the note correctly does not advertise it as an internal physical clock.

## W2: correct menu criterion, with its sign and carrier guards

A union DAG gives a common ordering and hence simultaneous unit-upper-triangular form. Conversely, traverse any simple directed cycle in the union support and select a factor containing each edge, in traversal order. The identity choices already contribute one to every product diagonal; the selected cycle contributes a strictly positive extra term to one diagonal. Nonnegativity prevents cancellation. Trace is then greater than the matrix dimension, so the product cannot be unipotent. A witness of length at most the carrier dimension follows by taking a simple cycle.

The expansion over ordered subsets of factor positions is valid without commutativity. A support-path bound truncates all terms above length h; a fixed finite menu gives a uniform entry bound and therefore polynomial growth in depth. These conclusions need the same fixed labelled carrier, nonnegative off-diagonal entries, and zero off-diagonal-part diagonal. They do not follow merely because each separate operator is individually unipotent or determinant one. The incoming note preserves these guards.

For the claimed hostile product, a nonnegative matrix with a diagonal entry greater than one has (A^d)_(i,i)>=A_(i,i)^d. Thus the stated exponential diagonal growth follows directly, not merely from an unsupported spectral heuristic.

## Independent literal native controls

`check_common_rank.py` uses standard-library integer arithmetic, no producer imports, and no `assert` guards. It builds the strong K4 seed with cycle 0→1→2→3→0 and chords 0→2, 1→3. Every permutation of its four vertices and each of three cuts is independently tested as an ordered pair of nonempty disjoint paths. This enumerates each ordered spanning two-path cover exactly once, including singleton paths. It then constructs M and the canonical matched B directly.

Both normal and optimized runs pass **251 explicit checks**, with byte-identical receipts:

- 22 literal ordered covers, 12 ordered port states.
- N nilpotence index 3; the general N^5=0 bound holds.
- Integer finite inverse and binomial powers through depth 10 agree with direct multiplication.
- All 24 vertex-induced re-markings independently preserve nilpotence.
- Trace histogram of B times its re-marked conjugate: **17 cases at 12, 6 at 16, 1 at 22**.
- The transposition swapping vertices 1 and 3 gives **trace 22** and a diagonal entry greater than one.
- **Seven trace-12 cases still have cyclic support union.** This strengthens the practical warning: trace equal to the dimension is genuinely insufficient even in this native census.
- For every cyclic union in the 24-case census, the checker extracts a simple support cycle, chooses factors edge by edge, and verifies the constructed finite word has trace greater than 12. It verifies the constructive W2 witness rather than relying solely on the two-factor trace.

Reproduce:

```sh
python3 -B audit_incoming/check_common_rank.py > audit_incoming/receipt.json
python3 -B -O audit_incoming/check_common_rank.py > audit_incoming/receipt.optimized.json
cmp audit_incoming/receipt.json audit_incoming/receipt.optimized.json
```

These finite controls do not replace the general proofs, and do not test the n=6,k=3 degree-9 sharpness fixture. Credit W1–W2 to the incoming concurrent packet and the underlying Session103 positive-grading mechanism; this audit provides independent confirmation and the trace-12 cyclic-union refinement.
