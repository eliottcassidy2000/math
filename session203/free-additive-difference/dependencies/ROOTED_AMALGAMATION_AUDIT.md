# Independent audit: positive ordered-root amalgamation and binary inverse

**PASS, no correction.** The positive retained-state representation, arbitrary-core binary formula, complete finite inverse and explicit repair have the stated scope. This report does not promote any later core-bank census or all-order coverage claim.

## 1. Snapshot

Source directory: `session197-odd-cycle-path-polynomials/`.

- `ROOTED_ORDERED_AMALGAMATION.md`: `0166889f66a6363506d336dcf7ae3d34a768f4ae79be5d387cfb426d78fe2883`.
- `verify_rooted_ordered_amalgamation.py`: `777daa85b00826349db0a29a9a36bb1539a4cb2103d782f4a921baf6b05b6a12`.
- Both frozen receipts: `3092333f52384fdf8afbe8570961fec501ce0e1fe136a31753f39594eb30f1b9`.
- Imported `verify_retained_root_sewing.py`: `62d824eb95905c6c408b38079eb7f56a3c7e8226539f9bc7588f3d1c3855beed`.
- Imported `verify.py`: `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.
- `RETAINED_ROOT_SEWING.md`: `c5a7a18611363985b49d4ab455ad42160150ab012157c7d4c01aa6f838332655`.
- `RETAINED_SEWING_SOURCE_SCOPE.md`: `e301d1085e9d66104ecfac6f514f4f2f2c1efeb9df40c5426aa237f82bc874c6`.

The general root-amalgamation h/R/L identities were also independently derived in the parallel C3 specialization, published at commit `6ed2d7a0f1e1bfb26b9cc865ddc9f35dd0ad4567`, `session197/lowcycle/root-amalgamation/`. That overlap is acknowledged rather than treated as a priority claim. The positive representation and arbitrary-core inverse are the distinct extension checked here. No frozen source was edited.

## 2. Proof review

The shared-root construction is associative as an actual graph operation because all nonroot blocks have the same total block order. Converse reverses that order. The four Hamilton-word cases are disjoint, and the subtraction of H from either reserved-cover count removes exactly the root-singleton reserved path. The two ordinary bins remain labelled throughout. The R formula correctly allows the other bin to be empty. The endpoint product requires the second factor to be non-singleton; its explicit identity exception is essential and is stated correctly.

Substitution of x=R-d and y=d-e gives the four displayed rows with nonnegative integer coefficients. Endpoint deletion, root insertion and the root-singleton cover supply the needed injections. The strict slope criteria are exact: a single in-neighbor or out-neighbor supplies a two-vertex reserved path, while its absence forbids every non-singleton reserved path of the corresponding orientation. The matrix anti-homomorphism follows from these identities and associativity; direct substitution, including the multiplicativity of the deletion's two-bin count, also verifies it on the whole four-dimensional space. The last two rows are proportional to the input d coordinate, so the nonidentity rank bound is correct. The negative minor is a valid separate guardrail against total nonnegativity.

The prefix reads q most-significant-first with O=1 and Z=0. The suffix reads t least-significant-first with Z=1 and O=0. These opposite conventions are necessary for the stated converse symmetry and are handled consistently. Leading prefix Z letters are neutral only for the singleton's retained four-coordinate state; terminal suffix O letters are neutral only for the final scalar. The proof removes no internal padding. Substituting the prefix state, core matrix and suffix readout gives the symmetric formula (8) with all increments after h_C nonnegative.

For the finite inverse, the interior bound follows from t>=1 and the positive beta*q*t term. The inner denominator alpha+beta*q is strictly positive. A zero-slope axis is exactly a universal-source or universal-sink root and reduces to H*F(index). Removing trailing binary ones preserves F; a remaining positive even index 2v satisfies F(2v)>=v+1, giving the complete finite normalized interval. This is a same-count graph replacement on a terminal axis, not permission to delete arbitrary internal vertices. The positive-slope axes and origin are also retained. The binary-length order bound follows from the product bound in the interior and the separate linear bounds on the axes. Neither uniqueness of q,t nor coverage of every odd integer is inferred.

The fixed core signature, nonnegative production matrix, negative minor -96 and scalar repair calculation all agree. The twenty-vertex graph uses the general-root core operation, so it does not contradict the separately restricted one-edge-core census. Its independently counted H is 462661.

## 3. Independent controls

The audit imports no target implementation. `literal.py` is the exact standalone directed-word counter and retained-label graph constructor from the separately reviewed C3 root-amalgamation packet. The target instead uses subset Hamilton-path DP. `check.py` builds its own positive matrices, binary histories and finite bounded images from those literal profiles.

Normal and optimized independent runs agree on:

- 8,236 actual retained-state pairs, covering every rooted labelled first factor through order three and non-singleton second factor through order four.
- 784 exact four-by-four matrix composition identities on nontrivial small factors.
- 1,331 literal associative triples.
- 1,152 actual binary-history graphs and every arc of their converse isomorphisms.
- All eighteen small-core scalar images through 200, matching the frozen receipt. This independent enumeration exhausts pairs within the proved total bit-length bound, rather than using the source's rectangular scan or its nested slope bounds.
- 4,096 zero-axis normalization controls.

The repair graph was independently reconstructed and matched every stored adjacency row. The separate C++ checker counts Hamilton paths and anchored cycles by inclusion-exclusion over walks in induced subsets, not subset path DP. Optimized and assertion/undefined-behavior-sanitized runs both give

    order 20, H=462661, C=3.

The additional C=3 readout is a control in this report, not a modification of the frozen theorem. It also follows from the cycle-cut product: the nonempty canonical prefix has root-start count one, and the core's root-ending count is three.

Copied-source normal/-O replays reproduce the frozen receipt exactly. Finite controls test the formulas and implementation; the arbitrary-core inverse and order bounds are proved in the packet.

## 4. Reproduction and hashes

From this audit directory:

```
python3 source/verify_rooted_ordered_amalgamation.py > source-normal.json
python3 -O source/verify_rooted_ordered_amalgamation.py > source-optimized.json
python3 check.py > independent-normal.json
python3 -O check.py > independent-optimized.json
g++ -std=c++20 -O2 walks.cpp -o /tmp/ordered-amalgamation-walks
/tmp/ordered-amalgamation-walks repair-graph.txt > walks-normal.txt
g++ -std=c++20 -O1 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined walks.cpp -o /tmp/ordered-amalgamation-walks-debug
/tmp/ordered-amalgamation-walks-debug repair-graph.txt > walks-debug.txt
```

`check.py` writes only the independent repair-graph input. Keep `literal.py` alongside it. The full `source/` subtree contains the exact two imported Python helpers and named mathematical context needed by the source replay. All implementations use standard libraries.

- `check.py`: `c952bc2c275a3f8658491698aa47f8debaed55cddb34c480c02c2cfcbd392e9d`.
- `literal.py`: `ccd82862bae7dd1b73be1412f74f6e3f4e1d339aaa3eba711ee1c998401c6bfc`.
- `walks.cpp`: `7b1fc1015fbd8e09ab49465b33a48e4fe9a5800e43abf55d15afbffa08e92b59`.
- `repair-graph.txt`: `5be42aac568a4073fa2013e0bcc03e7cee21fbcaa080e03a33ec27b942325036`.
- Both independent Python receipts: `d41685fbf3a040a846c67e0d8990de14e0a97ea5ca4c7d55dc9b98d68647318d`.
- Both walk receipts: `17c4982133a5ba497b7d48a08df75ded6341aac464be99d17d63c82118ce600f`.
- Both source replay receipts: `3092333f52384fdf8afbe8570961fec501ce0e1fe136a31753f39594eb30f1b9`.

No unresolved issue remains within the frozen packet's stated scope.
