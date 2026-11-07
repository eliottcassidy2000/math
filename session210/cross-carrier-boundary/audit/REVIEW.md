# Independent review: cover polynomial versus projective topology

**Verdict: PASS.** Independently checked the two explicit order-six tournaments, the complete proof, the standalone verifier, and both normal and optimized source replays. No correction is requested. The previously published projective-carrier packet was not altered.

## Exact checked source

- PROOF.md SHA256: `8cbf174f3b4b9be3036cc5ea7d460c0b4d00ff858207793c8a0acad650bcaa71`
- verify.py SHA256: `cf48ab0d9737783ff3b30372991702ce76d597a171d0b7ca60142aafa3121fb2`
- Both source replays exactly match RESULT.json, SHA256: `83fcb6da674f1957911e7910933944604333c62c411d3ca2534e02cae2d1482c`.

## Independent two-case computation

The independent checker uses the earlier projective worker's Hamilton subset-DP and cycle-mask machinery, not the source verifier's direct Hamilton-permutation counter. It independently enumerates the forward-adjacency histogram, expands its binomial transform over rational coefficients, and sums all bin assignments at q=0,...,4. It also directly counts Hamilton paths after zero, one, and two gate insertions.

For BOTH specified row arrays:

- Forward-adjacency histogram: `(17,105,238,238,105,17)`.
- Exact polynomial: `P(q)=q^6+8q^4+8q^2`.
- Direct q=0,...,4 bin totals: `(0,17,224,1449,6272)`.
- Direct nested Hamilton counts at depths 0,1,2: `(17,241,4001)`.

A has packing polynomial `1+8z`, eight isolated projective cells, Betti vector `(8,0)`, and packing number one. B has polynomial `1+6z+z^2`, exactly one compatible pair of cycles, Betti vector `(5,1)`, and packing number two. Its quotient is therefore a circle plus four isolated points. All cycle words are supplied in the independent RESULT.json.

## Proof and scope

The first-order-six claim is correct: below six vertices, no two odd cycles of length at least three can be disjoint. Thus the projective quotient is a discrete set with size determined by H=P(1). The explicit pair attains failure at six.

Equality at every nested gate depth follows from the proved polynomial gate operator, not from the three direct-depth checks. Equality of face-sink incidence polynomials follows from the previously audited forward-adjacency transform. Neither equality says the oriented permutohedra themselves are isomorphic. The proof correctly distinguishes all these scopes and makes no spectrum-completeness or literature-priority claim.

## Reproduction

`python3 check.py` replays only the two controls. `carrier_support.py` is an unchanged copy of the previously published projective checker (SHA256 `d3bc82ae26da32e0bbff6e7b1ee0cd4e7e011c3d9d9aaf3a052bafd00a962d95`); importing it does not run its census. No full census was rerun. `source.replay.json` and `source.optimized.json` record the separate normal/optimized source checks.
