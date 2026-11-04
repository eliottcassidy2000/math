# Independent review: twin-contact quartic response

**ACCEPT.** Reviewed the stable structure proof and implementation independently. The all-order contact formula is correct, the finite sharp minimum is correctly scoped, and the named quartic hostile is valid.

## Mathematical review

Only Hamilton words using the reversed internal arc can change. Such a word has one consecutive twin block, whose contraction is unique. Every adjacency to an old exterior vertex agrees because the twins are homogeneous. With one adjacent port, the incoming and outgoing boundary differences are both a_v−a_u (or c_v−c_u). With two adjacent ports, the two possible contracted blocks x,z,y and y,z,x give precisely the displayed quadratic factors. The four contracted classes are disjoint and exhaustive for noncancelling words. This proves the formula at every order, including endpoint positions and the singleton quotient.

Each port has at most two path incidences, so its attachment variables have degree at most two; total degree at most four follows. There is no global-converse argument in this contraction proof and no cancellation of its quartic terms is asserted.

The hostile adjacency is a valid tournament; the named twins are0,4 and the only changed edge is0→4. Its quotient has the directed spanning cycle0→3→2→1→0, proving strongness of both pair expansions. Fixed root0 covers all rooted quotient types after relabelling, which preserves polynomial degree; hence the finite census genuinely proves the minimum old-vertex order five.

The scope statements are correct: T4's global-converse cubic is unchanged, no full-response collision is produced, and the low-degree C2 sign can still distinguish twins. A quartic coefficient is an obstruction to reconstructing the complete ambient polynomial by falsely setting degree-four terms to zero, not an obstruction to every low-degree separation test.

## Independent computations

`independent.py` imports no author or repository code. It enumerates each Hamilton word on the two-port extension, determines its forced zero/one attachment bits, and adds one to every compatible Boolean assignment. This differs from both author's endpoint DP and author's symbolic product expansion. Möbius inversion is then applied to the literal truth tables.

- All2048 named attachment counts for the five-vertex hostile are rebuilt from5040 permutations per tournament.
- All16 frozen witness-query counts agree with the author's receipt.
- Every coefficient is retained independently in `full_coefficients.json`.
- The named coefficient at mask586, a_1 a_3 c_1 c_4, is exactly2.
- Nonzero coefficients by degree1,2,3,4 are4,20,56,32; all higher degrees vanish.
- The entire75-quotient census through quotient order four is independently repeated with literal truth tables. Its degree histogram is1×degree1, 2×degree3, 8×degree3, 64×degree4, grouped by quotient order1,2,3,4.
- Independent normal and optimized receipts byte-match.

The all-order proof is not inferred from these finite checks.

Reviewed stable source hashes:

    PROOF.md 346cb02dcb19d388d65d4282c6c69cc87fcf84128280da1fdd9de0678794013d
    verify.py dd3e6ff0f306befa75d031365ba968722de78c05ed1c8fb5add4195cb1536ac2
    normal.json = optimized.json
    a7e40bce8251c2d663e3acae7edd464815bb31098536cabd73d41f74cc85dabc

## Final added identity: independently accepted

The final source adds

    [a_i a_j c_i c_v] Δ = C_zj(Q−i) − C_jz(Q−i).

This proof is correct. The absence of twin a-bits in the requested monomial eliminates A_x and every relevant quadratic twin factor. The remaining P terms have at most two exterior incidences, whereas three distinct exterior attachment variables are required. Thus only c_v A_y contributes. Since y is already adjacent to z, c_i forces its other neighbor to be i; a_i and a_j force x between i,j. The only possible combined blocks are z,y,i,x,j and j,x,i,y,z. Their signed cubic attachment coefficients are respectively+1 and−1. Deleting the two ports and i gives exactly the two ordered-cover sets with untested z,j gap. Singleton paths are included, with no lost endpoint case or factorial.

The necessary degree≤3 condition is correctly stated as necessary only. No sufficiency claim is made.

The updated independent checker validates all400 ordered exterior-pair instances in the complete census. It compares literal-truth-table Möbius coefficients against a separately written subset partition convolution of endpoint-fixed path counts on Q−i. This differs from the author's permutation-cut enumeration of ordered covers. Normal and optimized checks agree. `census.json` records400 successful checks. The hashes above refer to this final updated proof and script, superseding the initial review hashes.
