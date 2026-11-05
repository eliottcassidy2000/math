# Independent audit: binary-module rebasing

**PASS, no correction required.** The count-preserving deletion/birth operations are literal, retain the four specified old ports, and connect the stated fixed-amplitude representations. The infinite strict one-final-gadget separation is valid. This is a rebuilding operation that uses deletion, not an induced extension or a universal target-coverage theorem.

## Frozen snapshot

Source: `session197-hybrid-constructions/binary-module-rebasing/`.

- PROOF.md: `75f242f5297f74525ced120a83cd50ff47321f950f6596c722cd16b50ae2c587`.
- core.py: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`.
- rebase.py: `12d1e67433a1d56a4069523aaf3b27bf566effa0456315e77f2f36778c396453`.
- verify.py: `d9f5ef7e78d0dfe7e6ccbd219c22dc5ce8ffcce18598e8319fc8ec27bc8b4a05`.
- VALIDATION.json and both source-mode receipts: `ca4087a7b4db480a83fe115f837673041fc69a0524ce03e04debcc35a39f6a4e`.

Normal and optimized source runs were executed in a copied snapshot and reproduce these receipt bytes. Originals were not changed. The ordinary final-gadget law is inherited from the already reviewed ROOTED/pair-entry packets, not newly asserted by this construction.

## All-order proof review

### Graph and count certificate

The seven-block quotient, in order a,U,b,X,Z,Y,R, has rows[14,4,32,118,103,67,7], with sizes[1,t,1,2,j,2,m]. This independently reproduces the orientation specification. The displayed cycle visits every block and uses legal crossing edges; positive module sizes make it a genuine Hamilton cycle. The marked last R vertex has exactly the t+2 outside vertices as out-neighbors, so different t values at fixed sum and m cannot be rooted-isomorphic.

For U, its crossing edges are all incident to b. For Z, they are all incident to a. For R, they are covered by a,b. Earlier and later vertices of the same transitive module create no new crossing edges: their order forbids the necessary return orientation. These covers therefore remain valid uniformly when the other modules have arbitrary positive sizes.

Deleting clone runs from a Hamilton word produces vertex-disjoint return fragments. A fragment starting in an out-neighbor and ending in an in-neighbor must cross an O→I edge. Choosing one such edge per fragment gives a matching, so a size-d cover permits at most d+1 runs. Assigning the ordered transitive clones to the named nonempty runs gives the stated surjection polynomial. Since the clone size is positive, the zero-base term vanishes, leaving exactly the claimed exponential span.

Separate fixed-variable spans with these uniform bases imply membership in their tensor product: apply the one-variable interpolation formula successively. The12-point grid is complete for the2×2×3 span, and its exact counts determine the coefficient30 on2^t2^j3^m and constant1. No zero-sized module is used. This is an all-order finite coefficient certificate after a structural run bound, not empirical extrapolation.

### Literal reversible replacement

Deleting the first Z vertex leaves the same block orientation with j reduced by one. A new right-born vertex beating only b has exactly the incidences of the last U vertex, so the explicitly specified relabeling gives T(t+1,j-1,m). Conversely deleting the last U vertex and adding a left-born vertex losing only to x0,x1 gives the first vertex of the enlarged Z block. All four old ports survive both deletions, and their identities map to the corresponding four ports afterward. The two operations are inverse up to those relabelings.

The positive-parameter requirements j>=2 or t>=2 keep the intermediate graph inside the strong family. The formula gives H_intermediate=(H+1)/2 and the following birth restores2H_intermediate−1. The proof does not silently use a count-preserving induced extension. It also correctly declines to claim that the same four ports suffice for arbitrary future R-cloning when U is large.

### Root response and strict next-target separation

For a path required to end at the outside root, the same U/Z run bounds still apply. The ending count is therefore fixed by its four-point grid. After adding w that beats only the root, paths ending at w contribute h−e, and paths continuing after w contribute exactly R: the first old path may be empty and otherwise cannot contain the root. This proves H(T+w)=h−e+R with the stated conventions.

The augmented U crossing cover is{b,w}; the Z cover remains{a}. Hence its count lies in the3×2 tensor span. The six exact evaluations determine the displayed augmented count, and subtraction gives the exact R formula. The marked root remains the same old vertex throughout this argument.

For the two equivalent-count representations, the displayed singleton-child values are correct. R>=h and e>=0 give2R_B>=F_B(1), while F_B(L)>=LR_B for L>=2. On the A side, F_A(1)<N_B<F_A(2) and N_B<3R_A exclude every L. All inequalities are strict at j=2 and become stronger with2^j. Thus neither child target is obtainable by one ordinary final gadget from the other representation. This excludes only one generation, not multiple gadgets or other graphs of the same H.

## Independent implementation

The standalone `check.cpp` imports no source code. It constructs the seven-block quotient above and expands its transitive blocks, rather than copying the source's pairwise orientation cases.

- 294 graphs with t,j=1,...,7 and m=1,...,6: all tournament incidences, the explicit Hamilton cycle, and all named crossing covers.
- 504 physical delete/birth maps: every retained adjacency and every retained-port identity compared under separately defined block labels.
- 56 actual graphs through order14: Hamilton words counted by literal directed-path DFS. This counter records path counts for every induced subset, not by endpoint DP.
- 15 response/augmented controls: the same literal words directly give ending counts and reserved two-cover counts, and count the added-w graph independently.
- The strict separation inequalities are checked through j=40, supplementing the symbolic inequalities above.

`check_tensor.py` performs exact rational elimination on the independently counted graph grid. It recovers all12 H coefficients, all4 ending coefficients, and all6 augmented coefficients, including the zero coefficients. Normal and optimized Python results agree. An optimized C++ build and a separate unoptimized build with library assertions and undefined-behavior sanitization produce identical receipts.

Hashes:
- check.cpp: `da0bb70a9777285f9338be2ed5f4f038419d18e8ee22d85f6bdc3b30a5151096`.
- check_tensor.py: `e5988e375f2f4b62812d3777d9b45756530e01877c629b766c54bb091b8bad04`.
- independent-normal.json and independent-debug.json: `177ab849afd844ca9983c6d1bc2038f6ce2a41e6e86ddf3bb27b518c4bca9024`.
- tensor-normal.json and tensor-optimized.json: `362a9bdc27198712af83f0b5438fe90ab11bfabf09ab918fba581a8d41d4cf47`.

Reproduce from this directory:

    g++ -O3 -std=c++17 check.cpp -o /tmp/s197-rebase-audit
    /tmp/s197-rebase-audit > independent-normal.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -fsanitize=undefined -std=c++17 check.cpp -o /tmp/s197-rebase-audit-debug
    /tmp/s197-rebase-audit-debug > independent-debug.json
    python3 check_tensor.py
    python3 -O check_tensor.py

The finite controls corroborate the exact coefficient and graph-map certificates. The scope remains positive modules, the two stated four-port rebuilding operations, the m=1 marked-response formulas, and the strict one-final-gadget consequence. No arbitrary-target coverage or universal ternary interface is inferred.
