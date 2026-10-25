# Session56 independent audit: valuation carrier, connection atlas, Cartesian fibres

Status: PASS AFTER TWO CORRECTIONS, verified 2026-10-25. Scope is the three frozen producer packets under `inputs/`, not a global certification of their inherited research programs. No producer files or remotes were changed by this auditor. Repository sources were read at immutable commit `4f1760fd1eff346767c032524c932d5e317401d2`; `PROVENANCE.json` records all producer hashes and 18 source-pin checks against the actual git blobs. Source access used git-show for sparsely materialized files. No claim about later remote freshness is made.

## Findings and repair verification

1. Atlas card 1 originally claimed the Cartesian map preserves the numerical minimum. That is false without the array/key sidecar: `[1,1]` and `[5,5]` have identical tree and argmin signatures but different values. Producer repaired the statement, explicitly retained the sidecar, and added the hostile to its checker. The corrected frozen snapshot passes.
2. Cartesian calibration originally said all cited values for n<=6 are smaller than Catalan. They are equal at n=0,1. Parent repaired this to “at most”, retaining the first strict obstruction at n=7. Corrected frozen snapshot verified.

No other blocking mathematical error found in the declared scopes. Improvements above are required for this verdict; audit snapshots contain both repairs.

## Valuation carrier

The exact theorem follows from equality of the first k least-significant digits iff congruence modulo p^k. Negative integers use Euclidean digits. Equal records and diagonal queries have infinite depth; a finite precision cap yields only a lower bound at saturation. Distinct values can be separated by p^H greater than their spread, then compressed by adjacent LCP depths in LOW-DIGIT lexicographic order. Numeric order fails on 0,1,2 at p=2. Empty/singleton and repeated-label separation conventions are correct.

CRT collision relations are intersections of named-pair relations. The actual modulus-35 arrays yield crossed nonlaminar clusters, so one same-leaf rooted hierarchy cannot contain both partition families. This does not prohibit a DAG/product or separate trees. Fixed-template Lo Shu thresholds are 3 at p=3 and 7 at p=2; neither trie counts root-vector lifts nor proves integer squarehood.

Collatz fuel was independently derived from the native shortcut orbit of 27, not assumed from the p-adic analogy: 59 steps, 37 odd steps, endpoint 23, b=1100931843921811423, D=126176846412426125. The rational fixed point has odd denominator and is nonintegral. The repeated-word cylinder gives K(n)=floor(v2(Dn-b)/59); a legal block consumes exactly one. Neither exhaustion nor arbitrary finite-cylinder depth gives universal coverage or a global Collatz rank. Source snapshots match the committed Session53 guarded-pumping and repetition packets byte for byte.

## Atlas source/status check

All ten cards retain correct source scopes after repair. Cards 1–2: argmin observer versus value sidecar, and idempotent AND only on fixed selected masks. Card 3: whole partition ancestry obstruction, not an impossibility of expression trees. Card 4: named boundary signatures, no persistence implementation claim. Card 5: THM-3503 innovation regrading v_m retained, no assumed isometry at raw depth. Card 6: THM-3382 linear versus exponential addresses, no identification of density with tree depth. Card 7: exact charge/carry identity, not color-only arithmetic. Card 8: STANDARD Collatz map distinct from shortcut map; THM-4592 all-odd-period bijection remains only FINITE-EXACT through period 22, consistent with its correction entry. Card 9: disjoint min-plus composition, with THM-2897 gate sufficient only. Card 10: abstract metric-monoid higher-arity hostile, not an exhibited knot symbiont or native tournament counterexample. All proposed applications remain OPEN; inherited PROVED labels mean scoped source claims rather than newly audited entire dependencies.

## Cartesian fibres and capacity

Independent stack-built Cartesian parents verify both stable-RMQ implications, including ties. Hook fibre formula follows by allocating n-1 ranks between left/right subtrees. Exact enumeration through n=8 verifies f(T)=n!/product s_v, Catalan support, and weighted mass sum 1. The random-BST statement is distributional only.

Tournament sidecar bound is an encoding-capacity bound for exactly n nodes with inorder labels and a fixed-width b-bit sidecar. It is not a lower bound on computing the scalar counting sequence or on canonicalization time. Integer rounding b>=max(0,binom(n,2)-floor(log2(C_n n!))) is correct. OEIS A000568 was directly read: offset 0, t7=456>C7=429; n0 and n1 are equal. This remains CITED, not an independently produced tournament census. ETH primary Bender/Farach-Colton LATIN2000 PDF confirms the RMQ-index/Cartesian-LCA context; the local equivalence proof handles stable ties explicitly. Fischer’s primary `Optimal Succinctness for Range Minimum Queries` PDF at https://ae.iti.kit.edu/download/fischer10optimal.pdf was also read; its abstract and §1 define RMQ as a minimum position, consistent with the repaired value-sidecar distinction.

## Independent controls and reproduction

Run from this directory:

    python3 -B check.py > normal.json
    python3 -O -B check.py > optimized.json
    cmp normal.json optimized.json
    python3 -B atlas_controls.py > atlas.normal.json
    python3 -O -B atlas_controls.py > atlas.optimized.json
    cmp atlas.normal.json atlas.optimized.json

All load-bearing checks raise explicit exceptions and survive optimization. No producer code is imported. The first checker covers 3,280 ternary arrays, all permutation fibres through n=8, 3,124 signed/repeated valuation lists, native Collatz replays for 48 cylinder lifts and 288 finite caps, 4,096 tropical triples, and integer capacity rounding through n=100. The second covers 10,000 exact-charge pairs, 8,000 carry cocycles, exhaustive small-mask intervals, named coloring ports, address growth identities, and 1,200 standard-Collatz steps. Finite controls supplement the mathematical proofs, not replace their universal quantifiers.

`provenance.py` creates closed producer snapshots, checks all producer normal/-O replays and saved certificates, and verifies historical source pins. Its git source path assumes the audit checkout exists; the frozen sources and input hashes permit offline inspection after publication. All three producer replays passed in both modes and were identical. SHA256SUMS seals this audit.

## Operational note

Bounded startup and relevant guardrails/protocol/core-source routes were read. The targeted processor check could not run because this read-only checkout has no `.machine-id`; no identity was fabricated. This is a coordination-tool limitation only, not a mathematical or source-verification blocker. Publication and remote synchronization remain with the parent.
