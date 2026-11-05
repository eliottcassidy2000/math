# Independent review: unit binary boundary cancellation forces a twin

**PASS, no correction required.** This is an all-order physical retraction theorem for rank-one root-clone contexts with unit exponential constant. It also proves the stated exact-from-time-zero fixed-mask consequence. It is not a scalar spectrum classification and does not cover a first transient, changing masks, multi-birth periods or nonunit offsets.

## Frozen source

Directory: `session197-hybrid-constructions/unit-boundary-twin/`.

- PROOF.md: `4f0b7d045e2e8dacd794399b9767631b397ef9c37f3b38a5e82506420e47034c`.
- verify.py: `e7e8e4b03a4a71ac4952dc2be15f652cd4b7bf10b309ca33fe673c6c875fd792`.
- core.py: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`.
- controls.json: `79c39b39f1a63a52b26acbf68497a875924371f130879102f1d5f8837179aa4a`.
- Both source receipts and both independent copied-source replays: `83c2adc0a68efd478b02b65a2ed579b1ce02d50ee83e189315f446ac80c27597`.

The six-file source manifest was checked. Normal and Python `-O` runs match. Original source files were left unchanged. The finite representative bank is corroboration, not the reason the theorem holds at every order.

Named mathematical dependencies are the independently reviewed rank-one boundary proof (SHA `b7935c6f082ff3b6840b71e414a109a15509da3daa8fe91f336796cbac8ebc3b`) and exact minimum-cover amplitude positivity (SHA `2a0451a8473cd17f38a9140c18788fca50325d6f07d8a2716281252f871fbfd6`). Their stated uses are appropriate.

## Structural proof review

1. A nonempty rank-one crossing graph has a one-vertex cover. Taking a converse preserves the scalar clone sequence and boundary and permits the cover center to lie in O. All I-to-(O minus center) incidences are then forward. The inherited signed-gap factorization is exactly the stated negative product, with the indicated endpoint conventions.
2. The possible starts of tournament Hamilton paths are precisely the source SCC, and the possible ends precisely the terminal SCC. The proof of positivity within a strong component is valid: a maximal path ending at the specified vertex cannot accept an outside vertex at its beginning or at any internal gap only when all its vertices dominate every outside vertex, contradicting strongness. No append-after-the-fixed-end step is needed.
3. If the first endpoint factor vanishes, the source SCC of O excludes the crossing center. It dominates every other O vertex, receives from all I, and receives from r. Thus its union with r is a homogeneous module ordered r→C. If the second factor vanishes, the terminal SCC of I has no vertex dominated by the crossing center. It beats that center and every other O vertex, receives from all earlier I components, and beats r. This gives the opposite module C→r. These arguments work when both factors vanish and for singleton SCCs.
4. Deleting C leaves the crossing graph unchanged, not merely of lower rank: no vertex of C was incident with a crossing edge. Therefore the quotient retains a positive rank-one amplitude. All labels and external arcs can be kept literally.
5. The arbitrary nonempty-module scalar law is correct. A rank-one quotient permits one or two nonempty module fragments. Their contributions are a1 H(D) and a2(B_D(2)−2H(D)), respectively. Empty bins occur only in the subtraction term, so no empty-fragment overcount is hidden. This yields a2 B_D(2)+beta H(D).
6. Ordered-join two-bin counts multiply: in each labelled path slot, the first component's path must precede the second component's path. Hence for D=Tr_m→C or C→Tr_m, H(D)=H(C) and B_D(2)=2^m B_C(2). Comparing the two independent exponential bases gives beta(T)=beta(Q)H(C), with integer beta(Q). Unit beta therefore forces H(C)=1, and a strong component with H=1 is a singleton. This produces an actual transitive twin, not merely a count-equivalent replacement.
7. Conversely, any transitive root twin is a two-vertex transitive module. Contracting it preserves the rank-one crossing graph. Its positive clone sequence is shifted by one clone, so A(T)=2A(Q), beta(T)=beta(Q), and the deletion value is A(T)+beta(T). This gives zero boundary without a unit hypothesis. Repeated forced deletion terminates by positive integer amplitude halving.
8. For the fixed-mask consequence, all newborns form one transitive module, so positive times are literal cloning of the first newborn. A nonconstant exact binary affine sequence has positive leading amplitude; positivity of the top run coefficient forces matching rank one. Including time zero in the recurrence is what gives zero boundary. The theorem supplies an old-vertex twin. If the newborn beats that old vertex, the new chain lies immediately before it; otherwise it lies immediately after it. This is a literal adjacency identification with cloning the old vertex, preserving all original vertices and arcs. It neither re-exposes a forgotten port nor authorizes a different mask.

The explicit H=27, A=12, beta=3 source has zero boundary but no root twin, so the unit condition is genuinely necessary. Its stated rows and run coefficients check exactly.

## Independent controls

The independent `check.cpp` imports no source code and does not use the source subset-path inventory. It exhausts all labelled tournaments through order six. Hamilton paths, side endpoint counts and root-run coefficients are counted by literal permutations and gap-factor products. A separate reachability calculation finds SCCs. Ordered two-bin covers are counted by splitting the vertex set into labelled slots and enumerating each slot's Hamilton paths.

The following pass in both optimized and assertion-enabled unoptimized builds:

- 33,867 labelled tournaments.
- 93,590 rank-one roots and 122,012 separate cover-center endpoint formulas.
- 35,928 zero-boundary roots, including 28,728 with unit constant.
- 480 zero-boundary roots without twins; all have nonunit constant.
- 51,696 exposed strong-module controls, checking both the exact amplitude factor A(T)=A(Q)B_C(2) and beta(T)=beta(Q)H(C), not only divisibility.
- 39,888 literal twin deletions with exact amplitude halving and preserved constant.
- 4,836 fixed-mask unit seeds and 827,736 all-entry adjacency comparisons for zero through three births, against explicit old-vertex transitive cloning with the correct placement before or after that vertex.

All positive-index run coefficients above degree two vanish in these rank-one controls. The root coefficient equals an independently enumerated full Hamilton count. Every available cover center is checked, including cases with more than one choice. The nonunit hostile is separately checked from its literal rows.

Independent checker SHA: `3a80dc2ea2838c70018250d5f60dadf563e5ea45a9478431853d81a84ac74e82`.
Both independent receipts: `fa928735c94e471315cdcc0bfbd857bf09abbe317ac6e5a1f3e2eb65308a49b0`.

## Reproduction

The unchanged source snapshot is under `source/`; run its verifier normally and with Python `-O`. Compile the independent `check.cpp` with any C++17 compiler and run with no arguments. It is standalone and regenerates every tested labelled tournament. Optimized and unoptimized outputs must be identical.

Conclusion: zero rank-one run boundary is equivalent to a physical root twin under the unit-constant hypothesis, with the exact smaller-source reduction and fixed-mask consequence as stated.
