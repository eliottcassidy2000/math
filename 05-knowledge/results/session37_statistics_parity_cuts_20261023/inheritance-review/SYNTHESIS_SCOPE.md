# Session 37 synthesis scope review

Status: PASS for the assigned mathematical scope, with small clarity recommendations below. No pending publication input or frozen audit artifact was edited.

Reviewed:
- `session37-root/SYNTHESIS.md`, SHA-256 `f38b21ac9cd28112e4198146dbc450232e2cd5417b886644eaee071bee9aee1a`.
- Current source tree at `12aa833cda`.
- Incoming `01-canon/theorems/THM-1970-H-sits-at-the-formula-sharp-P-edge.md`, SHA-256 `b8568559ee66ca13601220360603af3913b93d73a7a4284ecd07081f9f533ae8`.

## Section 2: adjacency versus skew

PASS. The matrices A, S=A-A^T, K=-S^2, vertex-label X and permutation action P_g are explicitly distinguished. The five-vertex H=3/15 example is skew/energy-cospectral, not an equality of adjacency characteristic/permanental data. Thus it does not contradict corrected recovery from the adjacency pair through order eight. The same distinction is stated explicitly in the arithmetic packet, Theorem D and README. No inference from equality of scalar spectra to equality of full Cayley matrices is made.

## Section 7: corrected inheritance

1. **Threshold and coverage: PASS.** Recovery through eight, failure first at nine, and source-extension propagation agree with the two complete finite computations and explicit witness. The count 2,097,152 is fixed-path presentations, not isomorphism classes; 3217 counts invariant fibers. The finite-exact/computer-assisted status at eight is stated, and no unexplained closed reconstruction formula is claimed.
2. **Rank versus functional recovery: PASS.** The true order-eight rank obstruction rules out affine reconstruction only. Three actual tournaments provide the concrete distinction. The synthesis does not infer an equal-coordinate collision from a carrier-space kernel.
3. **Realized order-nine kernel: PASS.** For the first witness minus the second, the audited packing counts give
   `C_T-C_U = z^9(u^3-u)+z^8(v^2-u^2)`,
   where z marks occupied vertices, u odd cycles and v even cycles. The D36/D45 difference cancels in this parity-resolved polynomial. It vanishes at permanent `(u,v)=(1,1)` and determinant `det(I+zA)` `(u,v)=(1,-1)`; equivalently it vanishes at the `det(I-zA)` cycle-sign specialization `(-1,-1)`. At `(z,u,v)=(1,2,0)` it equals two, matching H=2843 versus 2841. This is an actual packing-polynomial difference, not an unrealized formal trade.
4. **Bosonic THM-1810 repair: PASS.** The permanent/transitive versus edge-product/all-tournaments distinction and the minimal six-versus-eight witness are correct. Positive combinatorial Wick signs do not rule out coefficient cancellation. The surviving detection-depth statement is the already-repaired THM-1790 projective argument, not the discarded positivity slogan.

## Incoming THM-1970 attribution

PASS. Commit `515de698d65e892280f929e9dba1051b1afe313b` inserted the correction now present at remote `12aa833cda`. Its banner is dated 2026-10-22, while the commit timestamp is 2026-10-23. It explicitly withdraws the implication from counting hardness to fixed-resolution non-determination, explaining that a determining invariant may have a computationally hard decoder. It also retracts unsupported general fixed-k growth/reconstruction conclusions and distinguishes skew traces from full decks.

Calling this a concurrent incoming repair is accurate and does not claim it as this session's discovery. It supports the synthesis's distinction among information sufficiency, linear decoding and computational difficulty. This review does not promote the inherited finite tables or historical paragraphs below that correction banner to new audited claims.

## Suggested clarity-only changes before final synthesis publication

- In section 7 item 2, define the kernel variables inline: “where z marks occupied vertices, u counts odd cycles and v counts even cycles.” Call it the **parity-resolved cycle-packing difference** rather than an unspecified missing-information polynomial. Its present formula is correct, but otherwise the variables are undefined in the synthesis itself.
- In section 7 item 3, name “the already-corrected THM-1790 projective detection-depth argument” instead of merely “the valid earlier detection-depth argument.” This makes the surviving proof lineage unambiguous.
- Link the incoming correction by its full canonical filename `THM-1970-H-sits-at-the-formula-sharp-P-edge.md`, and optionally give its 2026-10-22 banner date. No need to include Git mechanics in the user-facing synthesis.
- Add ordinary spaces to compressed phrases such as “order8”, “all2,097,152”, “with3217”, “H2843” and “by6 versus8”. These are presentation fixes, not mathematical conditions.

No blocking mathematical error was found in the assigned sections. The first inheritance-correction checkpoint need not wait for these synthesis-only clarity edits or for further research.
