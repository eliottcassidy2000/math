# Independent audit: balanced source normal form and residual C1 graphs

Verdict: PASS on the bytes listed below. Isolated normal and optimized replays
pass and reproduce the supplied receipts byte-for-byte.

- PROOF.md SHA256: 00ddcec08dd88e6c91917ff2915bc52662c4e23fe4ba2e1d948bad5d1461a6b3
- verify.py SHA256: bd9ac59f7631f5919b6fc53116b567bfc7887550654b5e9941eb0fb943bb74f2
- core.py SHA256: dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042
- REPAIRS.json SHA256: 432b0334f3925e23adf2884e0ca01f812e643e9b11fa132397da078487d30184

The graph-level converse is sound. Balance forces f=p=0. The physical grammar
implies the distinguished root has outdegree one. A last Z leaves its newborn
as a sink after root deletion; a last J has an E=0 predecessor and hence a
terminal singleton component. Thus root deletion has a unique sink. Deleting
that sink leaves a physical rooted source, and every reconstructed edge agrees
with J(O(W)). The two canonically deleted vertices make the bijection injective
on rooted graph isomorphism classes. The singleton is correctly excluded from
the balanced domain but allowed as the inverse source.

The exact scalar projection D+2U, the three-coordinate indexed orbit, both
finite cutoffs and the all-order growth estimates are correct. Dropping d is
legitimate only because all indexed transitions are independent of d; each
projected word still has a physical source recipe. No inequality-only source
membership is asserted. The sufficient coverage lemma is clearly marked OPEN.

The J Z^k J versus J O^(k+1) rebase is an all-edge identity after changing the
root, with the claimed singleton boundary. Splitting an arbitrary grammar word
at its last non-Z letter yields the exact residual family J Z^k O(W). If W is
nontrivial and k>=1, its final root is the unique degree-one vertex, and deleting
it leaves no sink. Consequently no root is balanced. The displayed residual
scalar formula and its transitive specialization are exact. The graph dichotomy
does not imply that its scalar images are disjoint; the stated95 overlap is
consistent with that distinction.

Independent no-import controls:
- Raw four-coordinate O/Z/J enumeration under D+2U<=50,001 reproduces all
  908,226 source profiles, 788,690 projected triples, 24,958 values and42 gaps.
- Literal Hamilton-word/cover enumeration checks all1,695 rooted sources and
  all270 balanced roots through order5, including canonical reconstruction.
- All386 labelled C1 graphs through order5 are classified independently.
  Every one of the120 unbalanced cases canonically reconstructs as the stated
  residual family with nontrivial inverse source.
- Every newly supplied C1 repair is independently compiled from indexed M_L
  into J and low-bit O/Z moves and counted by a separate subset algorithm.
  All13 succeed, with maximum order13. All27 inherited C2/C3 records are
  unchanged from the previously audited monotone-C1 packet.

These controls are ../balanced-independent.py, ../BALANCED_INDEPENDENT.json,
../balanced-repairs-independent.py and ../BALANCED_REPAIRS_INDEPENDENT.json.
The bounded computation does not prove the cofinite lemma or universal odd
coverage. The conditional implication, using the40 finite repairs, is valid.
