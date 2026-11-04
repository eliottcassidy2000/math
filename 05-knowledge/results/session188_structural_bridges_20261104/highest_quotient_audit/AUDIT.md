# Independent highest-changed-quotient reduction audit

**PASS on the frozen files pinned below.** Reviewed the complete proof and checker, independently reconstructed the finite census, and checked every linear-node physical separator by literal integer Hamilton-path DP. The reduction does not assert full injectivity or that every residual prime converse collides.

## All-order proof

- The original strong-module tree is explicitly distinguished from strong connectivity and from the possibly different canonical tree of the target.
- Every changed edge has a changed quotient at the least original strong node containing its endpoints. A highest changed ancestor therefore contains the edge. Highest changed nodes are an antichain, hence pairwise disjoint, and all incidences outside their interiors are unchanged.
- The scalar-equality dependency gives hereditary equality of complete ordered-cover profiles. It is not silently replaced by equality of H alone, or by an assumption that a child's operation is one converse.
- At a highest changed linear quotient, the target is another transitive child order. The arbitrary-block cut theorem applies with the same physical names, equal child q-profiles and unchanged ambient incidences.
- Other highest changed nodes remain modules in the actual port extensions: their incidences to the chosen node are uniform, and the ports copy that chosen node's incidences. Replacing each other node by its equal-profile target is exactly Hamilton-count preserving. This restores the actual full target without treating distinct changes as additive response contrasts.
- The strict lower bound comes from positive run-substitution coefficients and the one-run coefficient H(R). No general exact scalar multiplication of response differences is assumed.
- Prime quotients of size three are retained. Linear descendants below a changed prime root are not highest and are not incorrectly excluded by this theorem.

A draft control `TT2 -> C3` was replaced before freeze because the linear parent flattens the proposed TT2 strong-tree node. The final control is `C3[TT2,C3,TT1]`, with genuine separate highest linear and prime children. This correction affects only the illustrative control, not the proof, and the final rows and original strong tree were independently checked.

## Independent computation

`independent.py` constructs original-module edge-incidence spans rather than iterating the producer's module-family operation loop. It builds all modules and original strong nodes directly, computes every highest changed quotient, checks hereditary child profiles by independent subset-cover convolution, and checks prime quotient structure directly. Its physical edge-bit convention is opposite the producer's enumeration convention; it compares labelled adjacency rows and aggregate census, not mismatched code identifiers.

Across all labelled tournaments through five old names:

- 132,041 source/target span presentations, including the empty family;
- 21,629 scalar-equal presentations, 20,530 nonidentity;
- 18,400 highest changed linear nodes, each with a verified physical cut query and quantitative lower bound;
- 2,250 highest changed prime nodes and prime-only presentations;
- 120 presentations with multiple highest changed nodes.

Every one of the 36,800 source/target query counts was checked both against the complete retained named response bank and by independent literal integer subset DP. The complete by-order census matches the producer. The two six-name controls are checked independently: highest masks {3,28} in the corrected mixed case; only root mask63 highest in the whole-converse C3[TT2,TT2,TT2] case, despite changed linear descendants3,12,48.

Normal and optimized independent Python receipts are byte-identical. The frozen producer normal and optimized receipts are also byte-identical. All gates are explicit exceptions. Finite validation supports the implementation; the all-order conclusion rests on the structural proof and the already reviewed arbitrary-block theorem.

## Frozen source pins

- `highest_quotient_reduction/PROOF.md`: 95437f44021ca708bfa2fde42637e53179ddb81a5cf93341c3bffa602e1b4678
- `highest_quotient_reduction/verify.py`: cf7aead4d1c32c41ab784fffdd66a76bdc12c755cb8a4125f5371f48ab8d338f
- Both producer receipts: 3dd401befb63c9f90121b99167154d90eaff85fba22e5e140b0ed8724d07e889
- `covered_laminar/ORDERED_BLOCKS.md`: 07e2c5f0ed81b57a8d476f61adc5d74da030d563ad6f3621fd7c8fcabc7f3a16
