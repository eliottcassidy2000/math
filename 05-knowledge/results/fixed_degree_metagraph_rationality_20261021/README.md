# Fixed-degree tournament metagraph rationality

New Session22 synthesis: rational counting series for every fixed fundamental black degree, general denominator bounds, explicit leading constants for degrees 3^r, and vanishing proportion of bounded-degree nodes under uniform merged-node sampling. This uses inherited finite strong-atom and center-divisor laws, not a complete wiring or connectivity classification. Formal executable controls are not tournament censuses; unknown finite atom polynomials remain unknown.

- [Proof](session22-fixed-degree/PROOF.md)
- [Independent review](session22-fixed-degree-review/REVIEW.md)

## Reproduce

Run from this packet directory, with Python 3 and SymPy installed:

```sh
python3 session22-fixed-degree/check.py
python3 -O session22-fixed-degree/check.py
python3 session22-fixed-degree-review/check.py
python3 -O session22-fixed-degree-review/check.py
```

Producer receipts record 463 checks and 9,236 formal words; independent receipts record 1,004 checks and 24,419 formal words. The reviewer script regenerates its receipt.json. Sibling names preserve its original relative paths. Files are copied byte-for-byte from the reviewed final packet, and the manifest pins all included files except itself.

## Inherited sources and novelty boundary

These repository sources were read at commit 622e4451384c6d6ec6ee56bba0f71dc50732512d:

- [Effective finite strong-atom bound](../effective_black_defect_order_bound_20261020/PROOF.txt)
- [General center-divisor law](../deficiency_three_20261020/centers.md)
- [Three-presentation classification](../metagraph_classification_20261021a/README.md)
- [THM-2454 pure-blue alphabet and center law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md)
- [Session21 reconstructed mathematics](../session21_recovery_20261021/RECONSTRUCTED_MATHEMATICS.md)
- [Session21 reconstructed metagraph proofs](../session21_recovery_20261021/RECONSTRUCTED_METAGRAPH_PROOFS.md)

The r=1/trivalent special case and center decomposition are inherited, not new discoveries. The proof does not depend on a stronger Session22 center cutoff. The reviewed proof SHA-256 is e428906db3b67dbf5c85eef86460a012be8a9574a06bcd62759e73778f23446d.
