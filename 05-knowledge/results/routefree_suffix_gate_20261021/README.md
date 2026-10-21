# Route-free all-future suffix gates

New Session22 proof and exact finite-scale certificates. Conditional route-free suffix bounds allow sources never to arrive; actual bounded native prefixes, full clocks, source ownership and mod15/mod35 guards remain explicit. No credit is recycled. Uniform 153N coverage over all N remains OPEN.

- [Base proof and native certificates](session22-residual/README.md)
- [Mixed-zonotope strengthening](session22-residual/MIXED-COROLLARY.md)
- [Independent proof and native-prefix review](session22-residual-review/REVIEW.md)

The producer certifies listed all-horizon scales N=4096,8192,16385,32768. Independent native replay covers N=4096 and8192 plus the N=7 hostile, not the two larger scales. Synthetic algebra controls are not a substitute for arbitrary native entry. The general mixed gate explicitly retains pending old sources; the simpler gate requires their arrival first.

## Reproduce

From this packet directory inside the repository:

```sh
python3 session22-residual/check.py --repo-root ../../.. --verify session22-residual/receipt.json
python3 -O session22-residual/check.py --repo-root ../../.. --verify session22-residual/receipt.json
python3 session22-residual/check.py --repo-root ../../.. --cases 4096:149 --verify session22-residual/moment-hostile.json
python3 session22-residual/mixed_gate.py --verify session22-residual/mixed_receipt.json
python3 -O session22-residual/mixed_gate.py --verify session22-residual/mixed_receipt.json
python3 session22-residual-review/audit.py
python3 -O session22-residual-review/audit.py
python3 session22-residual-review/joint_audit.py
python3 -O session22-residual-review/joint_audit.py
```

Review scripts regenerate their receipts. Original selected files and both original manifests are preserved without edits; no exploratory probes are included. The top manifest pins every included file except itself. Historical review references to exploratory files document the audit process, not an additional executable dependency.

## Inherited source pins

Repository dependencies were read at 777c7240f75d3c574c6d0b91e8f6e65a70b2fde6:

- [Frozen frontier and native ownership](../session18_native_annulus/README.md)
- [Native source charge and fixed-route transport](../native_cancellation_105_20261018g.md)
- [Recovered cross-scale source-owned credits](../cross_scale_frontier_credits_20261021/PUBLICATION.md)

Producer and independent checkers authenticate the exact frontier SHA-256. The independent audit source hash is 9d12d89a3c18dfeed1500408628f3840366478400001ec15f78d2368bf3636f8. This is new verified evidence, not a claim to have restored every lost Session21 packet.
