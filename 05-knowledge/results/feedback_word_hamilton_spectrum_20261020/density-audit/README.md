# Binary-product extension: publication and reproduction

20 October 2026. This independently reviewed extension preserves the earlier three-piece packet unchanged. It proves the all-real-cutoff presentation asymptotic W(X) = X log 2 + O(sqrt X), and only a one-sided upper bound for distinct one-apex/intercyclic values. It does not prove an exact distinct-value density, global forbidden values, or the complete odd tournament spectrum. No independent priority claim is made.

The [full proof](../BINARY_PRODUCT_DENSITY.md) and [independent audit](AUDIT.md) retain their frozen reviewed bytes. The imported verify.py and cited proofs.md/intercyclic_rays.md belong to the earlier packet and are pinned as dependencies in the publication inventory. This extension adds files only.

From the parent result directory, run:

```sh
python3 verify_binary_product.py
python3 -O verify_binary_product.py
python3 density-audit/independent_density.py
python3 -O density-audit/independent_density.py
sed 's@session16-hamilton-spectrum/@@' BINARY_PRODUCT_SHA256SUMS | sha256sum -c -
```

All four copied-layout checks reproduce the archived output bytes. The producer's numerical analytic controls are distinguished from the proof; the clean-room audit additionally uses rational logarithm enclosures. Historical source_snapshot commands and the prior-seal verification transcript record the original audit layout, not commands to run from this directory. External source snapshots and duplicate dependency copies are not included.
