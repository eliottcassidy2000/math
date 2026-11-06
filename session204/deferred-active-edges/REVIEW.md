# Standalone deferred-edge release audit

Status: PASS for `session204-deferred-release` payload frozen in its SHA256SUMS manifest (PROOF.md SHA256 begins 52b5fb4873d5).

All seven payload hashes match. Normal and optimized standalone verification reruns agree, with 65,712 regenerated sparse columns, 1,464 exhaustive literal transitions, 292 literal completions, and 48 larger permuted transitions. No cached matrices, NumPy arithmetic, external seed bank, or sibling-directory imports remain. Counts use arbitrary-precision Python integers and explicit exceptions.

The corrected proof resolves all prior scope issues: equivalence of terminal graph classes to ordinary induced-register constructions is explicit; boundary W and peak W+1 are distinguished; edge processing is restricted to forget/terminal batches; the root remains active for rooted readouts; and L/B require p>=2. The h,s,L,B definitions correctly account for unordered forest components without an extra ordering factor.

The regenerated sparse columns establish formal state closure and nonnegative integer coefficients. Independent literal subset/partition counts validate the complete vectors on the stated finite cases; the all-order correctness remains the unique inverse deletion/edge-removal proof, with monotone bulk-start counts ensuring safe truncation. Finite exploratory search totals are excluded.

No required changes. This certifies a physical partial-graph transfer representation, not numerical occupancy or an eventual threshold.

## Publication packaging addendum

The status-only PROOF.md edit to “PROVED and independently audited in the stated representation scope” is approved. Its new SHA256 is `9b3d4bbf4178deab95b36121eb551f25756eca4df08bd4d11779f3713832b16d`. Replacing only that line by the original reviewed status reproduces SHA256 `52b5fb4873d5fcc4d42d0b3ecb045c6f9fb74577c4e7caec5f0135bd4057c1d9`. The mathematical content is unchanged. The regenerated payload manifest identifies the new proof version; this review's opening hash refers to the original audited version.
