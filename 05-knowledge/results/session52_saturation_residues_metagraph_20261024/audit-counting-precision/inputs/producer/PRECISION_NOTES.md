# Narrow precision repair for THM-587 and THM-588

2026-10-24. This is a separate release supplement. The independently audited `session52-counting-lenses` packet remains unchanged. No shared repository source has been edited.

## Exact edits

- THM-587: replace the average of `dim(W^sigma)` by the average of `trace(sigma on W)`, with `W` explicitly the direct sum of the Fourier levels. A permutation's trace on the labelled-tournament function space counts fixed labelled tournaments; the dimension of its fixed subspace counts its cycles on that basis and is a different quantity.
- THM-587: consistently call `P_n(-1)` the graded antipodal trace and the fixed-class count. The historical proposed Kaczynski / Ky-Fan / computational-homology direction remains credited and referenced, but any actual complex and differential must be supplied separately before making Euler-characteristic, Betti-number or Borsuk-Ulam claims. The labelled complement's freeness is explicitly scoped to n>=2.
- THM-588: the exact adjacency eigenfunction is `c_0(T)=c(T)-binom(n,3)/4`, not uncentered cyclicity. Its mean uses uniform labelled tournaments, equivalently labelled-orbit-size weights on isomorphism classes. The proof adds the constant Fourier coefficient and shows `A c_0=(d-4)c_0`, `L c_0=4c_0`. Historical numerical correlation observations remain labelled historical.

The original theorem IDs, source attribution, count tables, spectral values, dependencies and result paths are preserved. These edits do not audit every other historical interpretation or performance statement in the files.

## Verification and publication boundary

`source-pins.json` specifies the exact before/after SHA256 of both targeted files. Publisher must check both before hashes before applying `precision.patch`, and both after hashes afterward. If a target differs, stop and reconcile; do not fuzz-apply or silently substitute a source version.

Run the non-mutating verifier with the actual repository path:

    python3 -B verify_precision.py --repo /path/to/repository
    python3 -B -O verify_precision.py --repo /path/to/repository --receipt /tmp/precision-optimized.json

The verifier copies only the two pinned sources to a temporary tree and runs `patch --fuzz=0` there. It checks the exact output hashes and rechecks that the source files are unchanged. It also checks the centered eigenfunction by literal arc flips on every labelled tournament of orders 3 through 6, the uniform mean, and the minimal raw-cyclicity hostile at n=3.

Normal and optimized runs passed **33,881 explicit checks** with byte-identical receipts. This adds finite controls for the displayed normalization; no numeric count or spectrum has changed. The separate historical SC warning is still the already-audited `session52-counting-lenses/historical-warning.patch`, not part of this patch.

`prepare_patch.py` records the exact string replacements used to construct the diff. It is provenance tooling, not a required publication step; use the pinned patch rather than regenerate it against changed inputs.
