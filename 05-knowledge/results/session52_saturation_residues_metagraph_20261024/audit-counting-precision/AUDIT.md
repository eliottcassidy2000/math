# Independent audit: narrow counting precision repair

Verdict: PASS for the pinned narrow edits, not a blanket certification of either historical theorem file.

## Candidate and inheritance

The exact audited patch is `inputs/producer/precision.patch`, SHA-256
`8880e6a357ecd0aa426378e6a9e2325d20c235124925eb9d05e9cddf66b587c0`.
All seven entries of the frozen producer manifest verify against the producer files. The two original theorem files are separately pinned by source SHA-256. Applying the patch with `patch --batch --fuzz=0 -p1` in temporary copies produces exactly both result hashes in `source-pins.json`.

The theorem IDs and source attribution, numeric tables, dependency IDs, related IDs and result paths are preserved. The HYP-3544 relationship remains present; its explanatory comment now properly marks the topological proposal as requiring additional structure. This is a precision repair to existing results, not two newly counted theorems, and does not replace or modify the independently audited counting-lenses packet. No repository or producer files were written.

## Mathematical review

THM-587 correctly replaces the average of fixed-subspace dimensions by the average of traces. For a permutation representation, a trace counts fixed basis vectors while its fixed-subspace dimension counts cycles of the permutation on the basis. Already at n=2 the nontrivial vertex permutation interchanges the two labelled tournaments, with trace zero but a one-dimensional fixed subspace. Averaging traces is the invariant-dimension projector formula; averaging fixed-subspace dimensions is not.

The complement involution acts on Fourier level k by (-1)^k and permutes tournament isomorphism classes. Its trace is therefore both the alternating sum of invariant level dimensions and the number of fixed classes. This is a graded antipodal trace. Without a chosen complex and differential it does not imply a topological Euler characteristic, Betti numbers or a Borsuk-Ulam obstruction. The repair states this boundary explicitly and retains the historical proposed direction rather than silently removing its attribution. Complement freeness on labelled tournaments requires n>=2 and does not imply freeness on isomorphism classes.

THM-588 correctly identifies c_0(T)=c(T)-binom(n,3)/4 as the adjacency eigenfunction. Each labelled triple is cyclic in two of its eight orientations. Thus the constant Fourier coefficient is binom(n,3)/4. Cyclicity has degree at most two; its invariant linear component vanishes by the arc transposition argument. The remaining component is nonzero for n>=3: a transitive tournament has no cyclic triple, while a tournament with one cyclic triple exists. Consequently c_0 is a nonzero level-two invariant, so A c_0=(d-4)c_0 and L c_0=4c_0. Uniform labelled weighting induces orbit-size weighting on classes, not uniform class weighting.

The independent checker counts cyclic triples directly from their three arc bits, unlike the producer's score-sequence formula. It checks every labelled tournament for n=3,4,5,6, respectively 8,64,1024,32768 states. The mean, nonzero centered mode and every one-arc-flip adjacency identity pass. The producer additionally checks that raw cyclicity fails the n=3 eigenfunction identity while its centered version passes. These finite checks support the elementary all-n argument; they are not a replacement for it.

## Reproducibility and dependency closure

Existing normal and optimized receipts were inspected and then independently reproduced from a real relocation. `relocation-manifest.json` lists the exact 11 files copied: the independent checker, the frozen producer manifest and its seven entries, and the two pinned original theorem files. No symlinks or undeclared files are copied. Execution takes place from `/tmp` using isolated Python `-I`, with and without `-O`, and with `PYTHONPATH` removed. Only Python's standard library and the system `patch` executable are required. `prepare_patch.py` is retained as provenance, never executed; its historical repository-relative regeneration path is not a runtime dependency.

The producer checker passes 33,881 active gates in each mode; its receipts are byte-identical with SHA-256 `6eada2ebabf8239e9a7fe41ac70fb410ce8cbe32fb976f30c5676ff0dd10a1fc`.
The independent checker passes 33,895 active gates in each mode; its receipts are byte-identical with SHA-256 `f4655f429684a9655ea98492381f462ea567b0e5ad1bbc61bc704ba44a550706`.

A modified pinned source is rejected in both modes before a success receipt is emitted. All original audit inputs retain their hashes after replay. `verify_relocation.py` itself was also run normally and optimized, producing byte-identical relocation receipts.

Reproduce with:

    python3 -B verify_relocation.py > relocation-receipt.json
    python3 -B -O verify_relocation.py > relocation-receipt-optimized.json
    cmp relocation-receipt.json relocation-receipt-optimized.json

## Limits and publication gate

This audit approves only the narrow diff. It does not revalidate untouched historical mixing-time, engineering, correlation or performance claims; it does not rerun unavailable historical scripts; and it does not revise any historical counts. A separate historical self-converse warning patch belongs to counting-lenses and is outside this patch.

Before applying to a publication checkout, verify both original source hashes. Stop on any mismatch rather than fuzz-applying or treating an already changed target as equivalent. Verify both exact output hashes afterward. No commit, push, publication, or producer mutation was performed by this audit.
