# Fixed-bank joins: an all-length finite-state obstruction

Chan–Pak completeness and the conjectured all-odd coverage by F2 remain OPEN.

`PROOF.md` proves a sharp restriction on constructive certificates: every finite-state transducer selecting ordinary binary witnesses (v,N) with N in a fixed finite-bank family a g(v)−b e(v)+c(v+1), a,c>0, can select only polylogarithmically many outputs up to X. It may be asynchronous, nondeterministic, partial, and omit the bank label. The theorem does not concern output-only recognition or unbounded-register algorithms.

Mechanism: repeating any mixed binary period adds a strictly positive linear drift to g, in addition to its exponential term. A pumped finite-state output has no such linear term. Consequently every useful input loop is monochromatic, every accepted witness has bounded run count, and output support is sparse.

`PROOF.md` also records the inherited Session216 curvature obstruction specialized to F2, an exact obstruction to uniform target-halving within the physical rooted bank: suffix substitution retains the odd g coefficient, which cannot become twice another odd physical coefficient.

Run `python3 verify.py`. The retained `RESULT.json` reports 24,312 exact recurrence checks across all 2,026 mixed words of lengths 2–10 and 100 alternating-word checks. These are regressions for the all-length elementary proofs, not a new coverage scan. No particular bank file is needed by the theorem or verifier.

Independent review is maintained by the session's audit worker. No shared repository file has been modified by this worker and no theorem ID has been reserved.
