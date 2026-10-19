# Three-color graceful release: exact bridge criterion

This package proves a connected-buffer release dichotomy for nonempty graceful completion fibers with three free difference colors. It also proves a general tight-kernel and tight-chain sufficient condition, using the correctly attributed strict-Rado connectivity theorem of Fon-Der-Flaass (1996).

Read [the proof](THREE_COLOR_RELEASE_DICHOTOMY.md). Its first-checkpoint dependency is [the release and matching-support note](../dynamics_labeling_entry_20261019/GRACEFUL_RELEASE_AND_MATCHING_SUPPORT.md).

The new 15-label witness has a nonforced color-1 release that fails: the original two components of size 3 become two components of size 6. Other single releases repair this witness; the existential one-edge-release question remains open.

Replay with Python 3, standard library only:

    python verify_three_color_release.py
    python -O verify_three_color_release.py

Both modes produce exactly the included receipt. The verifier checks both arithmetic witnesses using every permitted integer-difference edge, all single-edge releases, rank/tight-kernel data, and the bridge-case component product. It does not replace the structural proofs. Exploratory abstract support censuses are excluded.
