# Random archive walk: three exact tools for absolute window placement

I sampled nine Markdown reflections without replacement from the sorted path list under 07-reflections in archive snapshot 1b91742ce55a202b73c80d56561ba56d00b503a6, using Python Random(75466). Additional independent samples and seeds are in the two proof files. The root draw, in order, was:

1. lrc14-haar-product-discrepancy-tournament-tiling-codex-s165
2. the-safe-routes-reformulate-G-they-do-not-reduce-it-opus-S114
3. labelled-response-polynomials-gap-graphs-and-boundary-targets-20260825
4. the-completion-identity-is-bilinear-the-residual-is-multilinear-opus-S262
5. the-honest-state-of-LRC14-covering-is-a-reframing-decorrelation-is-the-core-macmini-S37
6. the-tail-floor-re-derived-on-the-longest-AP-axis-opus-S156
7. tiling-hamiltonian-ratios-s95
8. what-can-1over25-mean-a-tested-hypothesis-battery-klein-S413
9. lrc14-fixed-colour-long-row-capture-and-gauge-thickness-opus-20260725

What came out of the walk:

- **[A two-flank physical unit window](FLANKS.md).** A parity-balanced separator permits two independently chosen binary flanks without disturbing the exact +2 core repair. The new exact inverse in [inverse.py](inverse.py) selects source addresses from the target; it finds the same addresses (3533710,0) for 1000000007 and 1000000009, with different core contacts. It also decides when the family does not land. Rather than round a free phase to a count, it keeps both binary carry and literal equality. [Independent formula replay](replay_flanks.py).
- **[A physical radix and deterministic carry](SPIGOT.md).** Inspired by the sampled covering-certificate digit emitter, two small real tiles give a base-128 weighted digit recipe with a known carry at each length. The absolute target determines its only possible length; modular peeling then either constructs a physical source or certifies there is none in the family. It also computes exactly the neighboring reachable values, and the count of reachable values grows only like X^(1/7). A new alphabet-size obstruction explains why adding more modes with the same fixed profile does not yield universal coverage. [Replay](replay_spigot.py).
- **[Variable-size digits](STRETCH.md).** An extension chooses both the digit and its physical stretch by taking a 2-adic valuation of the exact target read. It increases the palette while keeping an exact inverse, but remains sparse.
- **[A Haar determinant repair](HAAR.md).** Two pieces have indistinguishable margins but different mixed boundary products. Reordering their marked copies moves the exact count by two at every return height. This is the physical counterpart of the selected old 2×2 Haar checkerboard, not just a metaphor. [Replay](replay_haar.py).

The proofs distinguish an absolute membership algorithm from a universal absolute placement theorem. The former now works on two substantially different families. Neither has been shown to contain every prescribed admissible target; the radix family is proved too sparse to do so.
