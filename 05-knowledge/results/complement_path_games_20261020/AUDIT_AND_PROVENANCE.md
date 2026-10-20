# Publication review and reproducibility

Reviewed 20 October 2026. The complete proof was independently reviewed for the four game thresholds, prescribed-start quantifiers, rooted join criteria, sharp red-edge costs (including orders six and eight), and the connected-both construction at every order at least seven. The original proof and verifier are preserved byte for byte. The README's original “Independent review pending” sentence records its pre-review state; this receipt records the subsequent review.

The connected-both construction was supplied during review and independently checked by the producer. Its minimal order seven combines an all-order construction with exhaustive exclusion through order six. No classification of all connected nontraceable complement pairs is claimed.

## Reproduction

From this directory, run:

```sh
python3 verify_path_games.py > replay.stdout
cmp replay.stdout verification.stdout
sed 's@session16-path-games/@@' SHA256SUMS | sha256sum -c -
```

The archived SHA256SUMS uses the original source-directory prefix; the command above removes only that prefix. The verifier writes verification.json beside itself. It requires only Python's standard library.

A copied-layout normal run passed all assertions and matched both archived outputs exactly. A separate `python3 -O verify_path_games.py` run also reproduced the same output. Optimized mode disables Python assertions; it is a portability/output control, not a substitute for the normal-mode validation.

The checks cover all 33,867 labelled simple graphs on one through six vertices, with a separate permutation endpoint calculation for every graph. They include 6,391 disconnected-red join cases and 37,982 prescribed-root checks. Witness templates are checked through order sixteen. Both archived outputs and both copied-layout outputs have SHA-256:

`b228e0395ecb8c7d301308c5c5fcdebbb2a5f311018eb82a48ce7ca83a75f5d0`

Original proof SHA-256: `8229d55ca56569c92d01d2cacf11836ee628453e794297ab648077007fc29d18`.
Verifier SHA-256: `6f98cf8ac2f2bea17fcd7c4f6e348e7f3b06efd482d03ea8954945cc32e5c2f4`.

## Provenance and scope

The primary Pósa and Balister–Győri–Lehel–Schelp PDF links in the [proof](README.md) were checked and resolve to the named papers. All-order arguments are included in the proof; no novelty claim is made for classical extremal facts.

Repository comparison routes are [THM-3121, path-cover substitution](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md), [THM-082, Hamilton-path edge use](../../../01-canon/theorems/THM-082-deletion-contraction-ham-paths.md), and [HYP-2202, the orientation warning](../../hypotheses/HYP-2202-unit-distance-tournament-hampath-flop.md). They are comparisons with explicit retained/lost data, not claims that undirected graphs are tournaments.

Concurrent repository work independently records the free-start thresholds and sharp budgets in [the complement-game note](../kernels_games_detours_20261020f/game.md), published at commit c402d288ac84f695e0fb355444ee67bfb7917d4f. It is preserved unchanged. The present packet also exposes the every-start quantifier, exact rooted join criteria, and the connected-both boundary.

No result here closes the tournament raw-count spectrum, black metagraph connectivity, or a torus/game equivalence.
