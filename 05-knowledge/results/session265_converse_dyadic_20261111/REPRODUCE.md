# Reproduction and status ledger

Run from the repository root with Python3. No nonstandard packages are needed.

    python3 -B 05-knowledge/results/session265_converse_dyadic_20261111/shadows/replay.py
    python3 -B 05-knowledge/results/session265_converse_dyadic_20261111/doubling/construct.py
    python3 -B 05-knowledge/results/session265_converse_dyadic_20261111/collatz/controls.py
    python3 -B 05-knowledge/results/session265_converse_dyadic_20261111/geometry/replay.py
    python3 -B 05-knowledge/results/session265_converse_dyadic_20261111/bridge.py
    python3 agents/check_docs.py

The first three programs regenerate their corresponding JSON receipts. The
last two print the deterministic `geometry/receipt.json` and
`BRIDGE_RECEIPT.json` contents. The same commands with `-O` reproduce the
matching outputs; the ordinary runs retain all assertion checks.

## Exact finite universes

- Shadows: all1098 labelled tournaments at orders2..5,48 seeded cases at6..9,
  and two six-vertex witnesses. All cuts, interfaces and selective-apex
  responses are checked. Witness extension counts use independent permutation
  enumeration. The Hamilton-cycle cut strengthening has additional small
  cycle-count controls.
- Fano switch: all128 line tables; labelled root/shadow certificates on seven
  vertices;121 small left/right context pairs; literal subset counts for all
  four17-vertex common-inventory outputs, independently matched to matrices.
- Source controller:28 full dyadic blocks through2^14;216 compiled targets
  through2^18; literal source recounts through14 vertices; seven prescribed
  Collatz depths;700 carry cases;30 prime-power top-digit controls;16 mixed
  3/7 controls; full residue orbits modulo7,21,63 with686,18,522,500,094 states.
  These are finite source-state universes, not isomorphism-class censuses.
- Electrical shadow:1354 tournaments, comprising all labelled orders2..5
  and128 seeded cases at each of6,7. Exact Bareiss tree counts check the floor
  and equality boundary. The spectral theorem is proved, not inferred from
  floating eigenvalues.
- Cross-lane bridge: eight tail/wrapper outputs at9..12 vertices with exact
  gate histograms; the strong ten-vertex1769/1771 pair is independently counted
  both by subset dynamic programming and actual path backtracking.

## Proof status

The general SCC-cut, electrical, binomial-amplifier and physical-reader
statements have direct proofs. Fano's specific input table is FINITE-EXACT;
the all-context sorting theorem follows algebraically from that certificate
and the inherited Session264 fusion law. The dyadic and prime-power theorems
are all-order proofs conditional only on the inherited physical transfer law.
The mixed3^a7^b theorem additionally uses the complete finite modulus21 base.
None of the all-order conclusions is inferred merely from a bounded census.

No mathematical claim from the supplied model texts or the pending THM-5133
draft is required. Old7/21 exclusions, Mersenne realizations, Camion's theorem,
odd-cycle counting and the established operation laws are accepted inputs.

## Publication and integrity

Checkpoint1: `17a27f63bbcaba9525e8ab7147ecc8267114662a`.
Checkpoint2: `19f07b39ae778092cfea0898c089e34e861f0bc2`.
Each remote commit's entire Git tree matched the corresponding verified local
tree before proceeding. Final synthesis and the explicit step2 reader are a
third checkpoint. Concurrent main changes were preserved; no force push.

`FINAL_CHECKS.json` records the final reproduction results. `MANIFEST.json`
hashes every delivered file in this packet except itself, plus its navigation
page. A final fetched-main comparison is required before reporting completion.
