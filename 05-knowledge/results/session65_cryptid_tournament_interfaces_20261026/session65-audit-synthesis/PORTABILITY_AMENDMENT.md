# Session 65 synthesis audit: portable publication amendment

2026-10-26. **PASS for portable replay; no mathematical amendment.**

The historical `REPORT.md`, `independent_check.py`, `SOURCE_PINS.json`, normal
and optimized receipts, `MANIFEST.json`, and `SHA256SUMS` remain unchanged.
They record the original local audit, where the scan was present. The original
checker is intentionally retained as historical evidence, not advertised as
runnable from the source-free publication bundle. Use the new entrypoint below.

## 1. Exactly one external nonruntime observation

The published packet does not redistribute this historical scan:

- Original local path: `session65-fifteens/sloane1973p98.jpg`.
- Source URL: <https://neilsloane.com/doc/persistence.2.jpg>.
- Recorded SHA-256:
  `2e2b51802a88dc7e3c785237f58605f06db5f7a8eb424bb465ed79f901071328`.
- Classification: **EXTERNAL_NONRUNTIME**.

The earlier fifteens audit inspected the scan. The portable checker does **not**
read, copy, download, display, or locally reverify it. Its receipt explicitly
sets `locally_verified: false`. The URL and digest identify a prior observation;
they are not current file-verification or redistribution evidence. No mathematical
function in the synthesis audit reads the image or depends on its pixels.

`EXTERNAL_NONRUNTIME_INPUTS.json` contains precisely this declaration. The
entrypoint hard-codes the exact declaration and validates its sealed local copy.
Its inherited-seal verifier recognizes only the exact path **and** expected hash.
There is no wildcard exclusion, optional missing-file branch, URL fetch, or
permission to externalize any other missing dependency. A changed historical
manifest or added exclusion is rejected.

## 2. Full local closure, including the authored bridge receipt

`PORTABLE_INPUTS.json` is the explicit 174-file local dependency allowlist.
The entrypoint pins its digest, then checks every listed path's presence,
byte length and SHA-256. It retains the 88 historical source-pin checks and
229 local inherited-seal checks. The previous count of 230 included one image
observation; the new output reports that one separately rather than claiming
230 current local verifications. The 174, 88 and 229 counts overlap: they are
not counts to add together.

The formerly unlisted file
`session65-weakcollatz/bridge_results_recheck.json` is **required and explicitly
approved for publication**, not skipped:

- Size: **1691 bytes**.
- SHA-256:
  `65f0489ce2c16784dae5866f0c36be27e960dbb9174fcd6d76cebd316e5a2bfc`.
- It equals the authored producer's normal and optimized result bytes.
- A frozen-package run of the inspected, authored `verify_bridge.py` in both
  modes reproduces it exactly. No external simulator or native long run is
  involved in that deterministic arithmetic/residue checker.

`PORTABILITY_MANIFEST.json` explicitly includes this additional source path,
the portable additions, and the entrypoint command. Other dependencies remain
covered by the existing producer/audit allowlists plus `PORTABLE_INPUTS.json`.
No repository file or historical seal was edited by this lane.

## 3. Mathematical body and consequences unchanged

`PORTABILITY_DIFF.patch` records the complete entrypoint difference. Changes
are confined to the runtime provenance/dependency gate, one tightly specified
external-observation classification, and receipt metadata. The entire body
from complete-carrier generation through the three action-result constructions
is byte-identical to the historical checker, including every mathematical
function. No pair orbit, orientation coefficient, finite geometry claim,
minimizer row, or hostile mathematical consequence was changed.

The frozen portable output's `carriers`, `results` and `scope` equal the
historical JSON exactly. Its normal and optimized outputs are byte-identical.
It therefore retains complete verification of the duad/syntheme polynomials
and actual literal minimizers, without rerunning the unrelated long audits.

## 4. Frozen closure and hostile verification

`portable_hostiles.py` creates a fresh temporary tree containing only the
174 explicitly approved local dependencies and the portable entrypoint/pin
manifest. The historical image, unshipped tinycryptids source directory,
original repository checkout and Git metadata are absent. Subprocesses use
isolated Python (`-I`), no ambient `PYTHONPATH`, and `PATH=/nonexistent`.
Licensed snapshots already in the approved evidence closure remain inert data;
no external source is executed or built.

Executed controls:

- Frozen normal and `-O` replays pass, agree byte-for-byte, and equal the
  stored portable receipt.
- Every one of the 174 local dependencies is separately removed and separately
  changed by a same-length first-byte mutation (or a nonempty mutation for an
  empty file), in both modes: **696 rejected hostile runs**.
- Two attempts to declare the missing bridge recheck as another external input,
  one through each declaration/allowlist file, are rejected in both modes:
  **4 additional rejected hostile runs**.
- Restoring every file restores the original successful result.
- The mathematical source-body equality and historical consequence equality
  are checked explicitly.
- Both frozen authored bridge replays reproduce the newly allowlisted receipt.
- The hostile runner itself also produces byte-identical output when launched
  normally and under `-O`.

## 5. Reproduction and seal

From a source-free packet root with the sibling layout preserved:

```sh
python3 -I -B session65-audit-synthesis/portable_check.py > portable-normal.json
python3 -I -B -O session65-audit-synthesis/portable_check.py > portable-optimized.json
cmp portable-normal.json portable-optimized.json
cmp portable-normal.json session65-audit-synthesis/portable-normal.json
python3 -I -B session65-audit-synthesis/portable_hostiles.py > portable-hostiles.json
python3 -I -B -O session65-audit-synthesis/portable_hostiles.py > portable-hostiles-optimized.json
cmp portable-hostiles.json portable-hostiles-optimized.json
```

The separately sealed amendment is `PORTABILITY_MANIFEST.json`; its companion
`PORTABILITY_SHA256SUMS` covers only new amendment files and does not replace
historical `SHA256SUMS`. All new code was read before execution. No copyrighted
image was copied, no external code or formal build was run, and no new research
claim is made. No permission denial or unresolved execution blocker occurred.
