# Publication amendment: portable authored replay versus historical source audit

This **additive amendment** does not alter the original `AUDIT.md`,
`independent_audit.py`, `normal.json`, `optimized.json`, or `manifest.json`.
Those retain their original historical meaning and hashes.

## Licensing decision and exact scope

The complete, untruncated 217-entry tree of `calr3/tinycryptids` at commit
`3c43f97768d4ea06d70e0eb154ba88060bc36dec` contains no path named as a license,
COPYING, copyright, or NOTICE file. The root README at blob
`05c6451c38f16563b16784a39c74f9932db2f1cf` contains no licensing grant.
This is a bounded source inspection, not a legal conclusion that no permission
could exist elsewhere. **No external source redistribution is authorized by this
audit. Do not vendor any tinycryptids source, binary, or the Rust excerpt.**
See `portability_license_check.json` for the inspection record.

The original `independent_audit.py` depends on the external local source
snapshot, including `tinycryptids/source-manifest.json`, the 15 targeted Lean
import-closure files, and root `TinyCryptids.lean`. Its original 30 local blob
checks also inspect whichever manifested source snapshots are present. Therefore
it is retained as **historical audit code, not a portable publication replay
target**. The original source observations are preserved in `normal.json`,
`optimized.json`, `manifest.json`, `remote_git_pin.json`, and `AUDIT.md`.
Retaining and hash-checking those receipts is not rerunning the source audit.

## New portable checker

`portable_math_audit.py` is a separate authored checker. It keeps the original
arithmetic, orbit, divisor-sum, cylinder, periodic-extension, orientation,
producer-output, optimization-gate, and receipt-hash checks. It also compiles
and actually runs the independently authored `native_mb.c` in a temporary
directory on every invocation, verifies all native controls, and compares its
stdout byte-for-byte with the archived native JSON. No downloaded code runs.
Python's standard library and a local `gcc` are required; an unavailable compiler
is an explicit failure, not a silently skipped check.

The following are **not re-performed by the portable checker**:

1. Git-blob hashing of the 30 downloaded external source files;
2. token scanning of the 15-file Lean closure;
3. resolving local Lean imports and inspecting the root import file;
4. the GitHub remote-tree request, earlier external Python/C tests, or any Lean
   or Rust build.

These remain clearly labeled historical observations. The portable checker
verifies their archived receipt hashes and internal pin consistency, then emits
`HISTORICAL_RECEIPTS_HASH_VERIFIED_NOT_SOURCE_CHECKS_REPERFORMED` for this portion.
It does not emit a new full-source-audit PASS. Seven unrelated root Lean imports
were absent from the original partial snapshot; no full build was performed.

## Exact publication allowlist

`PORTABILITY_MANIFEST.json` is the machine-readable exact allowlist and SHA-256
seal: copy each entry's `source` to the identical relative payload path, plus
the manifest itself. Entries comprise only authored files, their stored outputs,
metadata/pins, and the original producer's five sealed authored payload files.
Do not recursively copy either input directory. No external source is included.
Do not place the original `independent_audit.py` in an executable replay task.

The sibling directory names must remain `session65-audit-weakcollatz` and
`session65-weakcollatz`, as the scripts resolve inputs relative to `__file__`.

From their common parent, run:

```sh
python3 -I session65-audit-weakcollatz/portable_math_audit.py > /tmp/s65-portable-normal.json
python3 -I -O session65-audit-weakcollatz/portable_math_audit.py > /tmp/s65-portable-optimized.json
cmp /tmp/s65-portable-normal.json /tmp/s65-portable-optimized.json
cmp /tmp/s65-portable-normal.json session65-audit-weakcollatz/portable_normal.json
python3 -I session65-weakcollatz/verify_bridge.py > /tmp/s65-producer-normal.json
python3 -I -O session65-weakcollatz/verify_bridge.py > /tmp/s65-producer-optimized.json
cmp /tmp/s65-producer-normal.json /tmp/s65-producer-optimized.json
cmp /tmp/s65-producer-normal.json session65-weakcollatz/bridge_results.json
```

The C replay is embedded in each portable Python run. It is not necessary to
publish binaries or separately execute the archived external simulator. All
portable output is deterministic. Both modes also pass in a frozen temporary
closure containing only the allowlisted authored inputs, with no tinycryptids
directory: see `portable_closure_receipt.json`.

This is a portability repair, not new mathematics or an upgraded formalization
claim. All-seed halting, ordinary Collatz convergence, Lean/Rust build success,
and new Busy Beaver maximality remain unclaimed.
