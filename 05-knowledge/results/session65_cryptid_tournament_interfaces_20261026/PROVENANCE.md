# Provenance, source contexts and licensing

The approved lane/audit files were copied byte for byte. Their seals are
recorded in `FIRST_SOURCE_VERIFICATION.json`; the closed shipped file set is
`MANIFEST.json`. The initial mathematical context is repository commit
`6a8ad5d3bf01a232da4de41bcbff68390a9337ee`. Integration fast-forwarded cleanly
to `fcc112a69d7080db4cd619d6e7927be1992bc51e` before session-owned edits.
No initial source pin was silently replaced by a newer theorem.

The fifteen own-repository inputs and two runtime theorem snapshots have
byte-identical historical source checks in `DEPENDENCY_PINS.json`.
`SOURCE_LINK_CONTEXT.json` gives their original link contexts;
[SOURCE_LINK_INDEX.md](SOURCE_LINK_INDEX.md) supplies usable repository routes.
Historical documents retain their own scope/correction history and are not
promoted to new current claims merely by being included.

## External observations are not all bundled inputs

`session65-fifteens/source-manifest.json` records inspected inputs, not the
shipped-file manifest. Its `sloane1973p98.jpg` entry is deliberately omitted:
that copyrighted full-page scan is available at the recorded primary-source
URL and its inspected SHA256 is retained, but it is not redistributed.
No runtime checker needs that scan. Exploratory search inventories and caches
are also omitted. Thus not every provenance input is represented as bundled.

The single pinned Coq source snapshot
`session65-audit-cryptids/Bouncers.source.txt` is redistributed with the exact
`Bouncers.LICENSE.txt` copyright/permission notice (Copyright 2023 Maja
Kadziołka). The license blob is `015db163bb1176fc0b84e868b5030d06e88dacce`
at upstream commit `1c23203761602c199a7ac5ea6f7947fb80851d93`; source URL,
source hash and source-only audit limits are in `SOURCE_AUDIT.json`.
The license permits copying/distribution with its notice. No external Coq
code was executed or built, no Coq axiom/import audit was performed, and no
certificate for the three four-state machines is claimed.

External published theorems and verification bounds remain CITED, not newly
reproved or rerun. Community machine registries retain their independent-
verification qualifications. All five executed packet checkers are authored
Python with standard-library dependencies only.
