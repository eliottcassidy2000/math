# Fresh570 independent integration audit

Read `AUDIT.md` for disposition and `PROOF.md` for the precise mathematics.
The incoming return passes. Apply the narrow2-primary wording correction
from `NATIVE570_WORDING.patch`; the full group has a quantitatively certified
nontrivial odd part. Native result credit remains with the incoming packet.

Self-contained reproduction from this directory, with no Git or network:

```sh
python3 -B replay_native.py native570 > producer-native.normal.json
python3 -O -B replay_native.py native570 > producer-native.optimized.json
python3 -B replay_native.py audit570 > producer-audit.normal.json
python3 -O -B replay_native.py audit570 > producer-audit.optimized.json
python3 -B check.py --out normal.json
python3 -O -B check.py --out optimized.json
python3 -B replay_guard_completeness.py --out guard-completeness.normal.json
python3 -O -B replay_guard_completeness.py --out guard-completeness.optimized.json
python3 -B check_odd_part.py --out odd-part.normal.json
python3 -O -B check_odd_part.py --out odd-part.optimized.json
cmp normal.json optimized.json
cmp producer-native.normal.json producer-native.optimized.json
cmp producer-audit.normal.json producer-audit.optimized.json
cmp guard-completeness.normal.json guard-completeness.optimized.json
cmp odd-part.normal.json odd-part.optimized.json
python3 -B verify_manifest.py
python3 -O -B verify_manifest.py
```

The full filtration uses only stdlib exact integer and finite-field operations.
It is retained separately so native-menu checks cannot be mistaken for proof
that a selected character bank is complete.

The auxiliary raw-domain extension is reproduced separately:

```sh
python3 -B check_helper_raw.py --out helper-raw.normal.json
python3 -O -B check_helper_raw.py --out helper-raw.optimized.json
cmp helper-raw.normal.json helper-raw.optimized.json
```

Section7 of `PROOF.md` explains why the actual raw2393-to570 alternative
cannot bypass the23-clock boundary under the specified margin-two baseline.

## Archival source portability

The five source documents under `inputs/` are byte-identical plain-text
`.txt` archives of original Markdown, not live relocated documents. Their
embedded relative links remain original source text. No archival byte or
source SHA changed. `SOURCE_PINS.json` maps the archival filenames to the
verified original repository commit URLs below. The frozen
`inputs/native/MANIFEST.json` remains an archival original manifest, so its
original `.md` filenames are intentionally not rewritten. The publisher
patch still correctly targets the live `NATIVE570.md`.

- `inputs/native/NATIVE570.txt`: [original NATIVE570.md](https://github.com/eliottcassidy2000/math/blob/617210b43a9015ec308e048eb2c7e07d3db9d99b/05-knowledge/results/sixslot570_20261024/NATIVE570.md).
- `inputs/native/PORTS.txt`: [original PORTS.md](https://github.com/eliottcassidy2000/math/blob/617210b43a9015ec308e048eb2c7e07d3db9d99b/05-knowledge/results/sixslot570_20261024/PORTS.md).
- `inputs/original/injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.txt`: [original ELECTRICAL_CODES.md](https://github.com/eliottcassidy2000/math/blob/cf87f98a581089b3f84d3fee16324195ab976e1e/05-knowledge/results/injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md).
- `inputs/original/metagraph_order8_electrical_20261023/README.txt`: [original README.md](https://github.com/eliottcassidy2000/math/blob/cf87f98a581089b3f84d3fee16324195ab976e1e/05-knowledge/results/metagraph_order8_electrical_20261023/README.md).
- `inputs/guard-theorem/AUDIT.txt`: [original AUDIT.md](https://github.com/eliottcassidy2000/math/blob/cf87f98a581089b3f84d3fee16324195ab976e1e/05-knowledge/results/session45_injective_multisource_moves_20261023/audit-dyadic/AUDIT.md).

## Standalone delivery check

Run `python3 -B check_portability.py` (or `python3 -O -B check_portability.py`)
to copy only manifest-listed files and the root manifest into a renamed
temporary package, then replay all six normal/-O mathematical receipt pairs
using isolated Python. The frozen `PORTABILITY_REPLAY.json` records this
check. The original nested `inputs/native/MANIFEST.json` is required input
data and is included in the delivery manifest. The verifier additionally
rejects any source-pin snapshot omitted from that manifest.
