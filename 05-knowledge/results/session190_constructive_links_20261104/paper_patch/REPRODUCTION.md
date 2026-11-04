# Reproduce the constructive paper-patch packet

This packet preserves the reviewed source files and their original sibling-directory layout. The independent checker reads `../attached_papers/boundary_signature_demo.json`; no installed packages or startup helper are needed.

First, from the repository root, restore and verify all available checkpoint manifests in ascending order using the result root's `RESTORE_ARCHIVES.py`. Follow the current composite command in `../RECEIPTS.md`; on current main, using only CP01 would compare historical navigation/tool bytes against their newer versions. Then replay in a temporary copy, because the independent checker writes its receipt beside its script:

```sh
packet=05-knowledge/results/session190_constructive_links_20261104/paper_patch
work=$(mktemp -d)
cp -R "$packet" "$work/paper_patch"
python "$work/paper_patch/attached_papers/boundary_signature_demo.py" > "$work/demo-replayed.json"
python "$work/paper_patch/session190-papers-audit/independent_boundary_audit.py" > "$work/audit-replayed.log"
python - "$packet" "$work" <<'PY'
import json, sys
from pathlib import Path
packet, work = map(Path, sys.argv[1:])
original = json.loads((packet/'attached_papers/boundary_signature_demo.json').read_text())
replayed = json.loads((work/'demo-replayed.json').read_text())
if original != replayed:
    raise ValueError('Original demo receipt differs from replay')
original = json.loads((packet/'session190-papers-audit/independent_boundary_audit.json').read_text())
replayed = json.loads((work/'paper_patch/session190-papers-audit/independent_boundary_audit.json').read_text())
original.pop('input'); replayed.pop('input')
if original != replayed:
    raise ValueError('Independent mathematical receipt differs from replay')
print('Both mathematical receipts reproduce; retained files were not modified.')
PY
```

The independent receipt's informational `input` field contains the resolved local input path and consequently changes with the checkout location; this is not a mathematical-data difference. The original historical receipt stays untouched in the repository.

The independent implementation uses edge-subset enumeration and a disjoint-set cycle test, and imports no original demo code. The retained audit includes witnesses for the complete full and vertex-deletion state sets in these two finite patches. It does not certify the papers' full constructions, exhaustive global searches, or coloring configuration library.
