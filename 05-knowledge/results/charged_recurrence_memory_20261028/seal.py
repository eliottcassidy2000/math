"""Pin the packet, executable imports, and principal mathematical dependencies."""
from pathlib import Path
import hashlib
import json
import sys
import charged, jets, memory, renewal, frobenius, rooks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


if '--write' in sys.argv:
    files = {p for p in HERE.iterdir() if p.is_file() and p.name != 'seal.json'}
    for module in list(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            p = Path(name).resolve()
            if p.is_relative_to(ROOT) and p.suffix == '.py':
                files.add(p)
    names = [
        '05-knowledge/results/complement_port_crt_20261028/PORT_LAWS.md',
        '05-knowledge/results/complement_port_crt_20261028/CLOSURE.md',
        '05-knowledge/results/boundary_support_20261020_metagraph.md',
        '05-knowledge/results/session95_niche_recurrences_20261029/jets/POLE_RECURRENCE_MEMORY.md',
        '05-knowledge/results/session95_niche_recurrences_20261029/flows/SCHUR_BOUNDARY_RECURRENCES.md',
        '05-knowledge/results/stabilizer_necklace_kernels_20261028/KERNELS.md',
        '05-knowledge/results/affine_tree_recurrences_20261028/README.md',
        '01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md',
        '01-canon/theorems/THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure.md',
        '01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md',
        '01-canon/theorems/THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md',
        '01-canon/theorems/THM-4225-bad-owner-upper-zeta-successor-rook-hierarchy.md',
        '01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md',
    ]
    files.update(ROOT / name for name in names)
    data = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)}
    (HERE / 'seal.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    print('Sealed', len(data), 'files')
else:
    data = json.loads((HERE / 'seal.json').read_text(encoding='utf-8'))
    bad = [name for name, value in data.items() if not (ROOT / name).is_file() or sha(ROOT / name) != value]
    if bad:
        raise RuntimeError('Changed sealed files: ' + ', '.join(bad))
    print('Verified', len(data), 'file hashes')
