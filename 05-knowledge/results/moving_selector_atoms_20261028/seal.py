"""Explicit maintainer sealing step; verify never invokes this module."""
from pathlib import Path
from hashlib import sha256
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):
    data=path.read_bytes()
    if path.suffix!='.gz':data=data.replace(b'\r\n',b'\n')
    return sha256(data).hexdigest()


def main():
    previous=HERE.parent/'minimal_blue_selectors_20261028'
    pins=json.loads((previous/'sources.json').read_text(encoding='utf-8'))
    for name,wanted in pins.items():
        if digest(ROOT/name)!=wanted:raise RuntimeError(('inherited source changed; audit first',name))
    extra=['minimal_blue_selectors_20261028/blue_minimum.py','minimal_blue_selectors_20261028/THEORY.md',
           'strong_atom_continuation_vitali_20261026/strong-atlas.json','strong_atom_continuation_vitali_20261026/README.md',
           'session90_even_graphs_odd_descent_20261029/oddness/check_centralizer_path_dp.py',
           'session90_even_graphs_odd_descent_20261029/oddness/check_fixed_orbit_descent.py',
           'session90_even_graphs_odd_descent_20261029/trijection/check_tree_compression.py',
           'session90_even_graphs_odd_descent_20261029/trijection/SPARSE_EVEN_CARRIERS.md']
    for name in extra:
        path=HERE.parent/name;pins[path.relative_to(ROOT).as_posix()]=digest(path)
    (HERE/'sources.json').write_text(json.dumps(dict(sorted(pins.items())),indent=2)+'\n',encoding='utf-8')
    excluded={'manifest.json','verification.json','verification.optimized.json'}
    manifest={path.name:digest(path) for path in sorted(HERE.iterdir()) if path.is_file() and path.name not in excluded}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(pinned_sources=len(pins),sealed_artifacts=len(manifest))))


if __name__=='__main__':main()
