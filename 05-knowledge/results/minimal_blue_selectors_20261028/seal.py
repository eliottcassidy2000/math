"""Explicit maintainer step after reviewing edits; never called by verify."""
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
    previous=HERE.parent/'singleton_pair_reconstruction_20261028'
    pins=json.loads((previous/'sources.json').read_text(encoding='utf-8'))
    for name,wanted in pins.items():
        if digest(ROOT/name)!=wanted:raise RuntimeError(('inherited source changed; audit before resealing',name))
    for name in ('reconstruct.py','boundary_profile.py','rooted-seven.json','THEORY.md'):
        path=previous/name;pins[path.relative_to(ROOT).as_posix()]=digest(path)
    (HERE/'sources.json').write_text(json.dumps(dict(sorted(pins.items())),indent=2)+'\n',encoding='utf-8')
    exclude={'manifest.json','verification.json','verification.optimized.json'}
    manifest={path.name:digest(path) for path in sorted(HERE.iterdir()) if path.is_file() and path.name not in exclude}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(pinned_sources=len(pins),sealed_artifacts=len(manifest))))


if __name__=='__main__':main()
