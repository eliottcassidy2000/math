"""Write or check hashes of the complete research packet, excluding caches."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    files = {}
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name == 'manifest.json':
            continue
        if p.suffix not in ('.py', '.md', '.json', '.cpp', '.txt'):
            raise RuntimeError('Unexpected packet artifact: ' + str(p))
        files[str(p.relative_to(HERE))] = hashlib.sha256(p.read_bytes()).hexdigest()
    record = dict(scope='All packet artifacts except this manifest and Python caches', files=files)
    path = HERE / 'manifest.json'
    if args.write:
        path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    elif json.loads(path.read_text()) != record:
        raise RuntimeError('Packet manifest differs from current artifacts')
    print(json.dumps(dict(status='PASS', artifact_files=len(files), mode='write' if args.write else 'check')))


if __name__ == '__main__':
    main()
