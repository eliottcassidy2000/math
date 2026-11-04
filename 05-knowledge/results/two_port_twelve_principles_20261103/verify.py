"""Verify frozen evidence; optional replays write only to a temporary directory."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ENGINES = [('check.py', 'check'), ('scout.py', 'cover_scout'),
           ('structure.py', 'structure'), ('run_twelve.py', 'twelve'), ('audit.py', 'audit')]


def need(ok, what):
    if not ok:
        raise RuntimeError(what)


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'manifest.json').read_text(encoding='utf8'))
    actual_files = {p.name for p in HERE.iterdir() if p.is_file() and p.name != 'manifest.json'}
    need(actual_files == set(manifest['files']), 'packet inventory')
    for name, digest in manifest['files'].items():
        need(sha(HERE/name) == digest, 'packet hash: '+name)
    for name, digest in manifest['dependencies'].items():
        need(sha(ROOT/name) == digest, 'dependency hash: '+name)
    for script, stem in ENGINES:
        receipt = HERE/(stem+'.json')
        need(receipt.read_bytes() == (HERE/(stem+'.optimized.json')).read_bytes(),
             'normal/optimized receipt: '+stem)
        data = json.loads(receipt.read_text(encoding='utf8'))
        pins = data['source_sha256']
        if isinstance(pins, str):
            pins = {script: pins}
        for name, digest in pins.items():
            need(sha(HERE/name) == digest, 'receipt source: '+name)
    bank = ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    for stem in ('cover_scout', 'twelve'):
        data = json.loads((HERE/(stem+'.json')).read_text(encoding='utf8'))
        need(sha(bank) == data['bank_sha256'], 'representative bank')
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='twelve-replay-') as tmp:
            for script, stem in ENGINES:
                for optimized in (False, True):
                    print('Replaying', script, '-O' if optimized else 'normal', flush=True)
                    output = Path(tmp)/(stem+('.optimized' if optimized else '')+'.json')
                    command = [sys.executable]+(['-O'] if optimized else [])+[
                        str(HERE/script), '--output', str(output)]
                    subprocess.run(command, check=True)
                    need(output.read_bytes() == (HERE/(stem+'.json')).read_bytes(), 'replay: '+stem)
    print(json.dumps(dict(status='PASS', frozen_files=len(manifest['files']),
                          dependencies=len(manifest['dependencies']), receipt_pairs=len(ENGINES),
                          replayed=args.replay,
                          scope='Scoped T1-T5 algebra and finite families; universal reconstruction OPEN'),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
