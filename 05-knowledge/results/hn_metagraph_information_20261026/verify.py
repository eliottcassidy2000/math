"""Verify frozen publication bytes; optionally replay into fresh temporary files."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def files():
    return {p.relative_to(HERE).as_posix(): p for p in HERE.rglob('*')
            if p.is_file() and p.name != 'MANIFEST.json' and '__pycache__' not in p.parts}


def run(mode, script, args=()):
    cmd = [sys.executable, '-B'] + (['-O'] if mode else []) + [str(HERE / script)] + list(args)
    result = subprocess.run(cmd, cwd=ROOT, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--full', action='store_true', help='also regenerate the order-nine bank in both modes')
    args = parser.parse_args()
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    found = files()
    need(set(found) == set(manifest['files']), 'manifest file inventory changed')
    for name, path in found.items():
        data = path.read_bytes()
        need(sha256(data).hexdigest() == manifest['files'][name]['sha256'], 'hash: ' + name)
        need(len(data) == manifest['files'][name]['bytes'], 'size: ' + name)
    pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for pin in pins['local_sources']:
        need(sha256((ROOT / pin['path']).read_bytes()).hexdigest() == pin['sha256'], 'source pin: ' + pin['path'])
    # Snapshot every expected result BEFORE any producer invocation.
    frozen = {name: path.read_bytes() for name, path in found.items()}
    print(f'PASS frozen inventory: {len(found)} files; {len(pins["local_sources"])} source pins', flush=True)
    if not (args.replay or args.full):
        return
    stdout_jobs = [
        ('metric_transport/check.py', 'metric_transport/receipt.json'),
        ('six_color/check.py', 'six_color/receipt.json'),
        ('information/check.py', 'information/receipt.json'),
        ('six_color/check_information_audit.py', 'six_color/information_audit_receipt.json'),
        ('information/check_six_color_audit.py', 'information/six_color_audit_receipt.json'),
        ('audits/check_information_independent.py', 'audits/information_independent.json'),
        ('audits/check_six_color_independent.py', 'audits/six_color_independent.json'),
        ('audits/check_enumeration_independent.py', 'audits/enumeration_independent.json')]
    with tempfile.TemporaryDirectory(prefix='hn-metagraph-replay-') as folder:
        temp = Path(folder)
        for mode in (False, True):
            tag = 'optimized' if mode else 'normal'
            for script, expected in stdout_jobs:
                need(json.loads(run(mode, script)) == json.loads(frozen[expected]), 'replay: ' + script)
                print('PASS ' + tag + ' ' + script, flush=True)
            output = temp / 'metric-audit.json'
            run(mode, 'enumeration/check_metric_independent.py', ['--output', str(output)])
            need(json.loads(output.read_bytes()) == json.loads(frozen['enumeration/metric-independent-normal.json']), 'metric audit replay')
            print('PASS ' + tag + ' independent metric audit', flush=True)
            output = temp / 'counts.json'
            run(mode, 'enumeration/check_counts.py', ['--output', str(output)])
            need(json.loads(output.read_bytes()) == json.loads(frozen['enumeration/counts-normal.json']), 'count engine replay')
            print('PASS ' + tag + ' count engine through128', flush=True)
            for n in ([8, 9] if args.full else [8]):
                output = temp / f'order{n}.json'
                bank = temp / (f'bank{n}.json.gz' if n == 9 else f'bank{n}.json')
                run(mode, 'enumeration/enumerate.py', ['--max-order', str(n), '--output', str(output), '--bank', str(bank)])
                saved = 'enumeration/order9-normal.json' if n == 9 else 'enumeration/normal.json'
                saved_bank = 'enumeration/representatives9.json.gz' if n == 9 else 'enumeration/representatives.json'
                need(json.loads(output.read_bytes()) == json.loads(frozen[saved]), 'family receipt replay')
                need(bank.read_bytes() == frozen[saved_bank], 'literal family bank replay')
                print('PASS ' + tag + f' literal family bank through{n}', flush=True)
    # Also detect any mutation of expected artifacts while replay was running.
    need(all(path.read_bytes() == frozen[name] for name, path in found.items()), 'publication bytes mutated during replay')
    print('PASS all requested replays; frozen publication bytes unchanged', flush=True)


if __name__ == '__main__':
    main()
