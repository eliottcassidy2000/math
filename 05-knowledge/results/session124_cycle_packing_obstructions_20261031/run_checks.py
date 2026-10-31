"""Portable exact replay of all four checkers in normal and optimized mode."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SCRIPTS = ('check', 'pattern_check', 'check_witnesses', 'independent_check')

def verify(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = json.loads((HERE / 'BASELINE_RECEIPTS.json').read_text())['sha256']
    expected = {f'{script}.{mode}.json' for script in SCRIPTS
                for mode in ('normal', 'optimized')}
    if set(baseline) != expected:
        raise RuntimeError('Baseline filename set is not exactly the eight expected receipts')
    rows = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for script in SCRIPTS:
        outputs = []
        for mode in ('normal', 'optimized'):
            command = [sys.executable]
            if mode == 'optimized':
                command.append('-O')
            command.append(str(HERE / (script + '.py')))
            completed = subprocess.run(command, cwd=output_dir, env=env,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       check=False)
            if completed.returncode:
                raise RuntimeError(f'{script}/{mode} failed:\n' + completed.stderr.decode())
            data = completed.stdout
            parsed = json.loads(data)
            if parsed.get('status') != 'PASS':
                raise RuntimeError(f'{script}/{mode} did not return PASS')
            name = f'{script}.{mode}.json'
            digest = hashlib.sha256(data).hexdigest()
            if digest != baseline[name]:
                raise RuntimeError(f'{name}: output differs from its audited baseline')
            (output_dir / name).write_bytes(data)
            outputs.append(data)
            rows.append({'receipt': name, 'checks': parsed['checks'],
                         'sha256': digest, 'baseline_identical': True})
        if outputs[0] != outputs[1]:
            raise RuntimeError(f'{script}: normal/optimized byte mismatch')
    result = {'status': 'PASS', 'programs': len(SCRIPTS), 'executions': len(rows),
              'normal_optimized_byte_identical': True, 'receipts': rows,
              'scope': 'restricted spectrum, finite forbidden language, and word/isomorphism bounds only'}
    (output_dir / 'verification.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        help='Retain receipts here; default is an automatically removed temporary directory.')
    args = parser.parse_args()
    if args.output_dir is None:
        with tempfile.TemporaryDirectory(prefix='cycle-packing-replay-') as tmp:
            result = verify(Path(tmp))
    else:
        result = verify(args.output_dir.expanduser().resolve())
    print(json.dumps(result, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
