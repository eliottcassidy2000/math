"""Reproduce exact checks normally and with optimization; require identical outputs."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOBS = [
    'characters/check_characters.py',
    'root_deletion/check_deletion.py',
    'fano_codes/check.py',
    'support_reflections/verify.py',
    'audits/check_fano_independent.py',
    'audits/check_root_independent.py',
    'audits/check_character_support_independent.py',
]


def main():
    rows = []
    for name in JOBS:
        script = HERE/name
        runs = []
        for optimized in (False, True):
            args = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(script)]
            result = subprocess.run(args, cwd=ROOT, capture_output=True, check=False)
            if result.returncode:
                raise RuntimeError(f'{name}: exit {result.returncode}\n{result.stderr.decode(errors="replace")}')
            runs.append(result.stdout)
        if runs[0] != runs[1]:
            raise RuntimeError(f'{name}: optimized output differs')
        row = {'script': name, 'source_sha256': hashlib.sha256(script.read_bytes()).hexdigest(),
               'normal_optimized_equal': True, 'stdout_sha256': hashlib.sha256(runs[0]).hexdigest(),
               'stdout_bytes': len(runs[0]), 'return_codes': [0, 0]}
        rows.append(row)
        print(f'PASS {name}', flush=True)
    receipt = {'status': 'PASS', 'python': sys.version.split()[0], 'jobs': rows,
               'scope': 'exact finite replays plus separately written mathematical proof audits; no Lean claim'}
    (HERE/'REPLAY_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
