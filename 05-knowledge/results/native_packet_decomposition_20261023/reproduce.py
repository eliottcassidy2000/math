"""Replay all packet certificates in normal and optimized Python."""
import hashlib
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def main():
    for relative in ('interface/check_filters.py', 'algebra/check.py',
                     'workspace/check.py', 'native/check.py', 'native/four_move_invariant.py'):
        outputs = []
        for optimized in (False, True):
            cmd = [sys.executable] + (['-O'] if optimized else []) + ['-B', str(HERE/relative)]
            result = subprocess.run(cmd, cwd=ROOT, text=True, encoding='utf-8',
                                    capture_output=True, timeout=240)
            if result.returncode:
                raise RuntimeError(f'{relative} failed ({result.returncode}):\n'
                                   + result.stdout + result.stderr)
            outputs.append(result.stdout.replace('\r\n', '\n'))
        if outputs[0] != outputs[1]:
            raise RuntimeError(f'normal/optimized output differs: {relative}')
        digest = hashlib.sha256(outputs[0].encode('utf-8')).hexdigest()
        print(f'PASS normal/-O {relative}; stdout SHA256 {digest}', flush=True)
    print('PASS all five frozen certificates and all five optimized replays')


if __name__ == '__main__':
    main()
