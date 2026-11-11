#!/usr/bin/env python3
"""Replay the new packet's finite controls; no inherited large census."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
CHECKS = (
    'actuators/compile.py',
    'history/probe.py',
    'history/port_probe.py',
    'history/native_controls.py',
    'prime_layers/compile.py',
    'prime_layers/periodic.py',
    'prime_layers/phase_clocks.py',
    'tournament/probe.py',
    'tournament/birth_search.py',
    'tournament/one_birth_close.py',
    'tournament/composite_clock.py',
)


def main():
    results = []
    for name in CHECKS:
        print('Checking ' + name, flush=True)
        path = HERE / name
        run = subprocess.run([sys.executable, '-B', str(path)],
                             cwd=path.parent, capture_output=True, text=True)
        results.append(dict(script=name, returncode=run.returncode,
                            stdout=run.stdout, stderr=run.stderr))
        if run.returncode:
            print(run.stdout + run.stderr, flush=True)
            raise SystemExit(run.returncode)
    print('Checking native ordered-run certificate', flush=True)
    with tempfile.TemporaryDirectory(prefix='odd-table-native-') as directory:
        binary = Path(directory) / 'native'
        subprocess.run(['c++', '-std=c++17', '-O2',
                        str(HERE / 'history/native_repair.cpp'), '-o', str(binary)],
                       check=True)
        run = subprocess.run([str(binary)], cwd=HERE / 'history',
                             capture_output=True, text=True)
        results.append(dict(script='history/native_repair.cpp',
                            returncode=run.returncode, stdout=run.stdout,
                            stderr=run.stderr))
        if run.returncode:
            raise SystemExit(run.returncode)
    (HERE / 'REPLAY.json').write_text(json.dumps(results, indent=2) + '\n')
    print('All12 new-claim engines passed.', flush=True)


if __name__ == '__main__':
    main()
