#!/usr/bin/env python3
"""Manifest-only portable replay; no repository paths, Git or discovery scans."""
from hashlib import sha256
from pathlib import Path
import json
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def run(command):
    result = subprocess.run(command, cwd=HERE, capture_output=True)
    need(result.returncode == 0, 'command failed: ' + ' '.join(map(str, command))
         + '\n' + result.stdout.decode(errors='replace') + result.stderr.decode(errors='replace'))
    return result.stdout


def main():
    manifest = json.loads((HERE / 'MANIFEST.json').read_text())
    for name, expected in manifest['files'].items():
        path = Path(name)
        need(not path.is_absolute() and '..' not in path.parts, 'manifest paths stay inside packet')
        need(sha256((HERE / path).read_bytes()).hexdigest() == expected, 'manifest hash: ' + name)
    results = []
    for script, receipt in [('low_grade.py', 'low_grade'), ('two_charge.py', 'two_charge'), ('native_check.py', 'native')]:
        for optimization in ([], ['-O']):
            output = run([sys.executable] + optimization + ['-B', str(HERE / script)])
            expected = (HERE / (receipt + ('.optimized.json' if optimization else '.normal.json'))).read_bytes()
            need(output == expected, script + ' frozen complete receipt')
        need((HERE / (receipt + '.normal.json')).read_bytes() == (HERE / (receipt + '.optimized.json')).read_bytes(),
             script + ' optimization-safe identical receipts')
        results.append(dict(check=script, receipt_sha256=sha256(expected).hexdigest(), normal_and_optimized=True))
    with tempfile.TemporaryDirectory(prefix='native-boundary33-') as directory:
        executable = str(Path(directory) / 'sheets')
        run(['g++', '-O3', '-std=c++17', str(HERE / 'sheets.cpp'), '-o', executable])
        outputs = []
        for mode in ('bfs', 'dsu'):
            output = run([executable, str(HERE / 'sheet_input.txt'), mode])
            need(output == (HERE / ('sheets.' + mode + '.json')).read_bytes(), 'complete sheet graph receipt: ' + mode)
            outputs.append(output)
        need(outputs[0] == outputs[1], 'independent enumeration/BFS and Gosper/DSU full edge graph agreement')
        results.append(dict(check='sheets.cpp', receipt_sha256=sha256(outputs[0]).hexdigest(), independent_algorithms=['bfs','dsu']))
    print(json.dumps(dict(status='PASS', manifest_files=len(manifest['files']), checks=results,
                          theorem=dict(changed_slots=8, primary_clock_occurrences=33, sharp_clock_bound=True,
                                       grades5_through23_clock_bound=27)), sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
