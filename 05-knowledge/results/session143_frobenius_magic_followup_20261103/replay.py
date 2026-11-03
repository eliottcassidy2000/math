#!/usr/bin/env python3
"""Replay both frozen implementations, hostile controls and agreement offline."""
from hashlib import sha256
from pathlib import Path
import json
import os
import subprocess
import sys
from verify_manifest import verify

ROOT = Path(__file__).resolve().parent

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def main():
    verify()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    tasks = [
        ('producer/check_bridge.py', 'producer/result.json', 'producer/result.optimized.json'),
        ('producer/test_controls.py', 'producer/hostile_results.json', 'producer/hostile_results.optimized.json'),
        ('independent/check_independent.py', 'independent/result.json', 'independent/result.optimized.json'),
        ('independent/test_negative_controls.py', 'independent/negative_controls.json', 'independent/negative_controls.optimized.json'),
        ('independent/audit_producer.py', 'independent/producer_review.json', 'independent/producer_review.json'),
    ]
    executions = []
    for program, normal, optimized in tasks:
        for mode, expected in (('normal', normal), ('optimized', optimized)):
            command = [sys.executable] + (['-O'] if mode == 'optimized' else [])
            run = subprocess.run(command + [str(ROOT / program)], cwd=ROOT,
                                 env=env, capture_output=True, timeout=180)
            require(run.returncode == 0, 'failed ' + program + ': ' + run.stderr.decode())
            require(run.stdout == (ROOT / expected).read_bytes(), 'output mismatch: ' + expected)
            executions.append({'program': program, 'mode': mode,
                               'expected_output': expected, 'status': 'PASS',
                               'output_sha256': sha256(run.stdout).hexdigest()})
    for program in ('producer/verify_manifest.py', 'independent/verify_manifest.py', 'verify_source_pins.py'):
        for mode in ('normal', 'optimized'):
            command = [sys.executable] + (['-O'] if mode == 'optimized' else [])
            run = subprocess.run(command + [str(ROOT / program)], cwd=ROOT,
                                 env=env, capture_output=True, timeout=60)
            require(run.returncode == 0, 'integrity failed ' + program + ': ' + run.stderr.decode())
            result = json.loads(run.stdout)
            require(result['status'] == 'PASS', 'integrity status: ' + program)
            executions.append({'program': program, 'mode': mode, 'status': 'PASS',
                               'output_sha256': sha256(run.stdout).hexdigest()})
    a = json.loads((ROOT / 'producer/result.json').read_text())
    b = json.loads((ROOT / 'independent/result.json').read_text())
    rows_a = sorted((r['a'], r['b'], r['t'], r['n5'], r['n25'], r['n125']) for r in a['models'])
    rows_b = sorted((r['a'], r['b'], r['t'], r['N1'], r['N2'], r['N3']) for r in b['models'])
    require(rows_a == rows_b and len(rows_a) == 20, 'independent model-count disagreement')
    ca = sorted((tuple(r['representative']), r['models'], r['t']) for r in a['isomorphism_classes'])
    cb = sorted((tuple(r['representative']), r['models'], r['trace']) for r in b['isomorphism_classes'])
    require(ca == cb and len(ca) == 12, 'independent class disagreement')
    require(a['magic_arrangement_count'] == b['magic_arrangement_counts']['N3'] == 8, 'independent complete magic enumeration')
    require(a['class_trace_fibers'] == b['trace_fibers'], 'independent trace fibers')
    require(a['magic_square'] == b['magic_grid'], 'independent magic witness')
    require(a['line_sums'] == b['magic_line_sums'] == [378] * 8, 'independent magic lines')
    require(a['twist_orbits'] == b['twist_orbits'] == 6, 'independent twist quotient')
    tc = b['tournament_control']
    require(len(a['tournament_classes']) == tc['classes'] == 12, 'independent tournament classes')
    require(a['tournament_spectrum_fibers'] == tc['spectrum_fiber_partition'], 'independent spectrum fibers')
    require(a['tournament_converse_orbits'] == tc['converse_orbits'] == 10, 'independent converse quotient')
    require(a['tournament_selfdual_classes'] == tc['self_converse'] == 8, 'independent fixed classes')
    require(a['f25_control']['indexed_values'] == b['F25_values'], 'independent quadratic hostile')
    require(a['square_entries'] == [144], 'squarehood boundary')
    verify()
    print(json.dumps({'status': 'PASS', 'integrity': 'PASS', 'executions': executions,
                      'cross_check': {'models': 20, 'classes': 12, 'trace_states': 9,
                                      'magic_line_sum': 378, 'square_entries': [144],
                                      'twist_orbits': 6, 'tournament_converse_orbits': 10},
                      'scope': 'Finite F5 short-model census and degrees 1, 2, 3 only; no MSOS solution.'},
                     indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
