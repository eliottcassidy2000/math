#!/usr/bin/env python3
"""Frozen relocation and dependency hostiles; no new mathematical experiments."""
from pathlib import Path
from hashlib import sha256
import json, os, shutil, subprocess, sys, tempfile

HERE = Path(__file__).resolve().parent
BASE = HERE.parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def run(packet, optimized, script='session65-audit-synthesis/portable_check.py'):
    command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [str(packet/script)]
    return subprocess.run(command, cwd=packet, env={'PATH': '/nonexistent'},
                          capture_output=True, check=False, timeout=120)

inputs = json.loads((HERE/'PORTABLE_INPUTS.json').read_text())
with tempfile.TemporaryDirectory(prefix='session65-synthesis-portable-') as temporary:
    packet = Path(temporary)
    local_paths = [r['source'] for r in inputs['files']]
    extra = ['session65-audit-synthesis/PORTABLE_INPUTS.json',
             'session65-audit-synthesis/portable_check.py']
    for relative in local_paths + extra:
        source = BASE/relative
        target = packet/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    need(not (packet/'session65-repo').exists(), 'ambient repository copied')
    need(not (packet/'session65-fifteens/sloane1973p98.jpg').exists(), 'excluded image copied')
    need(not (packet/'session65-weakcollatz/tinycryptids').exists(), 'unshipped external source copied')
    need(not any(p.name == '.git' for p in packet.rglob('*')), 'Git metadata copied')
    positive = [run(packet, opt) for opt in (False, True)]
    need(all(p.returncode == 0 and not p.stderr for p in positive), 'frozen replay failed')
    need(positive[0].stdout == positive[1].stdout, 'frozen normal/-O disagreement')
    need(positive[0].stdout == (HERE/'portable-normal.json').read_bytes(), 'frozen result differs')
    result = json.loads(positive[0].stdout)
    historic = json.loads((HERE/'normal.json').read_text())
    math_keys = ('carriers', 'results', 'scope')
    need(all(result[k] == historic[k] for k in math_keys), 'historical mathematics changed')
    old = (HERE/'independent_check.py').read_text()
    new = (HERE/'portable_check.py').read_text()
    start = '# Generate all matchings by filtering every three-element subset of the full\n'
    end = 'print(json.dumps('
    old_body = old[old.index(start):old.index(end)]
    new_body = new[new.index(start):new.index(end)]
    need(old_body == new_body, 'mathematical body changed')
    hostile_rows = []
    for relative in local_paths:
        target = packet/relative
        original = target.read_bytes()
        modes_checked = []
        for mutation in ('missing_file', 'wrong_bytes_same_length'):
            try:
                if mutation == 'missing_file':
                    target.unlink()
                elif original:
                    target.write_bytes(bytes((original[0] ^ 1,)) + original[1:])
                else:
                    target.write_bytes(b'hostile')
                exits = []
                for opt in (False, True):
                    response = run(packet, opt)
                    need(response.returncode != 0 and not response.stdout,
                         'mutated dependency accepted: ' + relative)
                    exits.append(response.returncode)
                modes_checked.append({'mutation': mutation, 'normal_exit': exits[0],
                                      'optimized_exit': exits[1]})
            finally:
                target.write_bytes(original)
        hostile_rows.append({'path': relative, 'rejections': modes_checked})
    exclusion_attacks = []
    for relative, field in (('session65-audit-synthesis/EXTERNAL_NONRUNTIME_INPUTS.json', 'inputs'),
                            ('session65-audit-synthesis/PORTABLE_INPUTS.json', 'external_nonruntime_inputs')):
        target = packet/relative
        original = target.read_bytes()
        altered = json.loads(original)
        altered[field].append({'path': 'session65-weakcollatz/bridge_results_recheck.json',
                               'classification': 'EXTERNAL_NONRUNTIME',
                               'source_url': 'https://invalid.example/undeclared-exclusion',
                               'sha256': '65f0489ce2c16784dae5866f0c36be27e960dbb9174fcd6d76cebd316e5a2bfc'})
        missing = packet/'session65-weakcollatz/bridge_results_recheck.json'
        missing_original = missing.read_bytes()
        try:
            target.write_text(json.dumps(altered))
            missing.unlink()
            exits = []
            for opt in (False, True):
                response = run(packet, opt)
                need(response.returncode != 0 and not response.stdout,
                     'undeclared exclusion accepted: ' + relative)
                exits.append(response.returncode)
            exclusion_attacks.append({'altered_declaration': relative,
                                      'attempted_exclusion': 'session65-weakcollatz/bridge_results_recheck.json',
                                      'normal_exit': exits[0], 'optimized_exit': exits[1]})
        finally:
            target.write_bytes(original)
            missing.write_bytes(missing_original)
    # The extra receipt is authored deterministic output, not omitted provenance.
    producer = [run(packet, opt, 'session65-weakcollatz/verify_bridge.py') for opt in (False, True)]
    expected = (packet/'session65-weakcollatz/bridge_results_recheck.json').read_bytes()
    need(all(p.returncode == 0 and not p.stderr and p.stdout == expected for p in producer),
         'authored bridge recheck cannot be reproduced')
    need(run(packet, False).stdout == positive[0].stdout, 'restored package failed')
    print(json.dumps({'status': 'PASS_PORTABLE_FROZEN_CLOSURE_AND_HOSTILES',
        'frozen_environment': {'original_repository_present': False, 'git_on_PATH': False,
                              'ambient_PYTHONPATH_used': False, 'isolated_python': True,
                              'excluded_image_present': False, 'unshipped_tinycryptids_present': False},
        'local_dependencies_checked': len(local_paths),
        'inherited_local_seal_rows': result['inherited_seal_rows_checked'],
        'external_nonruntime_observations': result['external_nonruntime_observations'],
        'normal_optimized_byte_identical': True,
        'historical_mathematical_keys_identical': list(math_keys),
        'mathematical_body_byte_identical': True,
        'mathematical_body_sha256': sha256(old_body.encode()).hexdigest(),
        'missing_and_wrong_byte_hostile_runs': 4*len(hostile_rows),
        'dependency_hostiles': hostile_rows,
        'undeclared_exclusion_hostile_runs': 2*len(exclusion_attacks),
        'undeclared_exclusion_hostiles': exclusion_attacks,
        'authored_bridge_recheck_reproduced_both_modes': True,
        'authored_bridge_recheck_sha256': sha256(expected).hexdigest()},
        sort_keys=True, indent=2))
