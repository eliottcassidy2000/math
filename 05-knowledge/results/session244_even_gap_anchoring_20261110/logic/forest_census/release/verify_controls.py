#!/usr/bin/env python3
"""Reproduce strict builds, independent recounts, and hostile controls.

Uses private temporary copies; never modifies the published CSV files.
Run: python3 -I -B verify_controls.py > CONTROLS.json
"""
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
FLAGS = ['-std=c++17', '-O3', '-Wall', '-Wextra', '-Wconversion', '-Werror', '-pedantic']


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def run(command):
    result = subprocess.run(list(map(str, command)), stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    return result


def check_run(command):
    result = run(command)
    need(result.returncode == 0,
         f'command failed: {command}: {result.stderr.decode(errors="replace")}')
    return result


def main():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    work = Path(tempfile.mkdtemp(prefix='forest-census-controls-',
                dir=os.environ.get('FOREST_CONTROL_TMPDIR')))
    receipts = {'status': 'PASS', 'strict_compiler_flags': FLAGS,
                'small_run_order': 8, 'controls': []}
    for name, source in [('census', 'census.cpp'), ('audit', 'audit_all_witnesses.cpp'),
                         ('normal-form', 'verify_normal_form.cpp'),
                         ('shallow', 'shallow_census.cpp')]:
        for mode, extra in [('normal', []), ('ndebug', ['-DNDEBUG'])]:
            check_run(['g++', *FLAGS, *extra, ROOT / source, '-o', work / f'{name}-{mode}'])
    outputs = []
    for mode in ['normal', 'ndebug']:
        directory = work / mode
        directory.mkdir()
        result = check_run([work / f'census-{mode}', '8', directory])
        (directory / 'SUMMARY.jsonl').write_bytes(result.stdout)
        outputs.append(result.stdout)
    need(outputs[0] == outputs[1], 'normal/-DNDEBUG summaries differ')
    check_run(['patch', '--reverse', '--silent', '--output', work / 'census-running.cpp',
               ROOT / 'census.cpp', ROOT / 'ARGUMENT_VALIDATION.patch'])
    provenance = json.loads((ROOT / 'RUN_PROVENANCE.json').read_text())
    need(digest((work / 'census-running.cpp').read_bytes()) == provenance['running_source_sha256'],
         'reconstructed running-source hash mismatch')
    check_run(['g++', *FLAGS, work / 'census-running.cpp', '-o', work / 'census-running'])
    original = work / 'original-wrapper'
    original.mkdir()
    original_result = check_run([work / 'census-running', '8', original])
    need(original_result.stdout == outputs[0], 'running/published wrapper summary mismatch')
    for n in range(1, 9):
        filename = f'n{n}_witnesses.csv'
        need((original / filename).read_bytes() == (work / 'normal' / filename).read_bytes(),
             f'running/published wrapper CSV mismatch at order {n}')
    receipts['controls'].append({'name': 'reconstructed_running_source_and_wrapper_only_replay',
                                'result': 'PASS',
                                'running_source_sha256': provenance['running_source_sha256'],
                                'summary_sha256': digest(original_result.stdout)})
    csv_checks = []
    for n in range(1, 9):
        name = f'n{n}_witnesses.csv'
        a, b = [(work / mode / name).read_bytes() for mode in ['normal', 'ndebug']]
        need(a == b, f'normal/-DNDEBUG {name} differs')
        csv_checks.append({'file': name, 'sha256': digest(a)})
    receipts['controls'].append({'name': 'normal_ndebug_byte_equality', 'result': 'PASS',
                                'summary_sha256': digest(outputs[0]), 'csv_files': csv_checks})
    audit_outputs = []
    for mode in ['normal', 'ndebug']:
        result = check_run([work / f'audit-{mode}', '8', work / mode])
        audit_outputs.append(result.stdout)
    need(audit_outputs[0] == audit_outputs[1], 'normal/-DNDEBUG audits differ')
    receipts['controls'].append({'name': 'small_run_all_codes_recount', 'result': 'PASS',
                                'stdout_sha256': digest(audit_outputs[0])})
    python_outputs = []
    for flags in [['-I', '-B'], ['-I', '-B', '-O']]:
        result = check_run([sys.executable, *flags, ROOT / 'audit.py', work / 'normal'])
        python_outputs.append(result.stdout)
    need(python_outputs[0] == python_outputs[1], 'Python normal/-O audit outputs differ')
    receipts['controls'].append({'name': 'python_normal_optimized_byte_equality', 'result': 'PASS',
                                'stdout_sha256': digest(python_outputs[0])})
    for name in ['census', 'audit']:
        for bad in ['13', '4294967298', '2junk', '-2', '0']:
            result = run([work / f'{name}-normal', bad, work / 'normal'])
            need(result.returncode != 0, f'{name} accepts bad order {bad}')
            receipts['controls'].append({'name': f'{name}_rejects_bad_order', 'input': bad,
                                        'result': 'PASS', 'returncode': result.returncode,
                                        'stderr': result.stderr.decode().strip()})
    normal_form_outputs = [check_run([work / f'normal-form-{mode}']).stdout
                           for mode in ['normal', 'ndebug']]
    need(normal_form_outputs[0] == normal_form_outputs[1], 'normal-form -DNDEBUG mismatch')
    need(normal_form_outputs[0] == (ROOT / 'NORMAL_FORM_CHECK.json').read_bytes(),
         'normal-form stored certificate mismatch')
    receipts['controls'].append({'name': 'normal_form_independent_literal_replay_ndebug_equality',
                                'result': 'PASS', 'stdout_sha256': digest(normal_form_outputs[0])})
    for normalized in ['0', '1']:
        shallow_outputs = [check_run([work / f'shallow-{mode}', '8', work / 'normal', normalized]).stdout
                           for mode in ['normal', 'ndebug']]
        need(shallow_outputs[0] == shallow_outputs[1], 'shallow normal/-DNDEBUG mismatch')
        receipts['controls'].append({'name': 'restricted_shallow_normal_ndebug_byte_equality',
                                    'normalization': normalized, 'result': 'PASS',
                                    'report': json.loads(shallow_outputs[0]),
                                    'stdout_sha256': digest(shallow_outputs[0])})
    shallow_reject = run([work / 'shallow-normal', '13', work / 'normal', '1'])
    need(shallow_reject.returncode != 0, 'shallow tool accepts order13')
    receipts['controls'].append({'name': 'shallow_rejects_order13', 'result': 'PASS',
                                'returncode': shallow_reject.returncode})
    exception_outputs = [check_run([sys.executable, *flags, ROOT / 'verify_depth_exception.py']).stdout
                         for flags in [['-I', '-B'], ['-I', '-B', '-O']]]
    need(exception_outputs[0] == exception_outputs[1], 'depth exception normal/-O mismatch')
    need(exception_outputs[0] == (ROOT / 'DEPTH_EXCEPTION_CHECK.json').read_bytes(),
         'stored depth-exception certificate mismatch')
    receipts['controls'].append({'name': 'depth_exception_and_order20_join_normal_optimized_equality',
                                'result': 'PASS', 'stdout_sha256': digest(exception_outputs[0])})
    damaged = work / 'corrupted-private-copy'
    shutil.copytree(work / 'normal', damaged)
    path = damaged / 'n3_witnesses.csv'
    lines = path.read_text().splitlines()
    need(lines[-1].startswith('3,'), 'expected positive control missing')
    lines[-1] = '5,' + lines[-1].split(',', 1)[1]
    path.write_text('\n'.join(lines) + '\n')
    for name, command in [
        ('cpp', [work / 'audit-normal', '8', damaged]),
        ('python', [sys.executable, '-I', '-B', ROOT / 'audit.py', damaged]),
        ('python_optimized', [sys.executable, '-I', '-B', '-O', ROOT / 'audit.py', damaged])]:
        result = run(command)
        need(result.returncode != 0, f'{name} accepts corrupted CSV')
        receipts['controls'].append({'name': f'{name}_rejects_corrupted_private_csv',
                                    'result': 'PASS', 'mutation': 'n3: H=3 changed to H=5; witness unchanged',
                                    'returncode': result.returncode,
                                    'error': result.stderr.decode().strip().splitlines()[-1]})
    receipts['artifact_source_sha256'] = {
        name: digest((ROOT / name).read_bytes())
        for name in ['census.cpp', 'audit_all_witnesses.cpp', 'audit.py', 'verify_controls.py',
                     'verify_normal_form.cpp', 'shallow_census.cpp', 'verify_depth_exception.py']}
    print(json.dumps(receipts, indent=2))


if __name__ == '__main__':
    main()
