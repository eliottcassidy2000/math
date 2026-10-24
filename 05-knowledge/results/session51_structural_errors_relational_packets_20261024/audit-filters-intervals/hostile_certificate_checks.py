#!/usr/bin/env python3
"""Mutate literal certificate files and execute the copied verifier in both modes."""
import copy
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'copies/prime-offsets'
original = json.loads((SOURCE / 'certificates.json').read_text())
mutants = {}


def mutate(name, fn):
    data = copy.deepcopy(original)
    fn(data)
    mutants[name] = data


mutate('missing final cutoff index', lambda d: d['rows'].pop())
mutate('incorrect primorial product', lambda d: d['rows'][8].__setitem__('P', d['rows'][8]['P'] + 2))
mutate('omitted interior minimality witness', lambda d: d['rows'][5]['plus']['witnesses'].pop(4))
mutate('proper divisor is one', lambda d: d['rows'][0]['plus']['witnesses'][0].__setitem__('divisor', 1))
mutate('prime falsely declared composite by MR', lambda d: d['rows'][0]['plus']['witnesses'].__setitem__(-1, {'mr_base': 2}))
mutate('composite falsely declared prime', lambda d: d['rows'][5]['plus']['witnesses'].__setitem__(0, {'prime': True}))
mutate('missing recursive prime factor certificate', lambda d: d['primality_certificates'].pop('2'))
mutate('Lucas base equals modulus', lambda d: d['primality_certificates']['30047'].__setitem__('base', 30047))
mutate('wrong complete factorization', lambda d: d['primality_certificates']['30047']['factors'][0].__setitem__(1, 100))
mutate('Fermat pseudoprime defeats Fermat-only check', lambda d: d['primality_certificates'].__setitem__('341', {'base': 2, 'factors': [[2, 2], [5, 1], [17, 1]]}))
mutate('Carmichael composite defeats Fermat-only check', lambda d: d['primality_certificates'].__setitem__('561', {'base': 2, 'factors': [[2, 4], [5, 1], [7, 1]]}))
mutate('undefined L1 fabricated', lambda d: d['rows'][0].__setitem__('minus', {'offset': 2, 'witnesses': [{'nonpositive_or_one': True}]}))


def false_mr(d):
    for row in d['rows']:
        for side, sign in [('plus', 1), ('minus', -1)]:
            if row[side] is None:
                continue
            for offset, w in enumerate(row[side]['witnesses'], 2):
                if 'mr_base' in w:
                    w['mr_base'] = row['P'] + sign * offset - 1
                    return
    raise RuntimeError('no MR witness found')


mutate('passing MR base forged as compositeness witness', false_mr)
rejections = []
with tempfile.TemporaryDirectory(prefix='session51-independent-prime-mutants-') as temp:
    tmp = Path(temp)
    for index, (name, data) in enumerate(mutants.items()):
        packet = tmp / str(index)
        packet.mkdir()
        shutil.copy2(SOURCE / 'check.py', packet / 'check.py')
        (packet / 'certificates.json').write_text(json.dumps(data))
        for optimized in (False, True):
            command = [sys.executable] + (['-O'] if optimized else []) + [str(packet / 'check.py')]
            result = subprocess.run(command, cwd=tmp, capture_output=True, text=True)
            if result.returncode == 0:
                raise RuntimeError('false certificate accepted: ' + name)
            if 'ValueError:' not in result.stderr:
                raise RuntimeError('did not fail by explicit mathematical require: ' + result.stderr)
            rejections.append({'mutation': name, 'optimized': optimized, 'exit_code': result.returncode, 'reason': result.stderr.strip().splitlines()[-1]})
print(json.dumps({'status': 'PASS', 'distinct_mutations': len(mutants), 'literal_file_subprocess_rejections': rejections}, indent=2, sort_keys=True))
