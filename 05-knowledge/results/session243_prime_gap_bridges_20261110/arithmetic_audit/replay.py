#!/usr/bin/env python3
"""Replay all arithmetic verifiers and freeze the exact audited sources."""
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = PACKET.parents[2]


def demand(ok, message):
    if not ok:
        raise RuntimeError(message)


def run(job):
    lane, optimized, failure = job
    relative = f'{lane}/verify.py'
    argv = [sys.executable, '-I', '-B']
    if optimized:
        argv.append('-O')
    argv.append(str(PACKET/relative))
    if failure is not None:
        argv.append('--fail')
        if failure:
            argv.append(failure)
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    expected = 0 if failure is None else 1
    demand(result.returncode == expected, f'{job}: exit {result.returncode}, wanted {expected}: {result.stderr}')
    if failure is None:
        payload = json.loads(result.stdout)
        demand(payload['status'] == 'PASS', f'{job}: missing PASS')
        demand(not result.stderr, f'{job}: unexpected stderr')
    else:
        payload = None
        demand(not result.stdout, f'{job}: failure printed success')
        demand('deliberate' in result.stderr or 'intentional' in result.stderr,
               f'{job}: unrelated failure instead of intended gate')
    return {'lane':lane,'optimized':optimized,'failure':failure,'exit_code':result.returncode,
            'stdout_sha256':sha256(result.stdout.encode()).hexdigest(),
            'failure_message':result.stderr.splitlines()[-1] if result.stderr else None}, result.stdout, payload


def main():
    failures = {'brown':[''], 'provenance':['skip-one','wheel-prime','brown-packet'],
                'arithmetic_audit':['skip-one','fixed-offset','brown-quotient','sign-exponent']}
    frozen = {str(path.relative_to(PACKET)): sha256(path.read_bytes()).hexdigest()
              for lane in failures for path in (PACKET/lane).iterdir()
              if path.name in ['PROOF.md', 'verify.py', 'source_pins.json']}
    jobs = [(lane, opt, failure) for lane, gates in failures.items()
            for failure in [None]+gates for opt in [False,True]]
    with ThreadPoolExecutor(max_workers=3) as pool:
        runs = list(pool.map(run, jobs))
    for relative, digest in frozen.items():
        demand(sha256((PACKET/relative).read_bytes()).hexdigest() == digest,
               'source changed during replay: '+relative)
    positive = {}
    for lane in failures:
        normal = next(r for r in runs if r[0]['lane']==lane and not r[0]['optimized'] and r[0]['failure'] is None)
        optimized = next(r for r in runs if r[0]['lane']==lane and r[0]['optimized'] and r[0]['failure'] is None)
        demand(normal[1] == optimized[1], f'{lane}: normal/-O differ')
        positive[lane] = normal[2]
        if lane == 'arithmetic_audit':
            (HERE/'CHECK.json').write_text(normal[1])
            (HERE/'CHECK.optimized.json').write_text(optimized[1])
        else:
            demand((PACKET/lane/'CHECK.json').read_text() == normal[1], f'{lane}/CHECK.json: stale output')
            if (PACKET/lane/'CHECK.optimized.json').exists():
                demand((PACKET/lane/'CHECK.optimized.json').read_text() == normal[1],
                       f'{lane}/CHECK.optimized.json: stale output')
    independent = positive['arithmetic_audit']
    demand(positive['brown']['selection_count'] == independent['brown']['partition_count'], 'Brown count mismatch')
    demand(positive['brown']['brown_roots'] == [r[:2] for r in independent['brown']['roots']], 'Brown roots mismatch')
    p = positive['provenance']
    demand(p['Brown_solutions'] == independent['brown']['roots'], 'provenance roots mismatch')
    expected_subsets = 2*sum(r['partitions'] for r in independent['brown']['rows'] if r['n'] <= 40)
    demand(p['Brown_subsets_checked'] == expected_subsets, 'provenance subset count mismatch')
    demand(p['prime_square_unions_checked'] == independent['prime_square_brocard']['pair_count'], 'bracket count mismatch')
    for earlier, audit in zip(p['primorial_rows'], independent['prime_gaps']['rows']):
        for key in ['n','P','F','L','nearest_pair']:
            demand(earlier[key] == audit[key], f'primorial {key} mismatch')
        demand(earlier['crossing_gap'] == audit['gap'], 'crossing gap mismatch')
    pins = json.loads((PACKET/'provenance/source_pins.json').read_text())
    checked_pins = []
    for pin in pins['sources']:
        obj = f"{pins['ref']}:{pin['path']}"
        blob = subprocess.run(['git','-C',str(ROOT),'show',obj],capture_output=True,check=True).stdout
        oid = subprocess.run(['git','-C',str(ROOT),'rev-parse',obj],capture_output=True,text=True,check=True).stdout.strip()
        demand(sha256(blob).hexdigest() == pin['sha256'] and oid == pin['git_blob'], 'source pin mismatch '+pin['path'])
        checked_pins.append({'path':pin['path'],'sha256':pin['sha256'],'git_blob':oid})
    sources = []
    for lane in failures:
        for filename in ['PROOF.md','verify.py','source_pins.json']:
            path = PACKET/lane/filename
            if path.exists():
                sources.append({'path':str(path.relative_to(PACKET)),'sha256':sha256(path.read_bytes()).hexdigest()})
    report = {'status':'PASS','positive_replays':6,'negative_replays':16,
              'normal_optimized_outputs_identical':True,'producer_output_comparison':'PASS',
              'source_pins_ref':pins['ref'],'source_pins_checked':checked_pins,
              'audited_source_hashes':sources,'runs':[r[0] for r in runs]}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
