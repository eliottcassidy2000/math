#!/usr/bin/env python3
"""Run positive and deliberately corrupted checks, with and without Python -O.
Usage: python3 test_bridge.py [--sources-dir PATH]
Uses explicit exceptions, not optimization-sensitive assert statements.
"""
import argparse, json, pathlib, subprocess, sys
parser=argparse.ArgumentParser()
parser.add_argument('--sources-dir',type=pathlib.Path)
args=parser.parse_args()
checker=pathlib.Path(__file__).resolve().with_name('check_bridge.py')
results=[]
for optimized in (False,True):
    prefix=[sys.executable]+(['-O'] if optimized else [])+[str(checker)]
    tests=[('portable-positive',[],0),('bad-eigenvalue',['--mutate','bad-eigenvalue'],1),('bad-character',['--mutate','bad-character'],1)]
    if args.sources_dir:
        sourceargs=['--sources-dir',str(args.sources_dir)]
        tests.extend([('source-positive',sourceargs,0),('bad-transcription',sourceargs+['--mutate','bad-transcription'],1)])
    for name,extra,expected in tests:
        result=subprocess.run(prefix+extra,capture_output=True,text=True)
        if result.returncode!=expected:
            raise SystemExit(f'FAIL: {name}, optimization={optimized}, expected {expected}, got {result.returncode}: {result.stderr}')
        if expected==0:
            output=json.loads(result.stdout)
            if output.get('status')!='PASS' or output.get('checks_executed',0)<35:
                raise SystemExit(f'FAIL: missing successful explicit checks: {name}')
            evidence={'checks_executed':output['checks_executed']}
        else:
            if result.stdout or not result.stderr.startswith('FAIL: '):
                raise SystemExit(f'FAIL: mutation did not trigger a clean mathematical check failure: {name}')
            evidence={'rejection':result.stderr.strip()}
        results.append({'test':name,'optimized':optimized,'exit_code':result.returncode,**evidence})
print(json.dumps({'status':'PASS','tests':results,'boundary':'Negative controls certify rejection of named corruptions, not all possible checker faults.'},indent=2))
