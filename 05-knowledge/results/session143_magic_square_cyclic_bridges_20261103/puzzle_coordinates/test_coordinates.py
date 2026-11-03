#!/usr/bin/env python3
"""Positive and mutation regression under normal Python and -O (stdlib only)."""
import json,pathlib,subprocess,sys
checker=pathlib.Path(__file__).resolve().with_name('check_coordinates.py')
results=[]
for optimized in (False,True):
    prefix=[sys.executable]+(['-O'] if optimized else [])+[str(checker)]
    for mutation in [None,'bad-move','bad-conjugator','bad-adjacency']:
        result=subprocess.run(prefix+(['--mutate',mutation] if mutation else []),capture_output=True,text=True)
        expected=1 if mutation else 0
        if result.returncode!=expected:
            raise SystemExit(f'FAIL: {mutation}, optimized={optimized}, got {result.returncode}, wanted {expected}: {result.stderr}')
        row={'optimized':optimized,'mutation':mutation,'exit_code':result.returncode}
        if mutation:
            if result.stdout or not result.stderr.startswith('FAIL: '):
                raise SystemExit('FAIL: corruption did not hit an explicit validation gate')
            row['rejection']=result.stderr.strip()
        else:
            output=json.loads(result.stdout)
            if output.get('status')!='PASS' or output.get('checks_executed')!=609:
                raise SystemExit('FAIL: expected all 609 explicit checks')
            row['checks_executed']=output['checks_executed']
        results.append(row)
print(json.dumps({'status':'PASS','tests':results,'boundary':'Specified corruptions are rejected; this is not an exhaustive proof of checker correctness.'},indent=2))
