#!/usr/bin/env python3
import json,pathlib,subprocess,sys
checker=pathlib.Path(__file__).resolve().with_name('check_corollary.py');tests=[]
for optimized in (False,True):
    for mutate in (False,True):
        cmd=[sys.executable]+(['-O'] if optimized else [])+[str(checker)]+(['--mutate','bad-centralizer'] if mutate else [])
        r=subprocess.run(cmd,text=True,capture_output=True);expected=int(mutate)
        if r.returncode!=expected:raise SystemExit(f'FAIL: expected exit{expected}, got{r.returncode}: {r.stderr}')
        row={'optimized':optimized,'mutation':mutate,'exit_code':r.returncode}
        if mutate:
            if r.stdout or not r.stderr.startswith('FAIL: '):raise SystemExit('FAIL: mutation did not hit an explicit gate')
            row['rejection']=r.stderr.strip()
        else:
            obj=json.loads(r.stdout)
            if obj.get('status')!='PASS' or obj.get('checks_executed')!=77:raise SystemExit('FAIL: checker gates incomplete')
            row['checks_executed']=obj['checks_executed']
        tests.append(row)
print(json.dumps({'status':'PASS','tests':tests},indent=2))
