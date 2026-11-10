#!/usr/bin/env python3
"""Re-run the complete Session238 producer/audit suite with live failure gates."""
from pathlib import Path
import hashlib,json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[3]
P=Path(__file__).resolve().parent
rel=P.relative_to(ROOT)
programs=[rel/'cyclic_cost/check.py',rel/'cyclic_cost_audit/check.py',Path('04-computation/experiments/session238_graph_groups/verify.py'),rel/'graph_groups_audit/check.py',rel/'planar_actions/verify.py',rel/'planar_actions_audit/verify.py']
rows=[]
for path in programs:
    row={'program':str(path),'runs':[]};outputs=[]
    for opt in (False,True):
        for fail in (False,True):
            cmd=[sys.executable]+(['-O'] if opt else [])+[str(path)]+(['--inject-failure'] if fail else [])
            start=time.monotonic();r=subprocess.run(cmd,cwd=ROOT,capture_output=True)
            if (r.returncode==0)==fail:
                raise RuntimeError(f'Unexpected exit for {cmd}: {r.returncode}\n{r.stderr.decode()}')
            if not fail:outputs.append(r.stdout)
            row['runs'].append({'optimized':opt,'injected_failure':fail,'returncode':r.returncode,'seconds':round(time.monotonic()-start,3),'stdout_sha256':hashlib.sha256(r.stdout).hexdigest(),'stderr_sha256':hashlib.sha256(r.stderr).hexdigest()})
    if outputs[0]!=outputs[1]:raise RuntimeError(f'Optimized mismatch: {path}')
    row['normal_optimized_match']=True;rows.append(row)
    print(f'PASS {path}',flush=True)
result={'status':'PASS','programs':len(programs),'executions':4*len(programs),'scope':'Six producer/audit programs; inherited all-order exclusion and vertex-transitive certificates not rerun','rows':rows}
(P/'REPLAY.json').write_text(json.dumps(result,indent=2)+'\n')
print('Session238 combined replay PASS',flush=True)
