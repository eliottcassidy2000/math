"""Replay the nine new-session exact engines; no historical blanket audit."""
from pathlib import Path
import hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parent
OPT='--optimized' in sys.argv
exe=[sys.executable]+(['-O'] if OPT else [])+['-B']
engines=[
 ('endpoint_cut','endpoints/check.py','endpoints/results_optimized.json' if OPT else 'endpoints/results.json'),
 ('nine_state','seed_decision/check.py',None),
 ('independent_seed','seed_decision/audit_independent.py','seed_decision/audit_receipt_optimized.json' if OPT else 'seed_decision/audit_receipt.json'),
 ('strong_core','seed_decision/core_check.py','seed_decision/core_receipt_optimized.json' if OPT else 'seed_decision/core_receipt.json'),
 ('physical_grammar','grammar/compile.py',None),
 ('positive_seed_probe','grammar/seed_probe_check.py',None),
 ('zeckendorf_guard','grammar/history_bridge_check.py',None),
 ('arithmetic','arithmetic/check.py',None),
 ('literal_constructor','construct_seed.py',None)]
rows=[]
for name,path,output in engines:
 p=ROOT/path
 r=subprocess.run(exe+[p.name]+(['--self-test'] if name=='literal_constructor' else []),cwd=p.parent,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},capture_output=True,text=True)
 if r.returncode:raise RuntimeError(f'{name} failed\n{r.stdout}\n{r.stderr}')
 if output:(ROOT/output).write_text(r.stdout)
 rows.append({'engine':name,'path':path,'status':'PASS','stdout_sha256':hashlib.sha256(r.stdout.encode()).hexdigest()})
 print('PASS',name,flush=True)
for a,b in [('endpoints/results.json','endpoints/results_optimized.json'),('seed_decision/audit_receipt.json','seed_decision/audit_receipt_optimized.json'),('seed_decision/core_receipt.json','seed_decision/core_receipt_optimized.json')]:
 if (ROOT/a).read_bytes()!=(ROOT/b).read_bytes():raise RuntimeError('Normal/optimized mismatch: '+a)
out={'status':'PASS','optimized':OPT,'scope':'Nine new-session engines; inherited literature and corpus are inputs, not blanket reruns','engines':rows,'paired_receipts':'byte-identical'}
path=ROOT/('REPLAY_OPTIMIZED.json' if OPT else 'REPLAY.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print('PASS all nine engines; paired normal/optimized receipts byte-identical')
