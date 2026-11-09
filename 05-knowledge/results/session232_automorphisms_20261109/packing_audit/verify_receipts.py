#!/usr/bin/env python3
import hashlib,json,subprocess,sys
from pathlib import Path
D=Path(__file__).resolve().parent
P=D.parent/'packing_mass'
if not P.exists(): P=D.parent/'session232-packing-mass'
def must(x,msg):
    if not x: raise RuntimeError(msg)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
receipts=[]
for name,path in [('producer',P/'check.py'),('independent',D/'independent_check.py')]:
    outputs=[]; artifacts=[]
    for optimized in (False,True):
        cmd=[sys.executable]+(['-O'] if optimized else [])+[str(path)]
        r=subprocess.run(cmd,capture_output=True,text=True,check=False)
        must(r.returncode==0,f'{name} failed: {r.stderr}')
        must(not r.stderr,f'{name} stderr: {r.stderr}')
        output=D/(name.upper()+('.optimized' if optimized else '')+'.json')
        output.write_text(r.stdout)
        outputs.append(r.stdout)
        receipts.append({'command':cmd,'exit_code':r.returncode,'stdout_sha256':sha(output)})
        if name=='independent':artifacts.append({f:sha(D/f) for f in ['WITNESSES.json','RESTRICTIONS.json']})
    must(outputs[0]==outputs[1],f'{name} normal/-O disagreement')
    if artifacts:must(artifacts[0]==artifacts[1],'independent normal/-O witness disagreement')
p=json.loads((D/'PRODUCER.json').read_text()); i=json.loads((D/'INDEPENDENT.json').read_text())
for left,right in [(p['C3'],i['C3'])]+list(zip(p['tower'],i['tower'])):
    for k in ['n','H','path_orbits','coarse_packing_orbits','stabilizer_histogram']:
        must(left[k]==right[k],f'producer/independent disagreement on {k}')
    must(left['aut_or_subgroup']==right['group_order'],'group size disagreement')
must((D/'PRODUCER.json').read_bytes()==(P/'RESULTS.json').read_bytes(),'stored producer output drift')
report={'status':'PASS','runs':receipts,'producer_independent_scalar_and_histogram_agreement':True,
        'normal_optimized_stdout_and_witness_agreement':True,
        'inputs':{str(p):sha(p) for p in [P/'PROOF.md',P/'check.py',D/'independent_check.py']},
        'artifact_sha256':{f:sha(D/f) for f in ['WITNESSES.json','RESTRICTIONS.json']}}
(D/'RECEIPTS.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(json.dumps(report,sort_keys=True,indent=2))
