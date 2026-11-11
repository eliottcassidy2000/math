#!/usr/bin/env python3
"""Reproduce only Session269's new operations and controls."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
CASES=[('jacobian','verify.py'),('hadwiger','replay.py'),('hadwiger','equality.py'),('chanpak','reproduce.py'),('bridges','check.py')]
report={'status':'PASS','scope':'new session mechanisms only; no inherited-claim audit','cases':[]}
for folder,script in CASES:
    outputs=[]
    for opt in [False,True]:
        args=[sys.executable]+(['-O'] if opt else [])+[script]
        p=subprocess.run(args,cwd=ROOT/folder,capture_output=True,text=True,check=True)
        outputs.append(p.stdout)
    if outputs[0]!=outputs[1]:raise RuntimeError('optimization mismatch: '+folder+'/'+script)
    report['cases'].append({'script':folder+'/'+script,'normal':'PASS','optimized':'PASS','identical_stdout':True,'stdout_sha256':hashlib.sha256(outputs[0].encode()).hexdigest()})
    print('PASS',folder+'/'+script,flush=True)
report['research_status']={'JC2':'OPEN','Hadwiger_Nelson':'OPEN','Chan_Pak':'OPEN'}
(ROOT/'FINAL_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS all five new-claim programs in both modes')
