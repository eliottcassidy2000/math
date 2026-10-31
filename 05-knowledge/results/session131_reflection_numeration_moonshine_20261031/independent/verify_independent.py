#!/usr/bin/env python3
from pathlib import Path
import hashlib
import json
import subprocess
import sys

HERE=Path(__file__).resolve().parent
CASES=[('check_module.py','module'),('check_moonshine_polynomials.py','moonshine'),
       ('check_numeration_certificates.py','numeration'),('check_tournament_carries.py','carries'),
       ('check_square_sum.py','square_sum')]

def require(ok,msg):
    if not ok:raise RuntimeError(msg)

rows=[]
for script,stem in CASES:
    reference=(HERE/(stem+'.normal.json')).read_bytes()
    for optimized in (False,True):
        command=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/script)]
        run=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        require(run.returncode==0,script+' failed: '+run.stderr.decode())
        require(run.stdout==reference,script+' changed receipt')
        mode='optimized' if optimized else 'normal'
        require((HERE/(stem+'.'+mode+'.json')).read_bytes()==reference,'saved receipt disagreement')
        hostile=subprocess.run(command+['--test-failing-gate'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        require(hostile.returncode==1 and b'intentional failure control' in hostile.stderr,'failure gate inactive')
        rows.append({'script':script,'mode':mode,'status':'PASS','failure_control_exit':hostile.returncode,
                     'receipt_sha256':hashlib.sha256(run.stdout).hexdigest()})
print(json.dumps({'status':'PASS','runs':rows},sort_keys=True,indent=2))
