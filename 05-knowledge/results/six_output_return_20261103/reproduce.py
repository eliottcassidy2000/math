"""Bounded proof replay; optional C++ search commands are listed in AUDIT.md."""
import subprocess,sys,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent
steps=['ray_gate.py','port_rank.py','carry.py','substitution.py','return_audit.py']
records=[]
for step in steps:
    cmd=[sys.executable,'-B']+([] if __debug__ else ['-O'])+[str(HERE/step)]
    p=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8')
    if p.returncode:raise RuntimeError(step+'\n'+p.stdout+'\n'+p.stderr)
    records.append(dict(step=step,output=p.stdout.splitlines()[-2:]));print('PASS',step,flush=True)
outputs=['ray_gate.json','port_rank.json','carry.json','substitution.json','return_audit.json']
hashes={name:hashlib.sha256((HERE/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in outputs}
report=dict(steps=records,hashes=hashes)
(HERE/('verification.normal.json' if __debug__ else 'verification.optimized.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(hashes,sort_keys=True),flush=True)
