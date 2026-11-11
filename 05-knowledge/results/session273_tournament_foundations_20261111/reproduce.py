"""Replay this packet's structural controls, without external dependencies."""
from pathlib import Path
import subprocess,sys,json,hashlib
root=Path(__file__).resolve().parent
jobs=[('orbits/control.py','orbits/control.json'),('microscopic/check.py','microscopic/OUTPUT.json'),('wildcard/clock_lift.py','wildcard/controls.json'),('interfaces/check.py','interfaces/RESULTS.json')]
results=[]
for src,out in jobs:
    a=subprocess.run([sys.executable,str(root/src)],cwd=root,capture_output=True,check=True).stdout
    b=subprocess.run([sys.executable,'-O',str(root/src)],cwd=root,capture_output=True,check=True).stdout
    if a!=b:raise RuntimeError('normal/optimized output mismatch: '+src)
    observed=json.loads(a);saved=json.loads((root/out).read_text())
    if observed!=saved:raise RuntimeError('saved output mismatch: '+src)
    results.append(dict(script=src,status='PASS',normal_optimized_byte_equal=True,saved_json_equal=True,source_sha256=hashlib.sha256((root/src).read_bytes()).hexdigest()))
print(json.dumps(dict(status='PASS',controls=results,scope='structural controls, not global spectrum enumeration; optimized replay is packaging consistency, not independent proof'),indent=2))
