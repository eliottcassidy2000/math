#!/usr/bin/env python3
"""Run only Session247's small new-formula probes; no historical audit."""
import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
jobs=[('marked','probe.py'),('transfer','check.py'),('gluing','check_scale_gluing.py'),('carry','verify.py')]
results=[]
for directory,script in jobs:
 p=root/directory/script
 r=subprocess.run([sys.executable,'-B',script],cwd=p.parent,check=True,capture_output=True,text=True)
 results.append({'lane':directory,'script':script,'exit_code':r.returncode,'stdout_sha256':hashlib.sha256(r.stdout.encode()).hexdigest(),'stdout_bytes':len(r.stdout.encode())})
print(json.dumps({'scope':'Only newly derived formulas and stated finite physical examples. No inherited proof audit, no universal spectrum census.','checks':results},indent=2))
