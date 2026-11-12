#!/usr/bin/env python3
"""Fresh normal/optimized/stored comparison of every proof engine and scoped screen."""
from pathlib import Path
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent
JOBS=[
 ('bridge/verify.py','bridge/REPLAY.json',None,None),
 ('lengths/verify.py','lengths/REPLAY.json','lengths/REPLAY.json','lengths/REPLAY.optimized.json'),
 ('wildcard/verify.py','wildcard/REPLAY.json','wildcard/REPLAY.json','wildcard/REPLAY.optimized.json'),
 ('stable/verify.py','stable/replay.json','stable/replay.json','stable/replay.optimized.json'),
 ('stable/cones.py','stable/cones.json','stable/cones.json','stable/cones.optimized.json'),
 ('stable/normalforms.py','stable/normalforms.json',None,None),
 ('stable/explore.py','stable/explore.out',None,None),
 ('stable/clocks.py','stable/clocks.json',None,None),
 ('audit/three_clock.py','audit/three_clock.json',None,None),
 ('audit/independent.py','audit/REPLAY.json',None,None),
]
def main():
 out=[]
 for script,expected,normal_path,opt_path in JOBS:
  before=(HERE/expected).read_bytes();results=[]
  for opt,p in [(False,normal_path),(True,opt_path)]:
   cmd=[sys.executable]+(['-O'] if opt else [])+[str(HERE/script)]
   result=subprocess.run(cmd,cwd=HERE,capture_output=True,check=True)
   actual=(HERE/p).read_bytes() if p else result.stdout
   if actual!=before:raise RuntimeError(f'{script}: stored/fresh mismatch optimized={opt}')
   results.append(actual)
  if results[0]!=results[1]:raise RuntimeError(f'{script}: optimized mismatch')
  out.append({'script':script,'normal_optimized_stored':'MATCH','receipt_sha256':hashlib.sha256(before).hexdigest()})
 report={'status':'PASS','jobs':out,'scope':'Ten declared proof-engine/screen comparisons. Inherited full census is not re-enumerated. Optional wildcard discovery scans are not part of this release replay.'}
 text=json.dumps(report,sort_keys=True,indent=2)+'\n';(HERE/'REPLAY.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
