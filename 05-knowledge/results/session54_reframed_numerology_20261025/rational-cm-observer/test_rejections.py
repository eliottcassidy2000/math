#!/usr/bin/env python3
import pathlib,tempfile,subprocess,sys,json
source=pathlib.Path(__file__).with_name('check.py').read_text()
cases=[('wrong_odd_formula','//math.gcd(4,p-1)','//math.gcd(2,p-1)','odd prime-power formula'),('wrong_N3','(65520,{2:4','(65521,{2:4','degree-three control'),('wrong_affine_CM_control','a,b=85995,-52515','a,b=85995,-52514','H-15 affine positive control')]
out=[]
with tempfile.TemporaryDirectory() as td:
 for name,old,new,reason in cases:
  if source.count(old)!=1:raise RuntimeError('mutation target count')
  p=pathlib.Path(td)/(name+'.py');p.write_text(source.replace(old,new))
  r=subprocess.run([sys.executable,*(['-O'] if not __debug__ else []),str(p)],capture_output=True,text=True,timeout=15)
  if not r.returncode or reason not in r.stderr:raise RuntimeError('mutation not rejected: '+name+r.stderr)
  out.append({'case':name,'rejected':True,'reason':reason})
print(json.dumps(out,indent=2))
