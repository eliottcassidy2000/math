import subprocess,json
from pathlib import Path
root=Path(__file__).resolve().parent
out=[]
for n,codes in [(2,range(2)),(4,range(64)),(6,[(i*1931+7)%32768 for i in range(16)])]:
 for code in codes:
  d=json.loads(subprocess.check_output([str(root/'cubic'),str(n),str(code)]))
  if n==4 and d['nonzero_modular_cubic_coefficients']:raise RuntimeError('order-four vanishing')
  out.append({'order':n,'code':code,**{k:v for k,v in d.items() if k not in ('coefficients','rows')}})
print(json.dumps(out,sort_keys=True,indent=2))
