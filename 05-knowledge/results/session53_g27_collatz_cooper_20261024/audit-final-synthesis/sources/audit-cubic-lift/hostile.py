import json,shutil,subprocess,sys,tempfile
from pathlib import Path
src=Path(__file__).resolve().parent
out={}
with tempfile.TemporaryDirectory() as tmp:
 for flag in [[],['-O']]:
  for kind in ['missing','modified','extra','duplicate','unsafe']:
   dst=Path(tmp)/(str(len(flag))+kind);shutil.copytree(src,dst)
   if kind=='missing':(dst/'normal.json').unlink()
   if kind=='modified':(dst/'REPORT.md').write_text('wrong')
   if kind=='extra':(dst/'nested').mkdir();(dst/'nested'/'extra').write_text('extra')
   if kind=='duplicate':
    p=dst/'MANIFEST.json';s=p.read_text();p.write_text(s.replace('"schema": 1','"schema": 1, "schema": 1'))
   if kind=='unsafe':
    p=dst/'MANIFEST.json';m=json.loads(p.read_text());m['files']['../escape']='0'*64;p.write_text(json.dumps(m))
   r=subprocess.run([sys.executable,*flag,str(dst/'verify.py')],capture_output=True)
   if r.returncode==0:raise RuntimeError('accepted '+kind)
   out[str(len(flag))+kind]='REJECTED'
print(json.dumps(out,sort_keys=True))
