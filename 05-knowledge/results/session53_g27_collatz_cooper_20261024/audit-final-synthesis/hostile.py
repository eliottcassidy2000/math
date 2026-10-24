from pathlib import Path
import tempfile,shutil,subprocess,sys,json
P=Path(__file__).resolve().parent
out={}
for label in ['source-corruption','missing-source','invalid-link','arithmetic-weight']:
 with tempfile.TemporaryDirectory() as t:
  q=Path(t)/'packet';shutil.copytree(P,q)
  if label=='source-corruption':
   f=q/'sources/dyadic-mod4/REPORT.md';f.write_bytes(f.read_bytes()+b'wrong')
  elif label=='missing-source':(q/'sources/dyadic-mod4/REPORT.md').unlink()
  elif label=='invalid-link':
   f=q/'REVIEWED_README.md';f.write_text(f.read_text().replace('(dyadic-mod4/REPORT.md)','(missing/REPORT.md)'))
  else:
   f=q/'check.py';s=f.read_text();old='weights=[1,2,4,3,3,6,0,0,0]'
   if old not in s:raise RuntimeError('mutation target missing')
   f.write_text(s.replace(old,'weights=[1,2,5,3,3,6,0,0,0]'))
  errors=[]
  for flag in [[],['-O']]:
   r=subprocess.run([sys.executable,*flag,str(q/'check.py')],capture_output=True,text=True,cwd='/tmp')
   if r.returncode==0:raise RuntimeError('survived '+label)
   errors.append(r.stderr.splitlines()[-1])
  out[label]=errors
print(json.dumps(out,indent=2))
