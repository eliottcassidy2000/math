import json,pathlib,tempfile,subprocess,sys
r=pathlib.Path(__file__).resolve().parent
source=(r/'independent.py').read_text();original=json.loads((r/'inputs'/'producer-results.json').read_text());out=[]
for name in ['wrong_N3','omit_five_d1','wrong_passing_count','wrong_affine_normalization']:
 with tempfile.TemporaryDirectory() as td:
  p=pathlib.Path(td);(p/'inputs').mkdir();d=json.loads(json.dumps(original));code=source
  if name=='wrong_N3':d['bounds'][2]['N']+=1
  if name=='omit_five_d1':d['bounds'][0]['prime_exponents'].pop('5')
  if name=='wrong_passing_count':d['bounds'][5]['passing_count']+=1
  if name=='wrong_affine_normalization':code=code.replace('a,b=85995,-52515','a,b=171990,-52515')
  (p/'inputs'/'producer-results.json').write_text(json.dumps(d));(p/'independent.py').write_text(code)
  for opt in [[],['-O']]:
   z=subprocess.run([sys.executable,*opt,str(p/'independent.py')],capture_output=True,text=True)
   if z.returncode==0:raise RuntimeError('mutation accepted '+name)
   out.append({'case':name,'optimized':bool(opt),'rejected':True,'reason':z.stderr.splitlines()[-1]})
print(json.dumps(out,indent=2))
