#!/usr/bin/env python3
import pathlib,json,hashlib,subprocess,sys,tempfile,shutil,os
R=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def validate(r):
 m=json.loads((r/'MANIFEST.json').read_text())['files']
 ck({str(p.relative_to(r)) for p in r.rglob('*') if p.is_file()}==set(m)|{'MANIFEST.json'},'file set')
 for f,h in m.items():ck(hashlib.sha256((r/f).read_bytes()).hexdigest()==h,'digest '+f)
validate(R)
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
for opts,out in [(['-B'],'result.json'),(['-B','-O'],'result-O.json')]:
 p=subprocess.run([sys.executable,*opts,str(R/'audit.py')],cwd='/',env=env,capture_output=True,check=True)
 ck(p.stdout==(R/out).read_bytes(),'independent replay '+out)
ck((R/'result.json').read_bytes()==(R/'result-O.json').read_bytes(),'optimization invariance')
p=subprocess.run([sys.executable,'-B',str(R/'input/replay.py'),'--self-test'],cwd='/',env=env,capture_output=True,check=True)
for mutation in ('add','remove','content','input-manifest'):
 with tempfile.TemporaryDirectory() as tmp:
  r=pathlib.Path(tmp)/'audit';shutil.copytree(R,r);validate(r)
  if mutation=='add':(r/'EXTRA').write_text('hostile')
  elif mutation=='remove':(r/'result.json').unlink()
  elif mutation=='content':(r/'result.json').write_bytes((r/'result.json').read_bytes()+b' ')
  else:(r/'input/MANIFEST.json').write_bytes((r/'input/MANIFEST.json').read_bytes()+b' ')
  try:validate(r)
  except RuntimeError:pass
  else:raise RuntimeError('hostile accepted '+mutation)
print('PASS: closed audit manifest; independent normal/-O; pinned producer replay; four mutation hostiles')
