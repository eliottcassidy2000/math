"""Portable audit closure and independent relocated producer verification."""
import pathlib,json,hashlib,subprocess,sys,tempfile,os,shutil
D=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def validate(root):
 m=json.loads((root/'MANIFEST.json').read_text())['files']
 ck({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}==set(m)|{'MANIFEST.json'},'exact file set')
 for n,h in m.items():ck(hashlib.sha256((root/n).read_bytes()).hexdigest()==h,'digest '+n)
def run(root,script,opt=()):
 return subprocess.run([sys.executable,'-B',*opt,str(root/script)],cwd='/',env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,check=True).stdout
def main():
 validate(D)
 with tempfile.TemporaryDirectory() as t:
  root=pathlib.Path(t)/'relocated';root.mkdir()
  pm=D/'producer'/'MANIFEST.json'
  ck(hashlib.sha256(pm.read_bytes()).hexdigest()=='63c4f7bab73cefde83192b5ca65e71a3b139cb4cb674111b3dee45ba2b3bd35e','producer anchor')
  m=json.loads(pm.read_text())['files']
  # Relocate from only manifest-listed paths, not a blanket directory copy.
  for n in [*m,'MANIFEST.json']:
   p=root/n;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(D/'producer'/n,p)
  validate(root)
  pins=json.loads((root/'source-pins.json').read_text())
  for n,h in pins['inputs'].items():ck(hashlib.sha256((root/'inputs'/n).read_bytes()).hexdigest()==h,'inherited source '+n)
  for opt in ((),('-O',)):
   for script,expected in [('verify.py','result.json'),('extend_converse.py','converse.json')]:
    ck(run(root,script,opt)==(root/expected).read_bytes(),'producer replay '+script)
   ck(run(D,'audit.py',opt)==(D/'result.json').read_bytes(),'independent replay')
   ck(run(D,'hostile_controls.py',opt)==(D/'hostile-controls.json').read_bytes(),'semantic mutations')
  for name in ('extra','digest','missing'):
   p=root/'result.json';old=p.read_bytes()
   if name=='extra':(root/'UNEXPECTED').write_text('hostile')
   elif name=='digest':p.write_bytes(old+b' ')
   else:p.unlink()
   try:validate(root)
   except RuntimeError:pass
   else:raise RuntimeError('accepted '+name)
   if name=='extra':(root/'UNEXPECTED').unlink()
   else:p.write_bytes(old)
  validate(root)
 validate(D)
 print('PASS: exact audit closure; manifest-only producer relocation; all source pins; producer and independent normal/-O replay; semantic and extra/digest/missing-file hostiles')
if __name__=='__main__':main()
