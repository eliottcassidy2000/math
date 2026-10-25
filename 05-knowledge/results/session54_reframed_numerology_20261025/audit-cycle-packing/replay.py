#!/usr/bin/env python3
import pathlib,hashlib,json,subprocess,sys,os,tempfile,shutil
ROOT=pathlib.Path(__file__).resolve().parent
def ck(x,s):
 if not x:raise RuntimeError(s)
def validate(root):
 m=json.loads((root/'MANIFEST.json').read_text())
 ck({str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}==set(m['files'])|{'MANIFEST.json'},'exact file set')
 for f,d in m['files'].items():ck(hashlib.sha256((root/f).read_bytes()).hexdigest()==d,'digest '+f)
def main():
 validate(ROOT);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 outputs=[]
 for flags in (['-B'],['-B','-O']):
  r=subprocess.run([sys.executable,*flags,str(ROOT/'audit.py')],cwd='/',env=env,capture_output=True,check=True)
  ck(r.stdout==(ROOT/'result.json').read_bytes(),'audit replay');outputs.append(r.stdout)
 ck(outputs[0]==outputs[1]==(ROOT/'result-O.json').read_bytes(),'optimized identity')
 producer=next((ROOT/'inputs').glob('session54-cycle-*'))
 subprocess.run([sys.executable,'-B',str(producer/'replay.py'),'--self-test'],cwd='/',env=env,capture_output=True,check=True)
 if '--self-test' in sys.argv:
  with tempfile.TemporaryDirectory() as tmp:
   c=pathlib.Path(tmp)/'closed';shutil.copytree(ROOT,c);validate(c)
   for kind in ('add','change','remove'):
    f=c/('UNEXPECTED' if kind=='add' else 'result.json');old=f.read_bytes() if f.exists() else None
    if kind=='remove':f.unlink()
    else:f.write_bytes((old or b'')+b'!')
    rejected=False
    try:validate(c)
    except RuntimeError:rejected=True
    ck(rejected,'hostile accepted '+kind)
    if old is None:f.unlink()
    else:f.write_bytes(old)
   # A second relocated copy can rerun with no workspace references.
   r=subprocess.run([sys.executable,'-B',str(c/'replay.py')],cwd='/',env=env,capture_output=True,check=True)
 print('PASS: exact file set, digests, independent normal/-O, bundled producer replay, closed relocation, addition/mutation/deletion hostiles')
if __name__=='__main__':main()
