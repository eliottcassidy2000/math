#!/usr/bin/env python3
"""Fail-closed manifest replay, independent of cwd and external repository."""
import pathlib,hashlib,json,subprocess,sys,os,tempfile,shutil
ROOT=pathlib.Path(__file__).resolve().parent

def ck(p,m):
 if not p:raise RuntimeError(m)
def validate(root):
 manifest=json.loads((root/'MANIFEST.json').read_text())
 actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
 expected=set(manifest['files'])|{'MANIFEST.json'}
 ck(actual==expected,'file-set mismatch: '+str(sorted(actual^expected)))
 for name,digest in manifest['files'].items():
  ck(hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,'digest mismatch: '+name)
 return manifest

def main():
 validate(ROOT)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 for script,expected in [('verify.py','result.json'),('algebra.py','algebra.json'),('two_sheet.py','two-sheet.json')]:
  outputs=[]
  for flags in (['-B'],['-B','-O']):
   result=subprocess.run([sys.executable,*flags,str(ROOT/script)],cwd='/',env=env,check=True,capture_output=True)
   ck(result.stdout==(ROOT/expected).read_bytes(),'replay mismatch: '+script+' '+str(flags));outputs.append(result.stdout)
  ck(outputs[0]==outputs[1],'optimized mismatch')
 if '--self-test' in sys.argv:
  with tempfile.TemporaryDirectory() as tmp:
   copy=pathlib.Path(tmp)/'packet';shutil.copytree(ROOT,copy)
   validate(copy)
   (copy/'UNEXPECTED').write_text('hostile')
   try:validate(copy)
   except RuntimeError:pass
   else:raise RuntimeError('file-set hostile accepted')
   (copy/'UNEXPECTED').unlink()
   with (copy/'result.json').open('a') as f:f.write(' ')
   try:validate(copy)
   except RuntimeError:pass
   else:raise RuntimeError('digest hostile accepted')
 print('PASS: manifest, exact file set, arbitrary-directory normal/-O replay'+(', mutation hostiles' if '--self-test' in sys.argv else ''))
if __name__=='__main__':main()
