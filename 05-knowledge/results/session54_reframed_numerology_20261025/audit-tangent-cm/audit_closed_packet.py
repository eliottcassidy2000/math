import pathlib,tempfile,shutil,subprocess,hashlib,json,sys
source=pathlib.Path(sys.argv[1] if len(sys.argv)>1 else 'session54-tangent-cm-publication').resolve()
def req(c,m):
 if not c:raise RuntimeError(m)
results=[]
with tempfile.TemporaryDirectory() as td:
 root=pathlib.Path(td)/'packet';shutil.copytree(source,root)
 def run(script,opt=False):
  return subprocess.run([sys.executable,*(['-O'] if opt else []),str(root/script)],cwd=td,capture_output=True,text=True,timeout=90)
 for opt in [False,True]:
  for script,reference in [('verify_manifest.py',None),('check.py','results.normal.json'),('test_rejections.py','rejections.normal.json')]:
   r=run(script,opt);req(r.returncode==0,script+' failed '+r.stderr)
   if reference:req(r.stdout==(root/reference).read_text(),'output mismatch '+script)
   results.append({'script':script,'optimized':opt,'status':'PASS'})
 victim=root/'check.py';original=victim.read_bytes()
 for kind in ['corrupt','missing','extra','nested-extra']:
  if kind=='corrupt':victim.write_bytes(original+b'\n# mutation\n')
  if kind=='missing':victim.unlink()
  if kind=='extra':(root/'unexpected-file.txt').write_text('not in manifest')
  if kind=='nested-extra':(root/'inputs'/'unexpected-file.txt').write_text('not in manifest')
  for opt in [False,True]:
   r=run('verify_manifest.py',opt);req(r.returncode!=0,'manifest accepted '+kind)
   results.append({'hostile':kind,'optimized':opt,'status':'REJECTED'})
  if kind=='extra':(root/'unexpected-file.txt').unlink()
  elif kind=='nested-extra':(root/'inputs'/'unexpected-file.txt').unlink()
  else:victim.write_bytes(original)
print(json.dumps({'status':'CLOSED REPLAY PASS','manifest_sha256':hashlib.sha256((source/'MANIFEST.json').read_bytes()).hexdigest(),'checks':results},indent=2))
