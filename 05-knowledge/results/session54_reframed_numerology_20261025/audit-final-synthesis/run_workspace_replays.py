from pathlib import Path
import subprocess,concurrent.futures,json,sys,os,tempfile,shutil
A=Path(__file__).resolve().parent;R=A.parent
def run(name):
 p=R/('session54-'+name)/'replay.py';s=subprocess.run([sys.executable,'-B',str(p),'--self-test'],cwd='/tmp',capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),timeout=1200)
 if s.returncode:raise RuntimeError(name+s.stdout+s.stderr)
 (A/(name+'-replay.txt')).write_text(s.stdout);return name
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(run,['audit-larger-paley','audit-endpoint-response','diagrams']))
with tempfile.TemporaryDirectory() as t:
 p=Path(t)/'packet';shutil.copytree(R/'session54-audit-paley-endpoint',p)
 for opt in [[],['-O']]:
  for script,fixture in [('verify_manifest.py',None),('audit.py','normal.json'),('hostiles.py','hostiles.normal.json')]:
   s=subprocess.run([sys.executable,'-B',*opt,str(p/script)],cwd=t,capture_output=True,check=True,timeout=180)
   if fixture and s.stdout!=(p/fixture).read_bytes():raise RuntimeError('endpoint fixture '+script)
 victim=p/'AUDIT.md';old=victim.read_bytes()
 for kind in ['corrupt','missing','extra','nested-extra']:
  v=p/('inputs/extra' if kind=='nested-extra' else 'extra')
  if kind=='corrupt':victim.write_bytes(old+b'X')
  elif kind=='missing':victim.unlink()
  else:v.write_text('x')
  for opt in [[],['-O']]:
   s=subprocess.run([sys.executable,'-B',*opt,str(p/'verify_manifest.py')],cwd=t,capture_output=True)
   if s.returncode==0:raise RuntimeError('accepted '+kind)
  if kind in ['corrupt','missing']:victim.write_bytes(old)
  else:v.unlink()
rows.append('audit-paley-endpoint: relocated normal/-O exact fixtures and four closure mutations')
(A/'workspace-replays.json').write_text(json.dumps({'status':'PASS','checks':rows},indent=2)+'\n');print(rows)
