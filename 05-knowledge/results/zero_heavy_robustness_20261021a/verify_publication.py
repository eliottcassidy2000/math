#!/usr/bin/env python3
"""Replay unchanged packets against exact, explicitly historical Git inputs."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,tempfile
H=Path(__file__).resolve().parent
def need(v,s):
 if not v:raise RuntimeError(s)
def sh(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);a=p.parse_args();m=json.loads((H/'manifest.json').read_text());v=json.loads((H/'PUBLICATION_PROVENANCE.json').read_text())
 need(sh((H/'manifest.json').read_bytes())==v['original_source_manifest_sha256'],'original manifest')
 need(sh((H/'review/manifest.json').read_bytes())==v['independent_audit_manifest_sha256'],'independent manifest')
 for f in v['incoming_pins']:
  need(sh((a.repo/f['required_destination']).read_bytes())==f['sha256'],'current correction pin '+f['required_destination'])
 with tempfile.TemporaryDirectory(prefix='zero-heavy-historical-') as d:
  root=Path(d)
  for n,h in m['dependencies'].items():
   b=subprocess.check_output(['git','-C',str(a.repo),'show',m['base_commit']+':'+n]);need(sh(b)==h,'historical input '+n);q=root/n;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(b)
  env=dict(os.environ,MATH_REPO=str(root),PYTHONDONTWRITEBYTECODE='1');runs=[]
  for cmd,receipt in [([sys.executable,'-B',str(H/'verify.py')],H/'verify.out'),([sys.executable,'-B',str(H/'review/verify.py'),'--source',str(H),'--repo',str(root)],H/'review/verification.json')]:
   z=subprocess.run(cmd,cwd=H,env=env,capture_output=True);need(z.returncode==0,z.stderr.decode());need(z.stdout==receipt.read_bytes(),'exact replay '+str(receipt));runs.append({'receipt':str(receipt.relative_to(H)),'sha256':sh(z.stdout)})
 print(json.dumps({'status':'PASS','original_packets_unchanged':True,'historical_commit':m['base_commit'],'historical_dependencies':len(m['dependencies']),'replays':runs},indent=2,sort_keys=True))
if __name__=='__main__':main()
