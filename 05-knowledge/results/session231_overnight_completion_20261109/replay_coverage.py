"""Regenerate the complete fixed-core context census and its independent audit.
Uses a temporary copy, leaving tracked files unchanged. Needs about 200 MB RAM.
"""
import ast,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
rows=[]
with tempfile.TemporaryDirectory(prefix='session231-coverage-') as tmp:
 root=Path(tmp)
 for name in ('bare_coverage','bare_coverage_audit'):
  shutil.copytree(HERE/name,root/name,ignore=shutil.ignore_patterns('*.bin','*.bitset','__pycache__'))
 scripts=['bare_coverage/census.py','bare_coverage_audit/audit.py','bare_coverage_audit/controls.py','bare_coverage_audit/verify_receipts.py']
 for script in scripts:
  tree=ast.parse((root/script).read_text());need(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'optimizable gate '+script)
  outs=[]
  for opt in (False,True):
   cmd=[sys.executable,'-I','-B']+(['-O'] if opt else [])+[str(root/script)]
   if script.endswith('census.py'):cmd.append('100000000')
   r=subprocess.run(cmd,cwd=(root/script).parent,check=True,capture_output=True);outs.append(r.stdout)
   if 'audit/' in script:
    s=subprocess.run(cmd+['--false-gate'],capture_output=True)
    need(s.returncode!=0 and b'deliberate-false gate sentinel' in s.stderr,'false gate did not fail: '+script)
  need(outs[0]==outs[1],'normal/optimized output '+script)
  rows.append(dict(script=script,status='PASS',stdout_sha256=hashlib.sha256(outs[0]).hexdigest()))
  print('PASS '+script,file=sys.stderr,flush=True)
 for folder in ('bare_coverage','bare_coverage_audit'):
  for name in ('RESULT.json','CONTROLS.json'):
   need((root/folder/name).read_bytes()==(HERE/folder/name).read_bytes(),'frozen receipt '+folder+'/'+name)
 a=(root/'bare_coverage/UNION.bin').read_bytes();b=(root/'bare_coverage_audit/union.bitset').read_bytes()
 need(a==b and len(a)==6250000,'set equality or size')
 sha=hashlib.sha256(a).hexdigest();need(sha=='871cca7430a860219c3658c8672a18728cbecd01989f28e0ad3af4431bb99c77','union hash')
print(json.dumps(dict(status='PASS',normal_and_optimized=True,union_bytes=6250000,union_sha256=sha,checks=rows,scope='Complete target-bounded fixed-core singleton-context language only.'),indent=2))
