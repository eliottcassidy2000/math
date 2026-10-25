#!/usr/bin/env python3
"""Normal/optimized equality plus real CLI rejection of tampered payloads."""
import hashlib,json,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
BASE=[str(HERE/'verify_and_search.py'),'--bound','100000','--independent-bound','2000','--verify']
def run(extra=(),optimized=False):
 cmd=[sys.executable]+(['-O'] if optimized else [])+BASE+list(extra)
 return subprocess.run(cmd,capture_output=True,text=True)
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 a,b=run(),run(optimized=True)
 require(a.returncode==b.returncode==0,'normal or optimized replay failed')
 require(a.stdout==b.stdout,'normal and optimized replay outputs differ')
 results={'normal_exit_code':a.returncode,'optimized_exit_code':b.returncode,'identical_stdout':True,'stdout_sha256':hashlib.sha256(a.stdout.encode()).hexdigest(),'mutation_rejections':{}}
 data=json.loads((HERE/'certificates.json').read_text())
 with tempfile.TemporaryDirectory(prefix='parker-verifier-') as tmp:
  for name in ['changed_content','missing_key','extra_key','missing_file']:
   d=json.loads(json.dumps(data));p=Path(tmp)/(name+'.json')
   if name=='changed_content':d['integer_search']['solutions']=[[0]*9]
   elif name=='missing_key':d.pop('schema')
   elif name=='extra_key':d['unexpected']=True
   if name!='missing_file':p.write_text(json.dumps(d))
   r=run(['--output',str(p)],optimized=True)
   require(r.returncode!=0,'accepted mutation: '+name)
   expected='FileNotFoundError' if name=='missing_file' else 'certificate data mismatch'
   require(expected in r.stderr,'mutation rejected for unexpected reason: '+name)
   results['mutation_rejections'][name]={'exit_code':r.returncode,'expected_failure_seen':True}
 (HERE/'verifier_checks.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
 print(json.dumps(results,indent=2,sort_keys=True))
if __name__=='__main__':main()
