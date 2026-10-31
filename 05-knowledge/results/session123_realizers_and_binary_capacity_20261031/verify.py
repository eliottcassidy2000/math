"""Reproduce every declared finite check, normal and optimized, standard library."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def sha(raw):return hashlib.sha256(raw).hexdigest()
def need(ok,message):
 if not ok:raise RuntimeError(message)
def pipeline(opt,m,checkpoint):
 flags=['-B']+(['-O'] if opt else []);rows=[]
 for row in m['checks']:
  out=subprocess.run([sys.executable,*flags,str(ROOT/row['script'])],capture_output=True,check=True).stdout
  need(out==(ROOT/row['receipt']).read_bytes(),'output receipt mismatch: '+row['script'])
  if checkpoint:
   name=('optimized_' if opt else 'normal_')+row['script'].replace('/','_')+'.out';(checkpoint/name).write_bytes(out)
  rows.append(dict(script=row['script'],stdout_sha256=sha(out),status='PASS'))
  print('PASS',('optimized' if opt else 'normal'),row['script'],file=sys.stderr,flush=True)
 return rows

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--checkpoint-directory',type=Path);a=ap.parse_args()
 if a.checkpoint_directory:a.checkpoint_directory.mkdir(parents=True,exist_ok=True)
 m=json.loads((ROOT/'manifest.json').read_text())
 for name,h in m['sha256'].items():need(sha((ROOT/name).read_bytes())==h,'sealed file mismatch: '+name)
 with ThreadPoolExecutor(max_workers=2) as pool:
  x=pool.submit(pipeline,False,m,a.checkpoint_directory);y=pool.submit(pipeline,True,m,a.checkpoint_directory);normal=x.result();opt=y.result()
 need(normal==opt,'normal/optimized mismatch')
 print(json.dumps(dict(status='PASS',normal_and_optimized='identical',checks=normal,scope='All declared finite controls replayed. Uniform analytic proofs are separately audited. Classical realizability/Lenstra inputs and the inherited Session122 family theorem are explicit dependencies, not re-proved by computation. No Lean or external human review claimed.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
