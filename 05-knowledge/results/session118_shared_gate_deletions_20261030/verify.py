"""Verify the sealed research packet and rerun independent finite certificates."""
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent

def main():
 m=json.loads((ROOT/'manifest.json').read_text());rows=[]
 for name,digest in m['sha256'].items():
  if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise RuntimeError('File hash mismatch: '+name)
 for name,digest in sorted(m['expected_stdout_sha256'].items()):
  results=[]
  for opts in ([],['-O']):
   r=subprocess.run([sys.executable,*opts,str(ROOT/name)],capture_output=True,check=True)
   json.loads(r.stdout);results.append(r.stdout)
  if results[0]!=results[1]:raise RuntimeError('Optimized output mismatch: '+name)
  actual=hashlib.sha256(results[0]).hexdigest()
  if actual!=digest:raise RuntimeError('Receipt mismatch: '+name)
  rows.append({'script':name,'normal_and_optimized':'PASS','stdout_sha256':actual})
 print(json.dumps({'status':'PASS','checks':rows,'scope':'finite certificates corroborating separate all-size proofs; no full external theorem certification'},indent=2,sort_keys=True))
if __name__=='__main__':main()
