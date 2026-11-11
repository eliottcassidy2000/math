"""Run new-claim engines only, comparing exact saved and optimized outputs."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parent
PAIRS=[('integral/probe.py','integral/results.json'),('moments/probe.py','moments/probe.out'),('niche/probe.py','niche/results.json'),('niche/cost_probe.py','niche/cost_results.json'),('niche/cube_probe.py','niche/cube_results.json'),('tournaments/probe.py','tournaments/probe.out'),('lrc/probe.py','lrc/results.json'),('lrc/asymmetric.py','lrc/asymmetric.json'),('lrc/hostile.py','lrc/hostile.json')]
def run(pair):
 script,output=pair;expected=(ROOT/output).read_bytes()
 results=[]
 for flags in (['-B'],['-B','-O']):
  r=subprocess.run([sys.executable,*flags,str(ROOT/script)],cwd=ROOT,capture_output=True,check=True)
  if r.stdout!=expected:raise RuntimeError('Output mismatch: '+script+' '+str(flags))
  results.append(hashlib.sha256(r.stdout).hexdigest())
 return {'script':script,'saved_output':output,'normal_and_optimized_match':True,'sha256':results[0]}
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,PAIRS))
 print(json.dumps({'scope':'new claims only; inherited large scans not repeated','engines':results},indent=2,sort_keys=True))
