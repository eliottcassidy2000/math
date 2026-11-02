#!/usr/bin/env python3
"""Portable release replay. Inputs are frozen; no generator install is needed."""
from pathlib import Path
import hashlib,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def run(args):
 p=subprocess.run(args,cwd=HERE,text=True,capture_output=True)
 need(p.returncode==0,str(args)+'\n'+p.stdout+'\n'+p.stderr)
 return p.stdout
def manifest():
 f=HERE/'MANIFEST.sha256'
 if f.exists():
  for line in f.read_text().splitlines():
   digest,name=line.split('  ',1)
   need(hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,'manifest '+name)

def main():
 manifest()
 checks=[]
 generator=json.loads((HERE/'GENERATOR.json').read_text())
 for r in generator['inputs']:
  b=(HERE/r['path']).read_bytes();need(hashlib.sha256(b).hexdigest()==r['sha256'],'frozen input hash')
  need(len(b.splitlines())==r['classes'],'frozen class count')
 print('frozen class inputs PASS',flush=True)
 for script,out in [('valuation.py','valuation.json'),('exchange.py','exchange.json'),('escape_audit.py','escape_audit.json'),('regular_balance.py','regular_balance.json')]:
  run([sys.executable,script]);a=(HERE/out).read_bytes()
  run([sys.executable,'-O',script]);b=(HERE/out).read_bytes()
  need(a==b,script+' normal/-O agreement')
  checks.append(dict(script=script,normal_optimized_equal=True,output=out,sha256=hashlib.sha256(a).hexdigest()))
  print(script,'PASS',flush=True)
 # Replay the full producer using frozen inputs and compare both algorithms.
 run([sys.executable,'research.py','--max-order','8'])
 producer=json.loads((HERE/'module_census.json').read_text())
 with tempfile.TemporaryDirectory(prefix='prime-inverse-audit-') as td:
  for stem in ['local_escape','module_audit']:
   expected=None
   for opt in ['-O0','-O2']:
    exe=str(Path(td)/(stem+opt))
    run(['c++',opt,'-std=c++17','-Wall','-Wextra','-Werror',str(HERE/(stem+'.cpp')),'-o',exe])
    args=[exe]+([str(HERE/'inputs')] if stem=='module_audit' else [])
    b=run(args).encode()
    if expected is not None:need(b==expected,stem+' O0/O2 agreement')
    expected=b
   need(json.loads(expected)==json.loads((HERE/(stem+'.json')).read_text()),stem+' recorded output')
   checks.append(dict(script=stem+'.cpp',O0_O2_equal=True,sha256=hashlib.sha256(expected).hexdigest()))
   print(stem,'PASS',flush=True)
 a=json.loads((HERE/'module_audit.json').read_text())
 need(a['cover_vectors']==[r['covers'] for r in producer['blocks']],'independent cover vectors')
 total=0
 for row in a['orders']:
  ref=producer['by_order'][str(row['n'])]
  need(row['strong_classes']==ref['strong_classes'],'independent strong counts')
  need(row['minima']==[r['min'] for r in ref['minima']],'independent marked minima')
  total+=row['n']*row['strong_classes']
 need(total==50782,'marked input universe')
 need(a['symmetric_7p_below_30625']==sorted(map(int,producer['threshold_30625_symmetric_7p'])),'independent prime targets')
 escape=json.loads((HERE/'escape_audit.json').read_text());local=json.loads((HERE/'local_escape.json').read_text())
 need(escape['distance_histogram']==local['strict_minimum_distance_histogram'],'independent move distances')
 need(escape['native_complement_distance_histogram']=={'2':480,'3':288},'native distances')
 need(all(r['barrier_walk']['counts']==[15,13,23,17] for r in escape['orbit_representatives']),'uniform count trajectory')
 out=dict(status='PASS',checks=checks,marked_strong_quotients=total,independent_prime_targets=12,
          input_generation_scope='Frozen nauty complete-class banks, not two independent class generators',
          audit_scope='Same-agent proof review with independent counting/control algorithms; no subagent or formal audit')
 (HERE/'REPLAY.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');manifest();print('all release gates PASS')
if __name__=='__main__':main()
