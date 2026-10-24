#!/usr/bin/env python3
"""Full bounded two-mode replay into a caller-selected directory."""
import argparse,json,pathlib,subprocess,sys,time
here=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output-directory',required=True);a=p.parse_args()
out=pathlib.Path(a.output_directory).resolve();out.mkdir(parents=True,exist_ok=True)
records=[]
for mode,flags in [('O2',['-O2']),('O3',['-O3','-DNDEBUG'])]:
 binary=out/f'regular_trace_{mode}'
 subprocess.run(['g++',*flags,'-std=c++17','-Wall','-Wextra','-pedantic',str(here/'regular_trace.cpp'),'-o',str(binary)],check=True)
 for m,name,source in [(5,'seven','parents-5.txt'),(7,'nine','parents-7.txt'),(9,'eleven','parents-nine.txt')]:
  prefix=out/f'{name}-{mode}';start=time.monotonic()
  run=subprocess.run([str(binary),str(m),str(here/'sources'/source),str(prefix)],text=True,capture_output=True,check=True)
  for suffix in ('.json','-parents.txt','-states.tsv','-components.tsv'):
   if pathlib.Path(str(prefix)+suffix).read_bytes()!=(here/f'{name}-O2{suffix}').read_bytes():raise RuntimeError(f'nonidentical replay {name} {mode} {suffix}')
  records.append({'order':m+2,'mode':mode,'seconds':round(time.monotonic()-start,4),'stdout':run.stdout,'stderr':run.stderr,'byte_identical':True})
for optimized in (False,True):
 target=out/('degree-mass-optimized.json' if optimized else 'degree-mass.json')
 subprocess.run([sys.executable,*(['-O'] if optimized else []),str(here/'degree_mass.py'),'--output',str(target)],check=True)
 if target.read_bytes()!=(here/'degree-mass.json').read_bytes():raise RuntimeError('degree DP replay mismatch')
(out/'replay.json').write_text(json.dumps({'status':'PASS','replays':records},indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','replay_directory':str(out),'runs':len(records)}))
