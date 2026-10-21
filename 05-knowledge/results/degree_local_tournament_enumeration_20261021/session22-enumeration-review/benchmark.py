"""Matched-telemetry sensitivity check, not a production replacement."""
from pathlib import Path
import importlib.util,json,time,statistics,hashlib
ROOT=Path(__file__).resolve().parents[1]
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
source=(ROOT/'session22-enumeration/baseline.py').read_text()
source=source.replace('                    response_coordinate_updates += count','                    # Auxiliary coordinate telemetry omitted in matched control.')
source=source.replace('        tuple_cells = sum(max(0, (min(ell - 2, bound - d) - start) // 2 + 1)\n                          for d, _ in states)','        tuple_cells = 0 # Auxiliary response telemetry omitted.')
source=source.replace('sum(r.bit_length() for _, r in states)','0')
p=Path(__file__).with_name('baseline_light.py');p.write_text(source)
old=load('light',p);new=load('new',ROOT/'session22-enumeration/local_scale.py');rows=[]
for N,core in ((100,False),(100,True),(200,False)):
 times={'baseline_light':[],'local':[]}
 for rep in range(3):
  answers={}; states={}
  pairs=[('baseline_light',old.packed_dp),('local',new.local_dp)]
  for key,fn in (pairs if rep%2==0 else pairs[::-1]):
   t=time.perf_counter(); answers[key],stats=fn(N,core=core);times[key].append(time.perf_counter()-t); states[key]=[(r["cycle_length"],r["states"],r["transitions"]) for r in stats["rows"]]
  if answers['baseline_light']!=answers['local'] or states['baseline_light']!=states['local']:raise RuntimeError('different coefficient vectors')
 rows.append({'bound':N,'core':core,'seconds':times,'median_speedup':statistics.median(times['baseline_light'])/statistics.median(times['local'])})
print(json.dumps({'description':'Baseline auxiliary per-transition coordinate telemetry and per-state tuple/bit telemetry removed; setup, arithmetic, rows and output retained. Shared-host sensitivity experiment, not attribution of all remaining differences.','baseline_light_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'benchmarks':rows},indent=2))
