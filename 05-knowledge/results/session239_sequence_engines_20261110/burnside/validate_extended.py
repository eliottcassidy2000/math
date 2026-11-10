#!/usr/bin/env python3
"""Independent external table, inherited compiled engine, residue and boundary checks."""
import sys,subprocess,tempfile,time,json,importlib.util
from pathlib import Path
from statistics import median
from engines import integer_response,require
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
reference=HERE.parents[1]/'tournament_future_boundary_20261020/a000568-reference.txt'
ref={int(a):int(b) for line in reference.read_text().splitlines() if line and not line.startswith('#') for a,b in [line.split()]}
a,_=integer_response(max(ref));require(all(a[n]==v for n,v in ref.items()),'external reference')
report={'reference_terms':len(ref),'reference_max_n':max(ref),'boundaries':[],'benchmarks':[]}
for n in range(13):
 a,_=integer_response(n)
 for b in [1,2,4,8,16]:
  r,_=integer_response(n,binary_precision=b);require(r==[x%(1<<b) for x in a],'boundary residue')
 report['boundaries'].append(n)
source=subprocess.run(['git','show','HEAD:04-computation/a000568_gmp_enum_v4.c'],cwd=ROOT,capture_output=True,text=True,check=True).stdout
with tempfile.TemporaryDirectory(prefix='s239-gmp-') as tmp:
 c=Path(tmp)/'inherited.c';exe=Path(tmp)/'inherited';c.write_text(source)
 subprocess.run(['gcc','-O3',str(c),'-lgmp','-pthread','-o',str(exe)],check=True)
 for n in [80,100,120]:
  durations=[]
  for _ in range(3):
   t=time.perf_counter();p=subprocess.run([str(exe),str(n),'1'],capture_output=True,text=True,check=True);durations.append(time.perf_counter()-t)
  old=int(p.stdout.splitlines()[0].split(' = ')[1]);t=time.perf_counter();new,stats=integer_response(n);dt=time.perf_counter()-t
  require(old==new[-1],'compiled independent check')
  t=time.perf_counter();r,rs=integer_response(n,binary_precision=16);rt=time.perf_counter()-t
  require(r==[x%65536 for x in new],'binary precision check')
  row={'n':n,'compiled_single_term_seconds':median(durations),'integer_full_prefix_seconds':dt,'integer_prefix_stats':stats,'binary16_full_prefix_seconds':rt,'binary16_stats':rs}
  report['benchmarks'].append(row);print(json.dumps(row),flush=True)
HERE.joinpath('extended.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: external table, compiled engine, binary residues and boundaries')
