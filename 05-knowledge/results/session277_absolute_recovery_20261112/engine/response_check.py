from pathlib import Path
import json,time,sys,importlib.util
from dyadic_completion import residue_response,integer_interval,recover_from_interval
sys.set_int_max_str_digits(0)
R=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('prior',R/'session270_chanpak_enumeration_20261111/enumeration/counter.py'); prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
values=list(map(int,json.loads((R/'session270_chanpak_enumeration_20261111/enumeration/prefix.json').read_text())['values']))
controls=0
for n in range(41):
 for b in (0,1,2,3,5,8,13,21,55,144):
  r,_=residue_response(n,b);assert r==values[n]%(1<<b);controls+=1
rows=[]
for n in (40,64,80,100,120,150,200):
 lo,hi=integer_interval(n,3);b=(hi-lo).bit_length()
 start=time.perf_counter();r,stats=residue_response(n,b); elapsed=time.perf_counter()-start
 assert recover_from_interval(lo,hi,r,b)==values[n]
 row=dict(n=n,bits=b,count_bits=values[n].bit_length(),seconds=elapsed,**stats)
 print(row,flush=True);rows.append(row)
(Path(__file__).parent/'response_results.json').write_text(json.dumps(dict(controls=controls,rows=rows,scope='single run; existing exact Session270 prefix comparator; modular raw recurrence inherits Session270 response mechanism'),indent=2)+'\n')
