from pathlib import Path
import json, random
import numpy as np
import ctypes,subprocess,tempfile
HERE=Path(__file__).resolve().parent
TMP=tempfile.TemporaryDirectory(prefix='session288_stable_')
SO=Path(TMP.name)/'delta.so'
subprocess.run(['g++','-std=c++17','-O2','-shared','-fPIC',str(HERE/'delta.cpp'),'-o',str(SO)],check=True)
lib=ctypes.CDLL(str(SO))
lib.defect.argtypes=[ctypes.POINTER(ctypes.c_longlong),ctypes.c_int,ctypes.c_int]
lib.defect.restype=ctypes.c_longlong
def delta(rows,O):return lib.defect(rows.ctypes.data_as(ctypes.POINTER(ctypes.c_longlong)),len(rows),O)

def inflate(rows,O,lens):
 typ=[v for v in range(len(rows)) for _ in range(lens.get(v,1))]
 rr=[sum(1<<j for j,w in enumerate(typ) if (v==w and i<j) or rows[v]>>w&1) for i,v in enumerate(typ)]
 return np.array(rr,dtype=np.int64),sum(1<<j for j,w in enumerate(typ) if O>>w&1)
def val(rows,O,lens={}):
 r,o=inflate(rows,O,lens);return int(delta(r,o))
def run(limit=12000):
 bank=json.loads((HERE.parents[1]/'five_stress_connection_atlas_20261112/cyclic_port_receipt.json').read_text())['positive_screen_certificates']
 bank=[r for r in bank if len(r['individually_stable'])>=2 and r['n']>=6]
 rng=random.Random(288);seen=set();tested=0;seeds=0;passing=0
 while tested<limit:
  b=rng.choice(bank);rows=b['rows'];n=len(rows);cut=rng.randrange(1<<n);newO=b['selected']|(rng.randrange(2)<<n);key=(tuple(rows),cut,newO)
  if key in seen:continue
  seen.add(key);tested+=1
  rr=[rows[v]|((1-(cut>>v&1))<<n) for v in range(n)]+[cut]
  if val(rr,newO)!=1:continue
  seeds+=1
  ports=[v for v in range(n+1) if not newO>>v&1 and val(rr,newO,{v:2})==1 and val(rr,newO,{v:3})==1]
  if len(ports)<2:continue
  passing+=1
  from itertools import combinations
  for a,c in combinations(ports,2):
   d=val(rr,newO,{a:2,c:2})
   if d!=1:
    out={'rows':rr,'O':newO,'ports':[a,c],'mixed':d,'tested':tested,'seeds':seeds,'passing':passing}
    print(json.dumps(out),flush=True);(HERE/'candidate.json').write_text(json.dumps(out));return
  if tested%500==0: print(tested,seeds,passing,flush=True)
 print(json.dumps({'tested':tested,'seeds':seeds,'passing':passing,'candidate':None}),flush=True)
if __name__=='__main__':run()
