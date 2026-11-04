from pathlib import Path
import json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
def old(Q):
 T=[0]*20
 for i in range(5):
  for u in range(4):
   T[4*i+u]=((15^((1<<(u+1))-1))<<(4*i))|sum(15<<(4*j) for j in range(5) if Q[i]>>j&1)
 return T
def phys(a):return sum(15<<(4*i) for i in range(5) if a>>i&1)
R=json.loads((ROOT/'quotients.json').read_text());bats=[[(0,2)],[(0,1)],[(0,1),(2,18)]];lines=[];records=[]
for r,qq in zip(R,bats):
 for rev in [False,True]:
  Q=[31^(1<<i)^v for i,v in enumerate(r['rows'])] if rev else r['rows'];T=old(Q);queries=[(phys(a),phys(c)) for a,c in qq];lines.extend([f'20 {len(queries)}',' '.join(map(str,T))]+[f'{a} {c}' for a,c in queries]);records.append({'quotient_rows':Q,'reverse':rev,'old_rows':T,'block_queries':qq,'physical_queries':queries})
with tempfile.TemporaryDirectory() as td:
 exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);vals=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
at=0
for r in records:n=len(r['block_queries']);r['values_mod16']=vals[at:at+n];at+=n
print(json.dumps({'physical_count':at,'physical_order':22,'records':records},sort_keys=True,indent=2))
