import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
import json,itertools,math
from core import compose,states
w=json.loads(Path(__file__).with_name('SHORT12_UNIT_LOOP_WITNESS.json').read_text());core=states(w['adj'])[w['mask']];one=[(1,2,b,b,0,b,b) for b in (0,1)]
def fold(ps):
 s=ps[0]
 for z in ps[1:]:s=compose(s,z)
 return s
def need(x,m):
 if not x:raise RuntimeError(m)
contexts=0;values=0
for a in range(5):
 for b in range(5):
  for P in range(1<<(a-1),1<<a) if a else [0]:
   prefix=[] if not a else [one[int(c)] for c in f'{P:0{a}b}']
   for R in range(1<<b):
    suffix=[] if not b else [one[int(c)] for c in f'{R:0{b}b}'];rev=0 if not b else int(f'{R:0{b}b}'[::-1],2);counts=[]
    for k in range(5):
     vals=[]
     for j in range(k+1):
      parts=prefix+[one[1],one[0]]+[one[int(c)] for c in w['X']]*j+[core]+[one[int(c)] for c in w['Y']]*(k-j)+[one[1]]+suffix
      z=fold(parts);N=2*z[2]+z[0]-z[4];vals.append(N);values+=1
      need(N>=5760*4096**k*2**(a+b),'lower bound')
     need(vals==[vals[0]-2*j for j in range(k+1)],'window')
     need((vals[0]-(2**(b+1)-rev-13*P+7+5*k))%45==0,'residue')
     counts.append(vals[0])
    need(all(counts[k+2]==4097*counts[k+1]-4096*counts[k]-264150 for k in range(3)),'recurrence')
    contexts+=1
print(json.dumps(dict(status='PASS',contexts=contexts,scalar_values=values,k_max=4,prefix_length_max=4,suffix_length_max=4,recurrence=[4097,-4096,-264150]),indent=2))
