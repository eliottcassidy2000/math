import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
import json
from verify import literal,need
root=Path(__file__).parent;w=json.loads((root/'SHORT12_UNIT_LOOP_WITNESS.json').read_text());core=literal(w['adj'])[w['mask']]
def matrix(s):
 h,q,A,G,T,S,E=s
 cats={(0,1):S,(0,2):h-S,(1,0):h-E,(2,0):E,(2,1):T,(1,1):A-S-T,(2,2):G-E-T,(1,2):q-2*h-A-G+S+E+T}
 need(min(cats.values())>=0,'category negative');M=[[0]*6 for _ in range(6)]
 for old in range(6):
  p,r=divmod(old,2)
  for (pp,rr),num in cats.items():
   if not r and rr==2:continue
   new=2*(pp or p)+int(bool(r or rr));M[new][old]+=num
 return M
one=[matrix((1,2,b,b,0,b,b)) for b in (0,1)];mc=matrix(core)
def run(ms):
 v=[1,0,0,0,0,0]
 for M in ms:v=[sum(a*b for a,b in zip(row,v)) for row in M]
 return sum(v[:4])
rows=[]
for k in range(16):
 expected=(245051264*4096**k+1602510*k+1118023)//24843
 need((245051264*4096**k+1602510*k+1118023)%24843==0,'fraction')
 for j in range(k+1):
  mats=[one[1],one[0]]+[one[int(b)] for b in w['X']]*j+[mc]+[one[int(b)] for b in w['Y']]*(k-j)+[one[1],one[0]]
  value=run(mats);need(value==expected-2*j,'six state count')
 rows.append(dict(k=k,first_count=expected,last_count=expected-2*k))
print(json.dumps(dict(status='PASS',method='literal core categories plus independent six-state DP',checks=rows),indent=2))
