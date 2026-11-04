import json,sys,importlib.util
from pathlib import Path
import core as q
def extend(T,c):
 n=len(T);return tuple(r|(((c>>i)&1)<<n) for i,r in enumerate(T))+(((1<<n)-1)^c,)

def main():
 rows=json.load(open(Path(__file__).parent/'quotients.json'));out=[]
 for case in rows:
  if any(case['p']):continue
  T=case['rows'];U=q.converse(T);n=len(T)
  Rs=[q.repeated_endpoint(T,j) for j in range(n)];Ru=[q.repeated_endpoint(U,j) for j in range(n)]
  starts=[[Ru[j][i]-Rs[j][i] for j in range(n)] for i in range(n)]
  d0=[];d2=[]
  for i in range(n):
   t=extend(T,1<<i);u=extend(U,1<<i)
   d0.append(sum(q.endpoint(t))-sum(q.endpoint(u)))
   d2.append([sum(q.repeated_endpoint(t,j))-sum(q.repeated_endpoint(u,j)) for j in range(n)])
  patterns=[]
  for V in range(1<<n):
   d=[(d0[i]+2*sum(d2[i][j] for j in range(n) if V>>j&1))%4 for i in range(n)]
   sig=[(case['S'][i]-case['E'][i]+2*sum(starts[i][j] for j in range(n) if V>>j&1))%4 for i in range(n)]
   a=[(2*sig[i]+2*d[i])%8 for i in range(n)]
   patterns.append(dict(V=V,first=d,second=a))
  blind=[r['V'] for r in patterns if not any(r['first']) and not any(r['second'])]
  out.append(dict(rows=T,start_run_skews=starts,one_port_base=d0,one_port_double=d2,patterns=patterns,blind=blind))
 json.dump(out,open(Path(__file__).parent/'anchor_exploration.json','w'),indent=2)
 print('38 cases; blind patterns distribution',[(r['rows'],len(r['blind']),r['blind']) for r in out if r['blind']])
if __name__=='__main__':main()
