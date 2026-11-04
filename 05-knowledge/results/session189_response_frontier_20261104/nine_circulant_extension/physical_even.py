from pathlib import Path
from itertools import product
from math import factorial
import json,subprocess
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def build(Q,sizes,a,c):
 sets=[];T=[]
 for size in sizes:sets.append(list(range(len(T),len(T)+size)));T.extend([0]*size)
 for i,S in enumerate(sets):
  for v in S:T[v]=sum(1<<u for u in S if u>v)|sum(1<<u for j,U in enumerate(sets) if Q[i]>>j&1 for u in U)
 aa=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);cc=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);n=len(T);T=[r|(((aa>>i)&1)<<n)|(((cc>>i)&1)<<(n+1)) for i,r in enumerate(T)]+[(((1<<n)-1)^aa)|(1<<(n+1)),((1<<n)-1)^cc]
 return dict(quotient_rows=Q,child_sizes=sizes,child_sets=sets,augmented_rows=T,physical_query=[aa,cc],port_names=[n,n+1])
def main():
 data=json.loads((D/'even_bank.json').read_text());records=[];cases=[]
 for qi,q in enumerate(data['queries']):
  for side in range(2):
   Q=op(data['rows']) if side else data['rows']
   for si,sizes in enumerate(q['shapes']):
    r=build(Q,sizes,*q['query']);r.update(query_index=qi,side=side,shape_index=si);records.append(r);cases.append(r['augmented_rows'])
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical_h.cpp'),'-o',str(D/'physical_h')],check=True);inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);raw=subprocess.check_output([str(D/'physical_h')],input=inp.encode()).decode();(D/'physical_even_values.txt').write_text(raw);values=list(map(int,raw.split()));need(len(values)==len(cases),'all literal physical controls')
 for r,H in zip(records,values):r['H']=H
 S=[[0]*5 for _ in range(5)];S[0][0]=1
 for n in range(1,5):
  for k in range(1,n+1):S[n][k]=S[n-1][k-1]+k*S[n-1][k]
 for qi,q in enumerate(data['queries']):
  shapes=[tuple(s) for s in q['shapes']]
  for side,key in [(0,'source_counts'),(1,'target_counts')]:
   hs={tuple(r['child_sizes']):r['H'] for r in records if r['query_index']==qi and r['side']==side};recovered={}
   for sizes in sorted(shapes,key=lambda s:(sum(s),s)):
    remaining=hs[sizes]
    for smaller in product(*(range(1,n+1) for n in sizes)):
     if smaller==sizes:continue
     weight=1
     for n,k in zip(sizes,smaller):weight*=factorial(k)*S[n][k]
     remaining-=weight*recovered[smaller]
    diagonal=1
    for n in sizes:diagonal*=factorial(n)
    need(remaining%diagonal==0,'exact multivariate Stirling inversion');recovered[sizes]=remaining//diagonal
   need([recovered[s] for s in shapes]==q[key],'all run counts independently recovered from actual physical graphs')
 (D/'physical_even_records.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',physical_graphs=len(cases),exact_recovered_run_counts=len(cases),maximum_augmented_order=max(map(len,cases)))))
if __name__=='__main__':main()
