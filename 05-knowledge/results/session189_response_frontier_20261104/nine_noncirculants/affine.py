from pathlib import Path
from itertools import product
import random,subprocess,json
D=Path(__file__).resolve().parent

def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def aug(Q,a,c):
 n=len(Q);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(Q)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def solve(equations,n=18):
 piv={};inconsistent=False
 for row,rhs in equations:
  while row:
   k=row.bit_length()-1
   if k not in piv:piv[k]=(row,rhs);break
   a,b=piv[k];row^=a;rhs^=b
  if not row and rhs:inconsistent=True
 if inconsistent:return dict(rank=len(piv),inconsistent=True,solutions=[])
 free=[i for i in range(n) if i not in piv];solutions=[]
 for mask in range(1<<len(free)):
  x=sum(((mask>>i)&1)<<j for i,j in enumerate(free))
  for k,(row,rhs) in sorted(piv.items()):
   bit=rhs^((row&x).bit_count()%2);x|=bit<<k
  solutions.append(x)
 return dict(rank=len(piv),inconsistent=False,solutions=solutions)
def main():
 Qs=[[192,289,99,71,15,153,306,94,189],[320,193,83,71,43,15,432,61,190],[448,209,195,71,45,15,432,312,62],[384,177,147,135,105,77,15,368,126]];masks=[m for m in range(512) if m.bit_count()<=2];rng=random.Random(1899);queries=[(0,7),(8,264),(1,2),(1,3),(3,1),(1,256),(256,1)]+rng.sample(list(product(masks,repeat=2)),200);queries=list(dict.fromkeys(queries));cases=[aug(Q if not side else op(Q),a,c) for Q in Qs for a,c in queries for side in range(2)]
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'runs.cpp'),'-o',str(D/'runs')],check=True);inp=str(len(cases))+'\n'+''.join('9\n'+' '.join(map(str,R))+'\n' for R in cases);raw=subprocess.check_output([str(D/'runs')],input=inp.encode()).decode();(D/'run_values.txt').write_text(raw);rows=[list(map(int,l.split())) for l in raw.splitlines()];records=[];at=0
 for Q in Qs:
  bank=[];equations=[];safe_equations=[]
  for a,c in queries:
   S,U=rows[at:at+2];at+=2;d=[s-u for s,u in zip(S,U)];baseline=(d[0]+2*sum(d[1:10])+4*sum(d[19:55]))%8
   if baseline%4:raise RuntimeError('nonaffine baseline residue')
   if any(x%2 for x in d[10:19]):raise RuntimeError('odd triple-visit contrast')
   bits=[x%2 for x in d[1:10]]+[(x//2)%2 for x in d[10:19]];row=sum(x<<i for i,x in enumerate(bits));rhs=baseline//4;safe=(d[55]%2==0);eq=(row,rhs);equations.append(eq)
   if safe:safe_equations.append(eq)
   bank.append(dict(query=[a,c],source=S,target=U,baseline=baseline,linear_bits=bits,cycle_safe=safe))
  records.append(dict(rows=Q,queries=bank,all_queries=solve(equations),safe_queries=solve(safe_equations),safe_query_count=len(safe_equations)))
 (D/'analysis.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps([{k:v for k,v in r.items() if k!='queries'} for r in records],indent=2))
if __name__=='__main__':main()
