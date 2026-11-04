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
