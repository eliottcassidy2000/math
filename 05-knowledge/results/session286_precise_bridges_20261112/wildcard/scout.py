import sys,json,itertools,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session285_five_web_20261112/tournaments'))
from explore import coefficients

def cyclic(n,S): return [sum(1<<((i+j)%n) for j in S) for i in range(n)]
def corr(mask,n):
 C=[i for i in range(n) if mask>>i&1]
 return tuple(sum((mask>>((i+d)%n))&1 for i in C) for d in range(1,n//2+1))
def canon(m,n):
 return min(sum(1<<((s*i+t)%n) for i in range(n) if m>>i&1) for s in [1,-1] for t in range(n))

def run(n,S):
 r=cyclic(n,S);h,Q,u,v,a,b,B,s,t=coefficients(r)
 W=[[B[i][j]+B[j][i] for j in range(n)]for i in range(n)]
 bins={};best=None;twins=0
 for k in range(3,n//2+1):
  for C in itertools.combinations(range(n),k):
   m=sum(1<<i for i in C)
   if canon(m,n)!=m:continue
   bins.setdefault((k,corr(m,n)),[]).append(m)
 for (k,cc),ms in bins.items():
  if len(ms)<2:continue
  twins+=1
  for a,b in itertools.combinations(ms,2):
   for sh in range(n):
    bb=sum(1<<((j+sh)%n)for j in range(n)if b>>j&1)
    for v in range(n):
     if (a>>v&1)!=(bb>>v&1):continue
     diff=sum(W[v][j]*((a>>j&1)-(bb>>j&1))for j in range(n))
     if not diff:continue
     if best is None or abs(diff)<best[0]:best=(abs(diff),k,a,bb,v,diff)
     if abs(diff)==2:return dict(n=n,S=S,h=h,Q=Q,weights=W[0],best=best,twins=twins)
 return dict(n=n,S=S,h=h,Q=Q,weights=W[0],best=best,twins=twins)
if __name__=='__main__':
 for n,S in [(9,[1,2,3,4]),(9,[1,7,6,5]),(11,[1,2,3,4,5]),(13,list(range(1,7))),(13,[1,2,3,4,5,7]),(15,list(range(1,8)))]:
  print(json.dumps(run(n,S)),flush=True)
