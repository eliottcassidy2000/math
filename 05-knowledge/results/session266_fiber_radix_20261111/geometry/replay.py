"""New electrical core-padding constructors, exact rational arithmetic."""
import itertools,json
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).parent

def join(A,B):
 a,b=len(A),len(B)
 return tuple(x|(((1<<b)-1)<<a) for x in A)+tuple(x<<a for x in B)
def path_data(T):
 n=len(T); full=(1<<n)-1
 pred=[sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n)]
 E=[[0]*n for _ in range(1<<n)];S=[[0]*n for _ in range(1<<n)]
 for v in range(n): E[1<<v][v]=S[1<<v][v]=1
 for mask in range(1,1<<n):
  if mask&(mask-1)==0:continue
  for v in range(n):
   if not mask>>v&1:continue
   rem=mask^(1<<v);x=rem&pred[v];y=rem&T[v]
   while x:
    bit=x&-x;x-=bit;E[mask][v]+=E[rem][bit.bit_length()-1]
   while y:
    bit=y&-y;y-=bit;S[mask][v]+=S[rem][bit.bit_length()-1]
 W=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  rem=full^mask
  for u in range(n):
   if not E[mask][u]:continue
   x=T[u]&rem
   while x:
    bit=x&-x;x-=bit;v=bit.bit_length()-1
    p=E[mask][u]*S[rem][v];W[u][v]+=p;W[v][u]+=p
 return sum(E[full]),W,S[full],E[full]
def solve(A,b):
 a=[[F(x) for x in row]+[F(y)] for row,y in zip(A,b)];n=len(a)
 for i in range(n):
  j=next(j for j in range(i,n) if a[j][i]);a[i],a[j]=a[j],a[i]
  q=a[i][i];a[i]=[x/q for x in a[i]]
  for j in range(n):
   if j!=i:
    q=a[j][i];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
 return [row[-1] for row in a]
def resistance(W,s,t):
 vs=[v for v in range(len(W)) if v!=t]
 L=[[sum(W[u]) if u==v else -W[u][v] for v in vs] for u in vs]
 x=solve(L,[int(u==s) for u in vs]);return x[vs.index(s)]
def atom(T):
 H,W,start,end=path_data(T);X=join(join((0,),T),(0,));HX,WX,_,_=path_data(X)
 rho=HX*resistance(WX,0,len(X)-1)
 if HX!=H:raise RuntimeError('observer count')
 if not 2<=rho<=len(T)+1:raise RuntimeError('resistance bounds')
 if (rho==2)!=(start==end):raise RuntimeError('balanced equality')
 return dict(H=H,n=len(T),rho=str(rho),defect=str(len(T)+1-rho),start=start,end=end)
def main():
 C3=(2,4,1);C4=(2,4,9,3)
 R5=tuple(sum(1<<((u+j)%5) for j in (1,2)) for u in range(5))
 P7=tuple(sum(1<<((u+j)%7) for j in (1,2,4)) for u in range(7))
 cases={'C3':C3,'C4':C4,'R5':R5,'P7':P7,'transitive3':(6,4,0)}
 rows={k:atom(T) for k,T in cases.items()};checks=3*len(cases)
 composed=[]
 for names,gaps in [(('C3',),(1,1)),(('C3',),(3,2)),(('C4',),(1,3)),(('C3','C3'),(2,1,2)),(('C3','C4'),(1,2,1)),(('C4','C3'),(2,1,1))]:
  T=();k=len(names)
  for i in range(k+1):
   for j in range(gaps[i]):T=join(T,(0,))
   if i<k:T=join(T,cases[names[i]])
  H,W,_,_=path_data(T);scaled=H*resistance(W,0,len(T)-1)
  D=sum((F(rows[name]['defect']) for name in names),F(0))
  if scaled!=len(T)-1-D:raise RuntimeError(('series law',names,gaps,scaled,D))
  composed.append(dict(cores=names,gaps=gaps,n=len(T),H=H,H_times_R=str(scaled),defect=str(D)));checks+=1
 # Independent path enumeration on the atom networks only.
 for T in cases.values():
  H,W,_,_=path_data(T);paths=[p for p in itertools.permutations(range(len(T))) if all(T[u]>>v&1 for u,v in zip(p,p[1:]))]
  V=[[0]*len(T) for _ in T]
  for p in paths:
   for u,v in zip(p,p[1:]):V[u][v]+=1;V[v][u]+=1
  if len(paths)!=H or V!=W:raise RuntimeError('path check')
  checks+=1
 out=dict(status='PASS',active_checks=checks,atoms=rows,composites=composed)
 (ROOT/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
