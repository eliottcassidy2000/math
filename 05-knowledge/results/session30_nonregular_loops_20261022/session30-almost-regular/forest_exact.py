"""Independent signed-path/subset-partition forest algorithm; integer arithmetic."""
import json,math,pathlib
from fractions import Fraction
from experiment import frombits
ROOT=pathlib.Path(__file__).resolve().parent

def forest(S):
 n=len(S);size=1<<n;dp=[[0]*n for _ in range(size)];path=[0]*size
 for m in range(1,size):
  if m&(m-1)==0:dp[m][m.bit_length()-1]=1
  else:
   for j in range(n):
    if m>>j&1:
     p=m^(1<<j);dp[m][j]=sum(dp[p][k]*int(S[k,j]) for k in range(n) if p>>k&1)
  path[m]=sum(dp[m])
  if m.bit_count()%2==0:assert path[m]==0
 F=[[0]*n for _ in range(size)];F[0][0]=1
 for m in range(1,size):
  first=m&-m;rest=m^first;sub=rest
  while True:
   block=sub|first;b=block.bit_count();w=path[block]
   if w:
    for ell,val in enumerate(F[m^block]):
     if val:F[m][ell+b-1]+=w*val
   if sub==0:break
   sub=(sub-1)&rest
 return F[-1]

def logcoeff(phi):
 out=[Fraction(0)]*len(phi)
 for k in range(1,len(phi)):
  out[k]=phi[k]-sum(Fraction(j,k)*out[j]*phi[k-j] for j in range(1,k))
 return out

def main():
 rows=[json.loads(x) for x in open(ROOT/'data.jsonl')];chosen=[r for r in rows if not r['name'].startswith('exhaustive') and r['n']<=12]
 classes=json.load(open(ROOT/'small_classes.json'))
 for n,group in classes.items():
  for i,r in enumerate(group):chosen.append(dict(r,n=int(n),name=f'exact-small-{n}-{i}'))
 with open(ROOT/'forest_controls.jsonl','w') as f:
  for r in chosen:
   n=r['n'];S=frombits(n,r['bits']);N=forest(S);val=sum(Fraction(v,math.factorial(n)//math.factorial(n-L)) for L,v in enumerate(N));assert val==Fraction(int(r['H'])*2**(n-1),math.factorial(n));assert N[2] if n>2 else True
   if n>2:assert N[2]==n*(n-1)-sum(int(s)**2 for s in S.sum(axis=1))
   phi=[Fraction(v,n**L) for L,v in enumerate(N)];lc=logcoeff(phi)
   out=dict(name=r['name'],n=n,bits=r['bits'],H=r['H'],N=N,log_phi=[str(x) for x in lc]);f.write(json.dumps(out)+'\n');print(r['name'],N,flush=True)
if __name__=='__main__':main()
