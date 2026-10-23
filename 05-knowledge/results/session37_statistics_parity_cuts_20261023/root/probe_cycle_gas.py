import itertools,json
from fractions import Fraction

def test(ok,msg):
 if not ok:raise RuntimeError(msg)
def make():
 n=9;A=[[0]*n for _ in range(n)]
 for i,j in itertools.combinations(range(n),2):
  if i//3!=j//3: win=i
  else:win=i if (j-i)%3==1 else j
  lose=j if win==i else i;A[win][lose]=1
 A[0][8]=0;A[8][0]=1
 return A

def cycles(A):
 n=len(A);out={}
 for mask in range(1,1<<n):
  vs=[i for i in range(n) if mask>>i&1]
  if len(vs)<3 or len(vs)%2==0:continue
  root=vs[0];num=0
  for rest in itertools.permutations(vs[1:]):
   p=(root,)+rest
   if all(A[p[i]][p[(i+1)%len(p)]] for i in range(len(p))):num+=1
  if num:out[mask]=num
 return out

def main():
 A=make();C=cycles(A);full=(1<<9)-1
 # cycle packings recurrence at lowest vertex: leave uncovered or include a cycle containing it.
 memo={0:(1,)}
 def poly(mask):
  if mask in memo:return memo[mask]
  v=mask&-mask;p=list(poly(mask^v))
  for cm,w in C.items():
   if cm&v and cm&mask==cm:
    q=poly(mask^cm)
    if len(p)<len(q)+1:p +=[0]*(len(q)+1-len(p))
    for k,a in enumerate(q):p[k+1]+=w*a
  memo[mask]=tuple(p);return tuple(p)
 p=poly(full);test(p==(1,49,7,1),str(p))
 n=len(A); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not mask>>j&1 and A[i][j]:dp[mask|1<<j][j]+=dp[mask][i]
 test(sum(dp[-1])==135,'independent Hamiltonian count')
 disc=18*7*49-4*7**3+7**2*49**2-4*49**3-27
 test(disc<0,'not complex')
 print(json.dumps(dict(coefficients=p,discriminant=disc,H_at2=sum(a*2**k for k,a in enumerate(p)),outmasks=[sum(A[i][j]<<j for j in range(9)) for i in range(9)],odd_cycle_counts={k:sum(c for m,c in C.items() if m.bit_count()==k) for k in (3,5,7,9)}),indent=2))
if __name__=='__main__':main()
