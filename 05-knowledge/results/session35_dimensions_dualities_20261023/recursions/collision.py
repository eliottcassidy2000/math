import sympy as s
from functools import lru_cache
z=s.symbols('z')
rows=[[58,4,9,18,38,14],[22,12,8,17,38,15]]
for out in rows:
 n=len(out);M=s.zeros(n)
 for i in range(n):
  for j in range(n):
   if i!=j:M[i,j]=1 if out[i]>>j&1 else -1
 if M!=-M.T:raise RuntimeError('not tournament')
 @lru_cache(None)
 def hp(mask,last):
  if mask==(1<<last):return 1
  return sum(hp(mask^(1<<last),j) for j in range(n) if j!=last and mask>>j&1 and out[j]>>last&1)
 H=sum(hp((1<<n)-1,j) for j in range(n))
 B=z*s.eye(n)-M;p=s.factor(B.det());q=s.factor((s.ones(1,n)*B.adjugate()*s.ones(n,1))[0])
 print('out',out,'H',H,'p',p,'q',q)
def path_coeff(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    for k in range(n):
     if not(mask>>k&1) and out[j]>>k&1:dp[mask|1<<k][k]+=dp[mask][j]
 return [1]+[sum(dp[m]) for m in range(1,1<<n)]
a,b=map(path_coeff,rows)
for k in range(1,7):
 diffs=[(m,a[m],b[m]) for m in range(1,64) if m.bit_count()==k and a[m]!=b[m]]
 print('degree',k,'differences',diffs[:5], 'count',len(diffs))
# E factorization uses new labels 0, n+1, original labels shifted one
for out in rows:
 n=len(out);e=[((1<<(n+2))-1)^1]+[(v<<1)|(1<<(n+1)) for v in out]+[0]
 c=path_coeff(out);ce=path_coeff(e)
 for mask in range(1<<(n+2)):
  inner=(mask>>1)&((1<<n)-1)
  if ce[mask]!=c[inner]:raise RuntimeError('E squarefree polynomial identity failed')
print('E squarefree factorization checked on every subset for both n6 witnesses.')
import itertools
for out in rows:
 cyc={}
 for k in [3,5]:
  sets=[]
  for sub in itertools.combinations(range(6),k):
   a=min(sub)
   for rest in itertools.permutations([i for i in sub if i!=a]):
    c=(a,)+rest
    if all(out[c[i]]>>c[(i+1)%k]&1 for i in range(k)):sets.append(sum(1<<i for i in c))
  cyc[k]=sets
 pairs=sum(not(a&b) for i,a in enumerate(cyc[3]) for b in cyc[3][i+1:])
 print('OCF counts c3,c5,disjoint_triangle_pairs,H:',len(cyc[3]),len(cyc[5]),pairs,1+2*(len(cyc[3])+len(cyc[5]))+4*pairs)
