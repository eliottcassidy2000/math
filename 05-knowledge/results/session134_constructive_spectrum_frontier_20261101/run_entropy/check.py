"""Finite controls for marked-spine run entropy; standard library only."""
import itertools,json,math
from collections import Counter

def need(x,m):
 if not x:raise RuntimeError(m)
def changes(w):return sum(a!=b for a,b in zip(w,w[1:]))
def words_bound(k,N,R):
 return sum((1<<math.comb(j,2))*(1+(1<<j)*sum(math.comb(N,r+1)*((1<<j)-1)**r for r in range(min(R,N-1)+1))) for j in range(k+1))
def rows(n,c):
 a=[0]*n;z=0
 for i in range(n):
  for j in range(i+1,n):
   if c>>z&1:a[i]|=1<<j
   else:a[j]|=1<<i
   z+=1
 return a
def induced(a,vs):return [sum(((a[v]>>u)&1)<<j for j,u in enumerate(vs)) for v in vs]
def h(a):
 n=len(a)
 if not n:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    prev=mask^(1<<v)
    if prev:dp[mask][v]=sum(dp[prev][u] for u in range(n) if prev>>u&1 and a[u]>>v&1)
 return sum(dp[-1])
def core_vertices(a):
 reach=[row|1<<i for i,row in enumerate(a)]
 for k in range(len(a)):
  for i in range(len(a)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return [i for i in range(len(a)) if any(j!=i and reach[i]>>j&1 and reach[j]>>i&1 for j in range(len(a)))]
def word(a,apices):
 chain=[v for v in range(len(a)) if v not in apices];cm=sum(1<<v for v in chain)
 chain.sort(key=lambda v:(a[v]&cm).bit_count(),reverse=True)
 if not all(a[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:]):return None
 return tuple(sum(((a[x]>>v)&1)<<i for i,x in enumerate(apices)) for v in chain)
def main():
 histchecks=0;projection=0;graphchecks=0;profiles=[]
 for q in (1,2,4,8):
  maxm=6 if q<=4 else 4
  totals=Counter()
  for m in range(1,maxm+1):
   hist=Counter(changes(w) for w in itertools.product(range(q),repeat=m))
   for r in range(m):
    need(hist[r]==q*math.comb(m-1,r)*(q-1)**r,'run histogram');histchecks+=1
   totals.update(hist)
  for R in range(maxm):
   need(sum(totals[r] for r in range(R+1))==q*sum(math.comb(maxm,r+1)*(q-1)**r for r in range(R+1)),'hockey-stick sum');histchecks+=1
 for m in range(7):
  for w in itertools.product(range(4),repeat=m):
   before=changes(w)
   for subset in range(1<<m):
    sub=[w[i] for i in range(m) if subset>>i&1]
    for keep in (0,1,2,3):
     need(changes([v&keep for v in sub])<=before,'deletion/projection changes increased');projection+=1
 for n in range(6):
  for c in range(1<<math.comb(n,2)):
   a=rows(n,c);C=core_vertices(a);b=induced(a,C)
   need(h(a)==h(b),'SCC kernel changed H')
   for k in range(min(2,n)+1):
    for A in itertools.combinations(range(n),k):
     w=word(a,A)
     if w is None:continue
     B=tuple(C.index(v) for v in A if v in C);new=word(b,B)
     need(new is not None and changes(new)<=changes(w),'marked core lost transitive spine or run budget')
     kept=sorted(set(C)|set(A));bb=induced(a,kept);AA=tuple(kept.index(v) for v in A);ww=word(bb,AA)
     need(ww is not None and changes(ww)<=changes(w),'retained-apex compression')
     need(h(bb)==h(a) and h(a)**3>=5**len(ww),'spine-specific Busch bound');graphchecks+=1
 # Fully rational certificates for elementary exponential inequalities used in prose.
 e4=sum((4**j)/math.factorial(j) for j in range(18))
 from fractions import Fraction
 exp4=sum((Fraction(4)**j)/math.factorial(j) for j in range(18))
 exp35=sum((Fraction(7,2)**j)/math.factorial(j) for j in range(15))
 eupper=sum(Fraction(1,math.factorial(j)) for j in range(7))+Fraction(1,math.factorial(7))*Fraction(8,7)
 need(eupper<Fraction(11,4) and exp4>54 and exp35>30,'logarithm certificates')
 for N in (5,10,20,40,80):
  profiles.append({'N':N,'two_apex_bounds':[words_bound(2,N,r) for r in range(min(6,N))]})
 return {'status':'PASS','word_histogram_checks':histchecks,'deletion_projection_checks':projection,
  'marked_tournament_kernel_checks':graphchecks,'finite_bounds':profiles,
  'rational_exp4_lower_bound':str(exp4),'rational_exp3_5_lower_bound':str(exp35),
  'scope':'Finite verification of encoding and monotone run resources; entropy consequences are proved in THEORY.md, not inferred from these checks.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
