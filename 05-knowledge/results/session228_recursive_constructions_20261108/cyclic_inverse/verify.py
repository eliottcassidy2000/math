"""Exact standalone verifier. No third-party packages; checks survive -O."""
import itertools,json,math

def check(x,msg):
 if not x: raise RuntimeError(msg)
def tournaments(n):
 edges=list(itertools.combinations(range(n),2))
 for code in range(1<<len(edges)):
  g=[0]*n
  for k,(i,j) in enumerate(edges):
   if code>>k&1:g[i]|=1<<j
   else:g[j]|=1<<i
  yield tuple(g)
def hsub(g):
 n=len(g); N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for v in range(n):dp[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:dp[S][v]=sum(dp[T][w] for w in range(n) if T>>w&1 and g[w]>>v&1)
  h[S]=sum(dp[S])
 return h,dp

def qvals(g,m=3):
 h,_=hsub(g);full=(1<<len(g))-1;v=[0]*(full+1);v[0]=1;out=[v[full]]
 for _ in range(m):
  w=[]
  for S in range(full+1):
   c=0;T=S
   while True:
    c+=v[S^T]*h[T]
    if not T:break
    T=(T-1)&S
   w.append(c)
  v=w;out.append(v[full])
 return out

def cyc(A,B):
 a=len(A);b=len(B);n=a+b+1;u=n-1;g=[0]*n
 for i,row in enumerate(A):g[i]=row|(((1<<b)-1)<<a)
 for i,row in enumerate(B):g[a+i]=(row<<a)|(1<<u)
 g[u]=(1<<a)-1
 return tuple(g)
def state(g,u):
 h,dp=hsub(g);full=(1<<len(g))-1;rest=full^(1<<u);rev=tuple(sum(1<<j for j in range(len(g)) if g[j]>>i&1) for i in range(len(g)));_,end=hsub(rev)
 U=0;S=rest
 while True:
  U+=h[S]*end[full^S][u]
  if not S:break
  S=(S-1)&rest
 return (h[full],dp[full][u],U,h[rest])
def pair(q):return(q[1],q[2]-2*q[1])
def formula(p,q):
 a,r=p;b,s=q;d=a*b;D=3*d+a*s+b*r+r*s
 return(D,d,D+b*r,d)
def inverse(x,bank):
 D,E,U,d=x
 if D<=0 or d<=0 or E!=d:return []
 t=U-D;z=2*D-U-3*d
 if t<0 or z<0:return []
 out=[]
 for b in range(1,math.isqrt(math.gcd(d,t))+1):
  if math.gcd(d,t)%b:continue
  for bb in sorted({b,math.gcd(d,t)//b}):
   a=d//bb;r=t//bb
   if bb*z%(d+t):continue
   s=bb*z//(d+t)
   if (a,r) in bank and (bb,s) in bank:out.append(((a,r),(bb,s)))
 return sorted(out)
def main():
 small=[];bank={};counts=[]
 for n in range(1,7):
  c=0
  for g in tournaments(n):
   h,_=hsub(g);full=(1<<n)-1;Q=sum(h[S]*h[full^S] for S in range(full+1));p=(h[full],Q-2*h[full]);bank.setdefault(p,g);c+=1
   if n<=3:small.append(g)
  counts.append(c)
 literal=0
 for A in small:
  for B in small:
   g=cyc(A,B);qa=qvals(A,4);qb=qvals(B,4);qc=qvals(g,3)
   for t in range(1,4):
    pred=t*(qa[t]*qb[t]+(qa[t+1]-qa[t])*(qb[t+1]-qb[t])+(qa[t]-qa[t-1])*(qb[t]-qb[t-1]))
    check(qc[t]==pred,'Q identity')
   check(state(g,len(g)-1)==formula(pair(qa),pair(qb)),'root state');literal+=1
 # Complete count image <=128: all possible factors have order<=6 by proof.
 smallbank={p:g for p,g in bank.items() if p[0]+p[1]<=126}
 scalar={};states={}
 for p in smallbank:
  for q in smallbank:
   x=formula(p,q)
   if x[0]<=128:scalar.setdefault(x[0],(p,q));states.setdefault(x,set()).add((p,q))
 for x,facts in states.items():
  check(set(inverse(x,smallbank))==facts,'complete inverse')
  p,q=sorted(facts)[0];g=cyc(smallbank[p],smallbank[q])
  check(state(g,len(g)-1)==x,'literal inverse witness')
 # Every vector in this box, physical or not: inverse iff in direct finite image.
 invchecks=0
 for D in range(1,130,2):
  for d in range(1,44,2):
   for U in range(D,2*D+1):
    x=(D,d,U,d);got=inverse(x,smallbank)
    check(bool(got)==any(formula(p,q)==x for p,q in got),'inverse consistency')
    if D<=128:check(set(got)==states.get(x,set()),'inverse hostile box')
    invchecks+=1
 check(13 not in scalar,'13 must fail')
 hostile=(20,1,2,7,14);hx=state(hostile,1)
 check(hx==(9,3,9,3) and not inverse(hx,bank),'physical inverse hostile')
 # Same scalar source H=1 does not determine response: TT1 vs TT2 against TT1.
 check(formula((1,0),(1,0))[0]==3 and formula((1,2),(1,0))[0]==5,'sidecar hostile')
 print(json.dumps({'status':'PASS','tournament_counts_by_order':counts,'distinct_pairs_through_6':len(bank),'literal_binary_substitutions_through_factor_order_3':literal,'polynomial_evaluations_per_literal':3,'rooted_inverse_states_H_le_128':len(states),'inverse_hostile_vectors':invchecks,'complete_binary_cyclic_counts_le_128':sorted(scalar),'missing_odds_le_127':[n for n in range(1,128,2) if n not in scalar],'first_missing_odd_other_than_1_7_21':13,'physical_state_outside_binary_image':hx,'example_11':scalar[11],'example_33':scalar[33]},indent=2,sort_keys=True))
if __name__=='__main__':main()
