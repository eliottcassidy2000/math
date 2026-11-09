import itertools,json
from collections import defaultdict
reports=[]

def check(ok,msg):
 if not ok: raise RuntimeError(msg)
def state(B,C):
 n=len(B); full=(1<<n)-1
 h=[0]*(1<<n); h[0]=1
 starts=[[0]*n for _ in h]; ends=[[0]*n for _ in h]
 # Literal path permutation enumeration, independent of producer subset DP.
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(B[x]>>y&1 for x,y in zip(p,p[1:])):
    U=sum(1<<x for x in p);h[U]+=1;starts[U][p[0]]+=1;ends[U][p[-1]]+=1
 s=[sum(r[v] for v in range(n) if C>>v&1) for r in starts]
 e=[sum(r[v] for v in range(n) if C>>v&1) for r in ends]
 return (h[full],sum(h[u]*h[full^u] for u in range(full+1)),sum(h[u]*s[full^u] for u in range(full+1)),sum(e[u]*h[full^u] for u in range(full+1)),sum(e[u]*s[full^u] for u in range(full+1)),s[full],e[full])
def hp(B):
 n=len(B);D=[{} for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for U,row in enumerate(D):
  for v,c in row.items():
   for w in range(n):
    if not U>>w&1 and B[v]>>w&1:
     r=D[U|1<<w];r[w]=r.get(w,0)+c
 return sum(D[-1].values())
def graph(B,C,a,p,b,r):
 n=a+len(B)+b;out=[0]*(n+2);cut=set()
 for i in range(n):
  for j in range(i+1,n):
   edge=(B[i-a]>>(j-a)&1) if a<=i<j<a+len(B) else 1
   out[i if edge else j]|=1<<(j if edge else i)
  out[i]|=1<<n
 for i in range(a):
  if p>>(a-1-i)&1:cut.add(i)
 cut|={a+i for i in range(len(B)) if C>>i&1}
 for i in range(b):
  if r>>(b-1-i)&1:cut.add(a+len(B)+i)
 out[n]|=1<<(n+1)
 for i in range(n):out[n+1 if i in cut else i]|=1<<(i if i in cut else n+1)
 return out

def strong(B):
 n=len(B)
 for v in range(n):
  R=1<<v
  while True:
   S=R
   for w in range(n):
    if R>>w&1:S|=B[w]
   if S==R:break
   R=S
  if R!=(1<<n)-1:return False
 return True

def f(x):
 bits=[i for i in range(x.bit_length()) if x>>i&1]
 return sum(2**(j-i-1)+(j==i+1) for i in bits for j in bits if i<j)
def formula(t,p,b,r):
 h,q,A,G,T,S,E=t
 rev=sum(((r>>i)&1)<<(b-1-i) for i in range(b));first=(r>>(b-1)&1) if b else 0
 return (q*p+A)*(2**(b+1)-rev)+h*(2*r+1-f(r))-p*G-h*f(p)-T-(p%2)*S-E*first
cases=[([24,17,3,6,12],3,(15,120,42,42,10,6,6)),([24,17,3,6,12],5,(15,120,42,42,8,6,6)),([48,49,3,7,12,28],9,(45,432,129,129,20,15,15)),([48,49,3,7,12,28],3,(45,432,129,129,22,15,15)),([40,37,17,6,11,28],3,(45,432,129,129,24,15,15))]
count=0
for B,C,t in cases:
 check(state(B,C)==t,'state');check(strong(B),'core strong')
 for total in range(4):
  for a in range(total+1):
   b=total-a
   for p in range(1<<a):
    for r in range(1<<b):
     g=graph(B,C,a,p,b,r);check(hp(g)==formula(t,p,b,r),'formula')
     check(strong(g)==(not a or p>=1<<(a-1)),'strongness iff for these strong cores');count+=1
reports.append({'status':'PASS','states':[list(x[2]) for x in cases],'literal_graphs':count,'universe':'five states, all prefix/suffix words of total length <=3, including illegal prefixes as hostile controls'})

def compose(u,v):
 h,q,A,G,T,S,E=u;k,l,B,J,U,V,W=v
 return(h*k,q*l,A*l+h*B,J*q+k*G,A*J+T*k+h*U+E*V,S*k,h*W)
H=[40,1,3,6,15,30];ss=[state(H,m) for m in range(64)];fib=defaultdict(list)
for m,s in enumerate(ss):fib[s[2]].append(m)
diffs=[]
for A,ms in sorted(fib.items()):
 if len(ms)==2:
  x,y=ms;d=tuple(v-u for u,v in zip(ss[x],ss[y]));diffs.append([ms,list(d)])
  for p in range(8):
   for b in range(4):
    for r in range(1<<b):
     e=(r>>(b-1)) if b else 0
     check(formula(ss[y],p,b,r)-formula(ss[x],p,b,r)==-p*d[3]-d[4]-(p%2)*d[5]-e*d[6],'H19 difference')
check(hp(graph(H,19,0,0,1,1))==345 and hp(graph(H,42,0,0,1,1))==343,'H19 twins')
C=[16,1,3,7,14];cs=[state(C,m) for m in (2,28,4,0)]
u=compose(cs[0],cs[1]);v=compose(cs[2],cs[3]);delta=tuple(b-a for a,b in zip(u,v));check(delta==(0,0,0,-3021,-558,9,-36),'carry')
B=cases[0][0];BB=[x|sum(1<<j for j in range(5,10)) for x in B]+[x<<5 for x in B]
counts=[]
for c,d in itertools.product((3,5),repeat=2):
 g=graph(BB,c|(d<<5),0,0,0,0);check(strong(g),'repeated strong');counts.append(hp(g))
check(sorted(set(counts))==list(range(min(counts),max(counts)+1,30)) and len(set(counts))==3,'repeated spacing')
reports.append({'H19_fibers':diffs,'carry_difference':delta,'two_switch_counts':counts})
Z=[26,44,57,16,34,9];zs=[state(Z,m) for m in (9,10,17,18)]
for p in range(5):
 for b in range(4):
  for r in range(1<<b):
   e=(r>>(b-1)) if b else 0
   vals=[formula(z,p,b,r)-formula(zs[0],p,b,r) for z in zs]
   check(vals==[0,12*p+6+2*e,36*p+10+6*e,48*p+8+8*e],'H29')
check(sorted(formula(z,0,0,0) for z in zs[1:])==[175,177,179],'H29 triple')
def all_states(B):
 n=len(B);F=(1<<n)-1;st=[[0]*n for _ in range(F+1)];en=[[0]*n for _ in range(F+1)];h=[0]*(F+1);h[0]=1
 for U in range(1,F+1):
  for v in range(n):
   if U>>v&1:
    W=U^(1<<v)
    if not W:st[U][v]=en[U][v]=1
    else:
     st[U][v]=sum(st[W][w] for w in range(n) if W>>w&1 and B[v]>>w&1)
     en[U][v]=sum(en[W][w] for w in range(n) if W>>w&1 and B[w]>>v&1)
  h[U]=sum(en[U])
 q=sum(h[u]*h[F^u] for u in range(F+1));result=[]
 for C in range(F+1):
  s=[sum(r[v] for v in range(n) if C>>v&1) for r in st];e=[sum(r[v] for v in range(n) if C>>v&1) for r in en]
  result.append((h[F],q,sum(h[u]*s[F^u] for u in range(F+1)),sum(e[u]*h[F^u] for u in range(F+1)),sum(e[u]*s[F^u] for u in range(F+1)),s[F],e[F]))
 return result
cyclic=[]
for n,expected in [(5,[8,10]),(7,[68,82,84]),(9,[1080,1196,1252,1258])]:
 B=[sum(1<<((v+j)%n) for j in range(1,(n+1)//2)) for v in range(n)];ss=all_states(B);fiber=[s for m,s in enumerate(ss) if m.bit_count()==2]
 check(sorted({s[4] for s in fiber})==expected,'cyclic T');check(len({s[:4]+s[5:] for s in fiber})==1,'cyclic interface')
 cyclic.append([n,list(ss[0][:2]),expected])
reports.append({'H29_states':zs,'cyclic_hostiles':cyclic})

print(json.dumps(reports,indent=2))
