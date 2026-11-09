import json,itertools

def gate(c,s):
 if not c:raise RuntimeError(s)
def graph(n,code):
 a=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return a
def dp(a):
 n=len(a);z=1<<n;e=[[0]*n for _ in range(z)];s=[[0]*n for _ in range(z)];h=[1]+[0]*(z-1)
 for mask in range(1,z):
  vs=[v for v in range(n) if mask>>v&1]
  for v in vs:
   t=mask^(1<<v)
   if t:
    e[mask][v]=sum(e[t][u] for u in vs if a[u]>>v&1)
    s[mask][v]=sum(s[t][u] for u in vs if a[v]>>u&1)
   else:e[mask][v]=s[mask][v]=1
  h[mask]=sum(e[mask])
 return h,e,s
def four(a,r):
 h,e,s=dp(a);f=(1<<len(a))-1
 return h[f],e[f][r],sum(s[t][r]*h[f^t] for t in range(f+1)),h[f^(1<<r)]
def hq(a):
 h,_,_=dp(a);f=len(h)-1
 return h[f],sum(h[t]*h[f^t] for t in range(f+1))
def ordinal(*blocks):
 n=sum(map(len,blocks));a=[];off=0
 for b in blocks:
  a.extend([(v<<off)|(((1<<n)-1)^((1<<(off+len(b)))-1)) for v in b]);off+=len(b)
 return a
def close(a,c):
 n=len(a)
 return [v|(1<<n)|((1<<(n+1)) if not c>>i&1 else 0) for i,v in enumerate(a)]+[1<<(n+1),c]
def graft(a,r,b):
 n=len(a);extra=((1<<len(b))-1)<<n
 return [v|(extra if i!=r else 0) for i,v in enumerate(a)]+[(v<<n)|(1<<r) for v in b]
def f(bits):
 return sum(2**(j-i-1)+(j==i+1) for i in range(len(bits)) for j in range(i+1,len(bits)) if bits[i] and bits[j])
def scalar(P,R,b):
 pb=list(map(int,bin(P)[2:])) if P else [];rb=[R>>(b-1-i)&1 for i in range(b)]
 rev=sum(v<<i for i,v in enumerate(rb))
 return (11520*P+7561)*(2**(b+1)-rev)+45*(2*R+1-f(rb))-8113*P-45*f(pb)-5258-45*(P%2)-45*(rb[0] if rb else 0)
C=[4,1,18,7,75,31,175,63];L=ordinal([0],[0],C,[0]);cut=1|(124<<2)|(1<<10)
gate(four(close(L,cut),12)==(9909,45,15167,45),'base')
contexts=0
for m in range(3):
 for a in range(m+1):
  b=m-a
  for P in (range(2**(a-1),2**a) if a else [0]):
   for R in range(2**b):
    pref=graph(a,2**(a*(a-1)//2)-1);suff=graph(b,2**(b*(b-1)//2)-1)
    pc=sum((P>>(a-1-i)&1)<<i for i in range(a));rc=sum((R>>(b-1-i)&1)<<i for i in range(b))
    whole=close(ordinal(pref,L,suff),pc|(cut<<a)|(rc<<(a+11)))
    D=four(whole,len(whole)-1)[0];gate(D==scalar(P,R,b),('context',P,R,b,D));gate(D>=9909 and D>=5760*2**m+45,'bounds');contexts+=1
# Exhaust all rooted source tournaments through order four on the claimed invariant slice.
# This tests root ports outside the special family as well.
grafts=0;eligible=0
for n in range(1,5):
 for code in range(2**(n*(n-1)//2)):
  a=graph(n,code)
  for r in range(n):
   D,E,U,d=four(a,r)
   if E!=d:continue
   eligible+=1;V=U-d
   for b in ([0],[2,0],[2,4,1]):
    H,Q=hq(b);got=four(graft(a,r,b),r)
    gate(got==(H*D+(Q-H)*V,H*d,Q*V+H*d,H*d),('graft',a,r,b,got));grafts+=1
# Exact endpoint and inverse-log inequalities over a broad integer range.
for k in range(101):
 R=4096**k;num=245051264*R+1602510*k+1118023
 gate(num%24843==0,'integral D');D=num//24843
 gate(24843*(D-2*k)>=245051264*R,'inverse bound')
print(json.dumps({'status':'PASS','literal_binary_contexts':contexts,'rooted_invariant_sources':eligible,'literal_grafts':grafts,'exact_depth_bounds':101,'caveat':'Unit-gap necessity only nondegenerate k>=1; O(N log N) counts candidates, not bit operations.'},indent=2))
