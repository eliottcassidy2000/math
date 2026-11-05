import json,itertools
T=(1,0,2,1)
def need(b,msg):
 if not b:raise RuntimeError(msg)
def move(v,a):
 D,E,U,d=v
 return {'I':v,'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[a]
def profile(word):
 v=T
 for a in word:v=move(v,a)
 return v
def birth(g,r,a):
 if a=='I':return g[:],r
 n=len(g);out=g[:]+[0]
 for j in range(n):
  if a!='O' and j==r:out[n]|=1<<j
  else:out[j]|=1<<n
 return out,n if a=='J' else r
def graph(word):
 g=[2,0];r=0
 for a in word:g,r=birth(g,r,a)
 return g,r
def inv(g,r=0):
 n=len(g);N=1<<n
 dp=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)]
 for j in range(n):dp[1<<j][j]=1
 sp[1<<r][r]=1
 for s in range(1,N):
  for j in range(n):
   if not(s>>j&1):continue
   avail=g[j]&~s&(N-1)
   while avail:
    b=avail&-avail;k=b.bit_length()-1;avail-=b
    dp[s|b][k]+=dp[s][j];sp[s|b][k]+=sp[s][j]
 hs=[sum(row) for row in dp];hs[0]=1
 return (hs[-1],dp[-1][r],sum(sum(sp[s])*hs[(N-1)^s] for s in range(1,N)),hs[(N-1)^(1<<r)]),sum(sp[-1]),sum(sp[-1][j] for j in range(n) if g[j]>>r&1)
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
def K(v,w):return sum(v[i]*Q[i][j]*w[j] for i in range(4) for j in range(4))
def NM(v,w,kind):
 D,E,U,d=v;g,s,A,b=w
 if kind=='N':return (D-E+U)*g-(D-E+d)*s+(D+3*U-d)*A-(U+E)*b
 return (D+2*U)*g-(U+d)*s+(2*D-E+6*U-d)*A-(U+E)*b
def c2(a,b,kind):
 g,r=graph(a);g,x=birth(g,r,'J');entries={r,x}
 if kind=='M':g,y=birth(g,x,'J');entries={r,y}
 h,q=graph(b);m=len(h);h=[((1<<m)-1)^(1<<j)^h[j] for j in range(m)]
 n=len(g);out=g[:]+[row<<n for row in h]
 for i in range(n):
  for j in range(m):
   if j==q and i in entries:out[n+j]|=1<<i
   else:out[i]|=1<<(n+j)
 return out
words=[w for l in range(3) for w in itertools.product('OZJ',repeat=l)]
for w in words:
 g,r=graph(w);actual,s,c=inv(g,r);need(actual==profile(w),('source',w,actual,profile(w)));need(s==1,('start',w))
checks=0
for a in words:
 for b in words:
  dag=tuple({'O':'Z','Z':'O','J':'J'}[x] for x in a[::-1]);g,r=graph(b+dag+('J',));actual,s,c=inv(g)
  need(actual[0]==K(profile(a),profile(b)) and c==1,('C1',a,b,actual[0],c))
  for kind in ('N','M'):
   g=c2(a,b,kind);actual,s,c=inv(g)
   need(actual[0]==NM(profile(a),profile(b),kind) and c==2,(kind,a,b,actual[0],c))
  checks+=1
for p in (2,3,5,7):
 valid=[]
 for es in itertools.product(range(p),repeat=4):
  if all((e*e-e)%p==0 for e in es) and all(es[i]*es[j]%p==0 for i in range(4) for j in range(i)) and (sum(es)-1)%p==0:valid.append(es)
 need(len(valid)==4,('selectors',p,valid))
print(json.dumps({'source_words':len(words),'source_pair_checks_each_C1_N_M':checks,'largest_graph_order':10,'selector_primes':[2,3,5,7],'result':'PASS'},sort_keys=True,indent=2))
