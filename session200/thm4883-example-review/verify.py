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

from math import gcd
result={}
for w in ('ZOJ','ZOJJ'):
 rows,root=graph(w);actual,starts,cycles=inv(rows,root)
 need(actual==profile(w) and starts==cycles==1,('actual strong source',w,actual,starts,cycles))
 result[w]={'rows':rows,'root':root,'profile':actual,'start':starts,'cycles':cycles}
u=profile('ZOJ');v=profile('ZOJJ');q=profile('Z')
need(u==(9,3,9,3) and v==(15,6,18,9),'claimed profiles')
D,E,U,d=v;alpha=D+4*U-2*d;beta=alpha+E-d;c=-D+U+d+2*E;gamma=gcd(alpha,d-E)
need((alpha,beta,c,gamma)==(69,66,24,3),'offset coefficients')
for kind,left,right in [('N','ZOJJ','ZO'),('M','ZOJJ','ZJ')]:
 rows=c2(left,right,kind);actual,starts,cycles=inv(rows)
 scalar=NM(v,profile(right),kind)
 need(actual[0]==scalar and cycles==2,('literal graft',kind,actual[0],cycles))
 result[kind]={'rows':rows,'order':len(rows),'H':actual[0],'C':cycles}
need(result['N']['H']==417 and result['M']['H']==867 and 867-2*417==33,'graft counts')
# Coefficient verification of the affine-in-(D,U), linear-in-v identity.
for i in range(4):
 v0=tuple(int(i==j) for j in range(4));A,B,V,b=v0
 need(move(v0,'J')[3]-move(v0,'J')[1]==B,'correct J identity')
 aa=A+4*V-2*b;bb=aa+B-b;cc=-A+V+b+2*B
 for DD,UU in [(0,0),(1,0),(0,1)]:
  q0=(DD,1,UU,1)
  need(NM(v0,move(q0,'J'),'M')-2*NM(v0,move(q0,'O'),'N')==aa*DD-bb*UU+cc,'full T4 formula')
for DD in range(1,100,2):
 for UU in range(1,100,2):need((69*DD-66*UU+24)%6==3,'odd-input residue obstruction')
result.update({'status':'PASS','alpha':alpha,'beta':beta,'constant':c,'gcd':gamma,'correct_J_gap_identity':'d(Ju)-E(Ju)=E(u)','offset':33,'offset_class':'3 modulo 6; offsets +1 and -1 impossible','scope':'Refutes the unqualified J-born example only. The T4 identity and gcd-qualified class survive.'})
print(json.dumps(result,sort_keys=True,indent=2))
