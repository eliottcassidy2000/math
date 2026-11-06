import itertools
from functools import lru_cache

def J(a):
 n=len(a)
 return relabel(tuple([a[0]]+[x|(1<<n) for x in a[1:]]+[1]),[n]+list(range(n)))
def relabel(a,p):return tuple(sum(1<<j for j,y in enumerate(p) if a[x]>>y&1)for x in p)
def gate(a,g):
 n=len(a)
 # input root 0, new vertices a,c,r as applicable, original nonroot beats every new.
 if g==1:
  aa,rr=n,n+1;b=list(a)+[1,1|(1<<aa)];b[0]&=~((1<<aa)|(1<<rr))
  for i in range(1,n):b[i]|=(1<<aa)|(1<<rr)
 else:
  cc,aa,rr=n,n+1,n+2;b=list(a)+[1|(1<<rr),(1<<cc),1|(1<<aa)];b[0]|=1<<aa
  for i in range(1,n):b[i]|=(1<<cc)|(1<<aa)|(1<<rr)
 return relabel(tuple(b),[rr]+list(range(rr)))

@lru_cache(None)
def profiles(a,K=3):
 n=len(a);M=1<<n;F=M-1
 hp=[[0]*n for _ in range(M)];sp=[[0]*n for _ in range(M)];H=[0]*M;H[0]=1
 for i in range(n):hp[1<<i][i]=sp[1<<i][i]=1
 for m in range(1,M):
  for i in range(n):
   if not(m>>i&1):continue
   for j in range(n):
    if m>>j&1:continue
    if a[i]>>j&1:hp[m|1<<j][j]+=hp[m][i]
    else:sp[m|1<<j][j]+=sp[m][i]
  H[m]=sum(hp[m])
 Cs=[[int(m==0)for m in range(M)],H]
 for k in range(2,K+1):
  cc=[0]*M
  for m in range(M):
   s=m
   while True:
    cc[m]+=Cs[-1][s]*H[m^s]
    if s==0:break
    s=(s-1)&m
  Cs.append(cc)
 ans=[]
 for k in range(1,K+1):
  L=sum(hp[m][0]*Cs[k-1][F^m]for m in range(1,M,2))
  R=sum(sp[m][0]*Cs[k][F^m]for m in range(1,M,2))
  ans.append((Cs[k][F],L,R,Cs[k][F^1]))
 return tuple(ans),sp[F][0]

def labelled(n):
 ps=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(ps)):
  a=[0]*n
  for i,(u,v) in enumerate(ps):
   if bits>>i&1:a[u]|=1<<v
   else:a[v]|=1<<u
  yield tuple(a)


def nine(a):
 p,s=profiles(tuple(a),2)
 (h,e,R,d),(C,L,T,b)=p
 return(h,s,e,L,R,d,C,b,T)

def check_tournament(a):
 n=len(a)
 if not n:raise ValueError('nonempty tournament required')
 for i in range(n):
  if a[i]<0 or a[i]>>n or a[i]>>i&1:raise ValueError('bad row')
  for j in range(i):
   if ((a[i]>>j&1)+(a[j]>>i&1))!=1:raise ValueError('invalid edge')
 return True

def literal(word):
 a=(2,0)
 for op in word:
  n=len(a)
  if op=='O':a=tuple([x|(1<<n)for x in a]+[0])
  elif op=='Z':a=tuple([a[0]]+[x|(1<<n)for x in a[1:]]+[1])
  elif op=='J':a=J(a)
  elif op=='K1':a=gate(a,1)
  elif op=='K2':a=gate(a,2)
  else:raise ValueError(op)
 return a
