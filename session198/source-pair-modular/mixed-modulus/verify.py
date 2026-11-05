#!/usr/bin/env python3
"""Standalone physical/modular controls; every check survives python -O."""
from pathlib import Path
import hashlib,json,math
P=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def gadget(L):
 n=(L-1).bit_length();out=[0]*(n+1)
 for i in range(n):
  x=i+1
  if (L-1)>>i&1:out[x]|=1
  else:out[0]|=1<<x
  for y in range(x+1,n+1):out[x]|=1<<y
 return out

def converse(T):
 full=(1<<len(T))-1
 return [full^(1<<i)^row for i,row in enumerate(T)]

def count(T,root=0):
 n=len(T);F=(1<<n)-1
 need(all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n)for j in range(i+1,n)), 'tournament')
 st=[[0]*n for _ in range(F+1)];en=[[0]*n for _ in range(F+1)];H=[1]+[0]*F
 for M in range(1,F+1):
  for x in range(n):
   if not M>>x&1:continue
   prev=M^(1<<x)
   if not prev:st[M][x]=en[M][x]=1
   else:
    st[M][x]=sum(st[prev][y] for y in range(n) if prev>>y&1 and T[x]>>y&1)
    en[M][x]=sum(en[prev][y] for y in range(n) if prev>>y&1 and T[y]>>x&1)
  H[M]=sum(en[M])
 R=sum(H[M]*st[F^M][root] for M in range(F+1))
 cyc=[[0]*n for _ in range(F+1)];cyc[1][0]=1
 for M in range(1,F+1,2):
  for x in range(n):
   if not cyc[M][x]:continue
   for y in range(n):
    if not M>>y&1 and T[x]>>y&1:cyc[M|1<<y][y]+=cyc[M][x]
 C=sum(cyc[F][x]for x in range(n)if T[x]&1)if n>=3 else 0
 return dict(H=H[F],s=st[F][root],E=en[F][root],U=R,d=H[F^(1<<root)],C=C)

def graft(L,kind):
 left=[10,12,1,4]if kind==2 else[2,12,1,5]
 G=gadget(L);out=left+[row<<4 for row in G]
 for x in range(4):
  for y in range(len(G)):
   if y==0 and x in (0,3):out[4]|=1<<x
   else:out[x]|=1<<(4+y)
 return out

def read_bits(word,mod=None):
 n=0;f=1
 for bit in word:
  if bit=='1':n=2*n+1
  else:f+=n+(n&1);n*=2
  if mod is not None:n%=mod;f%=mod
 return n,f

def gg(L):return read_bits(bin(L-1)[2:])[1]

def compiler(M,u,z):
 need(M>=2 and M%2==0 and z%2==1,'compiler input')
 k=0;q=M
 while q%2==0:k+=1;q//=2
 a=2
 while pow(2,a,q)!=1%q:a+=1
 delta=(pow(2,a)-1)*(pow(2,k)-1)+1
 need(math.gcd(delta,M)==2,'loop gcd')
 s=(u-1)%M;t=(z-gg(s+1))//2
 j=(t*pow((delta//2)%(M//2),-1,M//2))%(M//2)if M>2 else 0
 j+=M//2
 return ('1'*a+'0'*k)*j+format(s,'0'+str((M-1).bit_length())+'b')

def select(m,r,kind):
 R=r if r%2 else r+m
 need(R%2==1,'admissible residue')
 t=(R-3)//2;M=2*m
 u=(-6*t if kind==2 else-2*t)%M
 z=(1+(16 if kind==2 else 6)*t)%M
 return M,u,z

literal=[]
for L in range(2,65):
 G=gadget(L);w=count(converse(G));g=gg(L)
 need((w['H'],w['E'],w['U'],w['d'],w['s'])==(g,L%2,L,1,1),('source',L,w))
 for kind in(2,3):
  T=graft(L,kind);p=count(T)
  h=5*g+(13 if kind==2 else 14)*L-(3 if kind==2 else 2)*(L%2)-2
  need(p['H']==h and p['C']==kind,('graft',L,kind,p,h))
  literal.append(dict(L=L,cycles=kind,H=h,order=len(T)))
modular=0;maxbits=0;exact_small=0
for m in list(range(1,51))+[64,128,256]:
 for r in range(m):
  if m%2==0 and r%2==0:continue
  for kind in(2,3):
   M,u,z=select(m,r,kind);word=compiler(M,u,z)
   n,f=read_bits(word,M);L=(n+1)%M
   need(L==u and f==z and L%2==0,('joint',m,r,kind,L,f,u,z))
   h=(5*f+(13 if kind==2 else 14)*L-2)%m
   need(h==r,('target',m,r,kind,h))
   if m<=7:
    nn,ff=read_bits(word);actualL=nn+1
    need(actualL%2==0 and (5*ff+(13 if kind==2 else 14)*actualL-2)%m==r,'exact integer replay')
    exact_small+=1
   modular+=1;maxbits=max(maxbits,len(word))
for L in range(2,1025,2):
 need((3*gg(L)+6*L-2)%3==1,'N hostile')
 need((5*gg(L)+15*L-6)%5==4,'A hostile')
 need(5*gg(L)+13*L-2>=29 and 5*gg(L)+14*L-2>=31,'exact-value exclusion')
receipt=dict(status='PASS',physical_L_universe=[2,64],physical_source_checks=63,physical_graft_checks=len(literal),literal_graphs=literal,moduli=list(range(1,51))+[64,128,256],modular_checks=modular,max_compiled_word_bits=maxbits,exact_integer_replays=exact_small,hostile_even_L_cap=1024,dependencies={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((P/'dependencies').glob('*.md'))})
s=json.dumps(receipt,sort_keys=True,indent=2)+'\n'
(P/'VALIDATION.json').write_text(s);print(s,end='')
