import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
"""Exact certificates: source states by independent path permutations and subset DP.
Run python3 -B verify.py > CHECK.json (also -O). No dependency on discovery census.
"""
import json,itertools
from core import states,compose

def need(v,msg):
 if not v: raise RuntimeError(msg)

def paths(adj,vs):
 return [p for p in itertools.permutations(vs) if all(adj[u]>>v&1 for u,v in zip(p,p[1:]))]
def literal(adj,cut):
 n=len(adj);F=(1<<n)-1;pp=[paths(adj,[v for v in range(n) if s>>v&1]) for s in range(F+1)]
 h=len(pp[F]);q=A=G=T=0
 for s in range(F+1):
  for p in pp[s]:
   for r in pp[F^s]:
    q+=1;a=bool(r and cut>>r[0]&1);g=bool(p and cut>>p[-1]&1)
    A+=a;G+=g;T+=a*g
 return h,q,A,G,T,sum(cut>>p[0]&1 for p in pp[F]),sum(cut>>p[-1]&1 for p in pp[F])
def count(adj,cut=None):
 n=len(adj);N=1<<n;dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(N):
  for u,z in enumerate(dp[s]):
   if not z:continue
   choices=adj[u]&~s
   while choices:
    bit=choices&-choices;choices-=bit;dp[s|bit][bit.bit_length()-1]+=z
 return sum(dp[-1]) if cut is None else sum(z for i,z in enumerate(dp[-1]) if cut>>i&1)
def graph(C,cut,P,b):
 bits=list(map(int,bin(P)[2:]));l=len(bits);a=[0]*l
 for i in range(l):a[i]=sum(1<<j for j in range(i+1,l+8+b))
 a += [(v<<l)|sum(1<<j for j in range(l+8,l+8+b)) for v in C]
 for i in range(l+8,l+8+b):a.append(sum(1<<j for j in range(i+1,l+8+b)))
 mask=sum(v<<i for i,v in enumerate(bits))|(cut<<l)
 return a,mask

def closure(a,mask):
 n=len(a);a=a[:]+[1<<(n+1),mask]
 for i in range(n):a[i]|=1<<n;a[i]|=(0 if mask>>i&1 else 1<<(n+1))
 for v in range(n+2):
  seen=1<<v
  while True:
   new=seen
   for u in range(n+2):
    if seen>>u&1:new|=a[u]
   if new==seen:break
   seen=new
  need(seen==(1<<(n+2))-1,'strong')
 return a

def f(P):
 s=bin(P)[2:]
 return sum((1<<(j-i-1))+(j==i+1) for i,a in enumerate(s) if a=='1' for j in range(i+1,len(s)) if s[j]=='1')
C=[4,65,26,99,43,135,181,31];D=[68,9,50,21,35,139,62,95]
cores=[(C,189),(D,190)];out={'cores':[],'literal_closures':[],'large_parameter_controls':[]}
for core,cut in cores:
 st=states(core)[cut];need(st==literal(core,cut),'literal core disagrees');out['cores'].append({'adj':core,'cut':cut,'state':st})
 for b in [0,1]:
  a,m=graph(core,cut,25,b);g=closure(a,m);H=count(g)
  c=98861 if core==C else 98863
  need(H==(1<<(18+b))-c,'closure count');out['literal_closures'].append({'b':b,'vertices':len(g),'H':H})
 # Source q=A is audited by ordinary counts on the root extension and its first zero extension.
 a,m=graph(core,cut,25,0);n=len(a);root=a[:]+[m]
 for i in range(n):
  if not(m>>i&1):root[i]|=1<<n
 F=count(root);g=count(a,m);q=count(closure(a,m))-F-g
 need(q==1<<17,'source slope');need(q-F-g==(98861 if core==C else 98863),'source deficit')
 out['cores'][-1]['source']={'F':F,'g':g,'q':q,'h':323}
M=1<<84;B=(M-1)*459//637;need((M-1)%637==0,'period');need(pow(2,17,5096)==3672,'entry');need(len(bin(B)[2:])==84,'word length');r=int(bin(B)[2:][::-1],2)
P=25;ft=14
for t in range(8):
 m=17+84*t;need(P==((1<<m)-3672)//5096 and P%2==1,'prefix');need(f(P)==ft,'repeat cocycle');c=3672*P+323*ft+2539
 need((1<<(m-1))<c<c+2<(1<<m),'cost interval');need(P.bit_length()==m-12,'order')
 # Compose genuine singleton states and the actual core; no huge graph subset DP.
 for core,cut in cores:
  st=None
  for bit in bin(P)[2:]:
   v=int(bit);s=(1,2,v,v,0,v,v);st=compose(st,s) if st else s
  st=compose(st,states(core)[cut]);need(st[2]==1<<m,'physical A');need(st[4]-st[0]==c+(2 if core==D else 0),'physical deficit')
 out['large_parameter_controls'].append({'t':t,'m':m,'c':c,'prefix_length':P.bit_length()})
 ft+=P*r+f(B)+1;P=M*P+B
out['block']={'M':M,'B':B,'reverse_B':r,'f_B':f(B),'drift_numerator':5096*(f(B)+1)-3672*r}
# Hostile: suffix cannot create a pure dyadic A when prefix-core A is nondyadic.
for h in range(1,16,2):
 for Z in range(1,20):
  for b in range(6):
   for R in range(1<<b):
    A=(1<<b)*Z+h*R
    if A&(A-1)==0:need(R==0 and Z&(Z-1)==0,'suffix rigidity')
# Wrong source threshold would claim value32211 at k17; our source requires positive zero length.
need((1<<17)-98861==32211,'hostile threshold');out['scope']='Parameterized twin deficits and optimal dyadic source costs; no growing interval or universal coverage.'
print(json.dumps(out,sort_keys=True,indent=2))
