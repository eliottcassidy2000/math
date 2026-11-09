import json,itertools,math,hashlib,sys,argparse
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,default=Path(__file__).resolve().parent.parent/'graft')
parser.add_argument('--sentinel',action='store_true')
args=parser.parse_args()
SRC=args.source
def need(x,m):
 if not x: raise RuntimeError(m)
if args.sentinel: need(False,'intentional audit failure sentinel')
def literal(adj,mask):
 n=len(adj); hs=[0]*(1<<n); ss=hs.copy(); es=hs.copy(); hs[0]=1
 # Independent ordered tuples of vertices, not endpoint recurrence.
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(adj[a]>>b&1 for a,b in zip(p,p[1:])):
    s=sum(1<<v for v in p);hs[s]+=1;ss[s]+=mask>>p[0]&1;es[s]+=mask>>p[-1]&1
 f=(1<<n)-1
 return (hs[f],sum(hs[s]*hs[f^s] for s in range(f+1)),sum(hs[s]*ss[f^s] for s in range(f+1)),sum(es[s]*hs[f^s] for s in range(f+1)),sum(es[s]*ss[f^s] for s in range(f+1)),ss[f],es[f])
def product(x,y):
 h,q,a,g,t,s,e=x;H,Q,A,G,T,S,E=y
 return (h*H,q*Q,a*Q+h*A,g*H+q*G,t*H+h*T+a*G+e*S,s*H,h*E)
def singleton(b):return (1,2,b,b,0,b,b)
def fold(xs):
 it=iter(xs);v=tuple(next(it))
 for x in it:v=product(v,x)
 return v
def word(w):return fold(singleton(int(b)) for b in w)
def ordinal(a,b):return [v|(((1<<len(b))-1)<<len(a)) for v in a]+[v<<len(a) for v in b]
def close(a,mask):
 n=len(a);return [v|(1<<n)|((not(mask>>i&1))<<(n+1)) for i,v in enumerate(a)]+[1<<(n+1),mask]
def count(adj):
 n=len(adj);N=1<<n;dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   z=dp[s][v]
   if z:
    nxt=adj[v]&(~s)&(N-1)
    while nxt:
     b=nxt&-nxt;nxt-=b;dp[s|b][b.bit_length()-1]+=z
 return sum(dp[-1])
cores={h:json.loads((SRC/f'core{h}.json').read_text()) for h in (37,45)}
for h,c in cores.items():
 n=len(c['adj']);need(all(not(c['adj'][i]>>i&1) for i in range(n)),'loops')
 need(all(((c['adj'][i]>>j&1)+(c['adj'][j]>>i&1))==1 for i in range(n) for j in range(i)),'tournament')
 need(literal(c['adj'],c['mask'])==tuple(c['core']),'permutation core')
 need(len(c['X'])==len(c['Y'])==(1188 if h==37 else 12),'length')
 end=int(h==45);c['x']=word(c['X']);c['y']=word(c['Y'])
 d=[v-u for u,v in zip(fold([singleton(0),c['core'],c['y'],singleton(end)]),fold([singleton(0),c['x'],c['core'],singleton(end)]))]
 need(d==[0,0,0,0,-2 if h==37 else 2,0,0],'identity')
 need(c['X'][-1]=='0' and c['Y'][0]==str(end),'flanks')
 need(c['core'][1]%h==0,'saturation')
 # Independent actual closure count of each core with no protective prefix.
 state=c['core'];need(count(close(c['adj'],c['mask']))==2*state[2]+h-state[4],'literal closure')
# Small literal ordinal laws, all possible cut masks for a triangle + pair.
a=[2,4,1];b=[2,0]
for ma in range(8):
 for mb in range(4):need(literal(ordinal(a,b),ma|(mb<<3))==product(literal(a,ma),literal(b,mb)),'literal ordinal law')
# Independent exhaustive fiber-residue selector: choose smallest legal j.
def select(a,b,K,L,t):
 js=[j for j in range(min(a,L+1)) if (t-b*j)%a==0]
 need(len(js)==1,'unique residue');j=js[0]
 # Shift to legal i by a units of j. This loop differs from producer ceiling formula.
 while t-b*j>a*K:j+=a
 i=(t-b*j)//a
 need(0<=i<=K and 0<=j<=L and a*i+b*j==t,'bounds')
 return i,j
cases=0
for a in range(2,22):
 for b in range(2,22):
  if math.gcd(a,b)>1:continue
  for K,L in [(b-1,a-1),(b,a+4),(2*b,2*a)]:
   C=(a-1)*(b-1);M=a*K+b*L;vals={a*i+b*j for i in range(K+1) for j in range(L+1)}
   need(C-1 not in vals and M-C+1 not in vals,'sharp omitted boundaries')
   for t in range(C,M-C+1):need(t in vals,'coverage');select(a,b,K,L,t)
   cases+=1
# Iterate words rather than producer exponentiation.
def family(h,k,j):
 c=cores[h];return fold([singleton(1),singleton(0)]+[c['x']]*j+[c['core']]+[c['y']]*(k-j)+[singleton(int(h==45))])
def H(s):return 2*s[2]+s[0]-s[4]
records=[]
for K,L in [(36,44),(36,50),(36,100),(37,44),(40,50)]:
 left=[family(37,K,i) for i in range(K+1)];right=[family(45,L,j) for j in range(L+1)]
 base=H(product(left[0],right[L]));lo=base+3168;M=45*K+37*L
 for i in range(K+1):
  for u in range(L+1):
   s=product(left[i],right[L-u]);need(H(s)==base+2*(45*i+37*u),'entire physical rectangle')
   need(s[:4]+s[5:]==product(left[0],right[0])[:4]+product(left[0],right[0])[5:],'common interface')
 for t in range(1584,M-1584+1):
  i,u=select(45,37,K,L,t);need(H(product(left[i],right[L-u]))==base+2*t,'every selected target')
 need(lo%2==1,'oddness');records.append({'K':K,'L':L,'lower_sha256':hashlib.sha256(lo.to_bytes((lo.bit_length()+7)//8,'big')).hexdigest(),'width':M-3167,'vertices':1188*K+12*L+23})
source=json.loads((SRC/'CHECK.json').read_text())
for r,s in zip(records,source['records']):need(r['lower_sha256']==s['lower_sha256'],'producer receipt match')
print(json.dumps({'status':'PASS','permutation_cores':[37,45],'literal_ordinal_cases':32,'literal_core_closures':2,'coin_rectangles':cases,'records':records},indent=2))
