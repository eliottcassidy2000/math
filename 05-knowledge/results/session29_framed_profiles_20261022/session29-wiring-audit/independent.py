import itertools as it,json,hashlib
from pathlib import Path
import sympy as sp
D=Path(__file__).parent; P=D.parent/'session29-profiles'
def need(x,msg):
 if not x: raise RuntimeError(msg)
def graph(mask,n):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(it.combinations(range(n),2)):a[i][j]=bool(mask>>k&1);a[j][i]=not a[i][j]
 return a
def encode(a):return sum(int(a[i][j])<<k for k,(i,j) in enumerate(it.combinations(range(len(a)),2)))
def tiles(mask):
 a=graph(mask,5);out=set()
 for p in it.permutations(range(5)):
  if all(a[p[i]][p[i+1]] for i in range(4)):out.add(encode([[a[i][j] for j in p] for i in p]))
 return sorted(out)
def flip(a):
 n=len(a);return [[False if i==j else a[i][j] if abs(i-j)==1 else not a[i][j] for j in range(n)] for i in range(n)]
def mirror(t):
 a=graph(t,5);return encode([[False if i==j else not a[4-i][4-j] for j in range(5)] for i in range(5)])
def dp(a):
 n=len(a);end=[[0]*n for _ in range(1<<n)];out=[sum(int(x)<<j for j,x in enumerate(row)) for row in a]
 for v in range(n):end[1<<v][v]=1
 for s,row in enumerate(end):
  for v,x in enumerate(row):
   if not x:continue
   b=out[v]&~s
   while b:
    bit=b&-b;b-=bit;end[s|bit][bit.bit_length()-1]+=x
 return end
# States ordered S,E,I,active nonsingleton,active singleton.
def signature(a):
 n=len(a);f=(1<<n)-1;p=n-1;e=dp(a);h=list(map(sum,e));st=dp([list(x) for x in zip(*a)])[f][p];en=e[f][p];singleton=h[f^(1<<p)]
 two=sum(e[s][p]*h[f^s] for s in range(1,f) if s>>p&1)
 return [st,en,h[f]-st-en,two-singleton,singleton]
def atomgraph(a,state):
 q=1 if state<3 else 2;n=q+len(a);b=[[False]*n for _ in range(n)]
 for i in range(len(a)):
  for j in range(len(a)):b[q+i][q+j]=a[i][j]
 # atom 0 contains portal. Incoming block arcs exist except deleted first->portal.
 for atom in range(q):
  for v in range(len(a)):b[q+v][atom]=not(atom==0 and v==0 and state in (0,4))
 if state in (1,3,4):b[0][q]=True
 return b
def compressed(rows):
 # output map H=S+E+I; S=S; R=E+DE+DB
 cols=[sp.Matrix([sum(row[:3]),row[0],row[1]+row[3]+row[4]]) for row in rows]
 need(cols[3]==cols[4],'singleton collapse')
 need(cols[1]-cols[2]==cols[3],'E Q collapse')
 return sp.Matrix.hstack(cols[2],cols[0]-cols[2],cols[3])
def joined(ts):
 n=5*len(ts);aa=[graph(t,5) for t in ts]
 return [[False if i==j else aa[i//5][i%5][j%5] if i//5==j//5 else i//5<j//5 for j in range(n)] for i in range(n)]
letters=tiles(40)+tiles(41);need(len(set(letters))==18,'alphabet18')
raw={t:[signature(atomgraph(flip(graph(t,5)),s)) for s in range(5)] for t in letters}
K={t:compressed(raw[t]) for t in letters};G=sp.diag(1,-1,1);v=sp.Matrix([1,0,0]);tau={t:mirror(t) for t in letters}
for t in letters:
 sig=signature(flip(graph(t,5)));need(K[t]*v==sp.Matrix([sum(sig[:3]),sig[0],sig[1]+sum(sig[3:])]),'empty initial')
 need(K[tau[t]]==G*K[t].T*G,'mirror adjoint')
producer=json.loads((P/'transfer-results.json').read_text());need(all(raw[t]==producer['matrices'][str(t)] for t in letters),'producer raw equality')
rows=[]
for t,u in it.product(letters,repeat=2):
 a=flip(joined((t,u)));h=sum(dp(a)[-1]);pred=(K[u]*K[t]*v)[0];need(h==pred,'r2 exact DP')
 need(all(a[i][i+1] for i in range(9)) and a[9][0],'strong Hamilton cycle')
 rows.append({'tiles':[t,u],'H':h,'blue':u==tau[t]})
checks=[]
for inds in ((0,0,0),(0,5,0),(5,0,5),(17,17,17),(1,8,16),(4,12,7)):
 ts=[letters[i] for i in inds];h=sum(dp(flip(joined(ts)))[-1]);pred=(K[ts[2]]*K[ts[1]]*K[ts[0]]*v)[0];need(h==pred,'r3 exact DP');checks.append([ts,h])
A=sum(K.values(),sp.zeros(3));z=sp.symbols('z');F=sp.cancel((v.T*(sp.eye(3)-z*A).inv()*v)[0]-1)
need(F==sp.cancel((172*z+45144*z*z+50400*z**3)/(1-874*z-43344*z*z-50400*z**3)),'GF')
need(G*A==A.T*G,'aggregate adjoint')
# Independent tensor-square moment operator; no fitted recurrence.
x=sp.symbols('h s r');xx=sp.Matrix(x);mons=[x[0]**2,x[1]**2,x[2]**2,x[0]*x[1],x[0]*x[2],x[1]*x[2]]
def coeff(p):return [sp.Poly(sp.expand(p),*x).coeff_monomial(m) for m in mons]
L=sp.zeros(6)
for M in K.values():
 y=M*xx;L+=sp.Matrix([coeff(p) for p in [y[0]**2,y[1]**2,y[2]**2,y[0]*y[1],y[0]*y[2],y[1]*y[2]]])
w=sp.symbols('w');seed=sp.Matrix([1,0,0,0,0,0]);fixed=sum((K[t] for t in letters if tau[t]==t),sp.zeros(3));res=(sp.eye(6)-w*L).inv()*seed
even=sp.cancel((sp.Matrix([[1,-1,1,0,0,0]])*res)[0]);odd=sp.cancel((sp.Matrix([coeff((xx.T*G*fixed*xx)[0])])*res)[0])
need(int((sp.Matrix([[1,-1,1,0,0,0]])*L*seed)[0])==sum(r['H'] for r in rows if r['blue']),'blue r2')
result={'alphabet':letters,'fixed_letters':[t for t in letters if tau[t]==t],'matrix':A.tolist(),'GF':str(F),'two_block_checks':len(rows),'two_block_H':sum(r['H'] for r in rows),'two_block_blue_H':sum(r['H'] for r in rows if r['blue']),'three_block_checks':checks,'mirror_adjoint_checks':18,'moment_matrix':L.tolist(),'even_GF':str(even),'odd_GF':str(odd),'raw_matrices':raw,'compressed':{t:M.tolist() for t,M in K.items()}}
(D/'evidence.json').write_text(json.dumps(result,default=int,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('raw_matrices','compressed','moment_matrix')},default=int,indent=2))
# Exhaust all 72 mirror-fixed length-three words against the full 15-vertex graph.
blue3=0
for t in letters:
 for u in letters:
  if tau[u]!=u:continue
  ts=(t,u,tau[t]);actual=sum(dp(flip(joined(ts)))[-1]);pred=(K[tau[t]]*K[u]*K[t]*v)[0]
  need(actual==pred,'all72 blue r3');blue3+=actual
cost=sp.Matrix([coeff((xx.T*G*fixed*xx)[0])]);need(blue3==(cost*L*seed)[0],'blue3 moment')
result['all72_blue_r3_H']=blue3
(D/'evidence.json').write_text(json.dumps(result,default=int,indent=2)+'\n');print('ALL72 BLUE R3 PASS',blue3)
