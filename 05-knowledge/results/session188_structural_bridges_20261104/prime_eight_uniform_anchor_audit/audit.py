from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'prime_eight_uniform_anchor';Q=[238,4,240,118,195,146,34,74];QUERIES=[(1,0),(255,254)]
def need(c,m):
 if not c:raise RuntimeError(m)
def graph(sizes,rev,a,c):
 R=[255^(1<<i)^v for i,v in enumerate(Q)] if rev else Q;offs=[sum(sizes[:i]) for i in range(8)];n=sum(sizes);T=[0]*(n+2);T[n]=1<<(n+1)
 for i,m in enumerate(sizes):
  for u in range(m):
   v=offs[i]+u;T[v]=sum(1<<(offs[i]+w) for w in range(u+1,m))|sum(((1<<sizes[j])-1)<<offs[j] for j in range(8) if R[i]>>j&1)
   for mask,p in [(a,n),(c,n+1)]:
    if mask>>i&1:T[v]|=1<<p
    else:T[p]|=1<<v
 return T

def main():
 source=json.loads((SRC/'normal.json').read_text());shapes=[]
 for t in product(range(1,4),repeat=8):
  if sum(v>1 for v in t)<=2:shapes.append(t)
 shapes.sort(key=lambda t:sum(v<<(2*i) for i,v in enumerate(t)));need(len(shapes)==129,'complete run-shape domain');cases=[];keys=[]
 for a,c in QUERIES:
  for rev in [False,True]:
   for t in shapes:keys.append((a,c,rev,t));cases.append(graph(t,rev,a,c))
 lines=[]
 for T in cases:lines.extend([str(len(T)),' '.join(map(str,T))])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'interpolate.cpp'),'-o',str(exe)],check=True);values=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).decode().splitlines()]
 need(len(values)==516,'all actual physical cover profiles');lookup=dict(zip(keys,values));operators={1:({1:1},1),2:({1:-1,2:1},2),3:({1:2,2:-3,3:1},6)};kernel=[]
 for a,c in QUERIES:
  for t in shapes:
   sides=[]
   for rev in [False,True]:
    den=1;choices=[]
    for r in t:op,d=operators[r];choices.append(list(op.items()));den*=d
    vv=[0,0,0]
    for factors in product(*choices):
     ss=tuple(s for s,w in factors);w=__import__('math').prod(w for s,w in factors)
     for h,v in enumerate(lookup[a,c,rev,ss]):vv[h]+=w*v
    need(all(v%den==0 for v in vv),'exact tensor profile interpolation divisions');sides.append([v//den%8 for v in vv])
   code=sum(v<<(2*i) for i,v in enumerate(t));row=[a,c,code]+[sides[s][h] for h in range(3) for s in [0,1]];kernel.append(row)
 need(kernel==source['full_kernels'],'every2x3x129 kernel entry independently matches');symbols=[]
 for a,c in QUERIES:
  rows=[r for r in kernel if r[:2]==[a,c]]
  for h in range(3):
   P={}
   for row in rows:
    t=[(row[2]>>(2*i))&3 for i in range(8)];terms={0:(row[3+2*h]-row[4+2*h])%8}
    for i,k in enumerate(t):
     if k==1:continue
     factors=[(0,2),(1<<(2*i),4)] if k==2 else[(1<<(2*i),4),(1<<(2*i+1),2)];out={}
     for mask,v in terms.items():
      for m,w in factors:out[mask|m]=(out.get(mask|m,0)+v*w)%8
     terms=out
    for m,v in terms.items():P[m]=(P.get(m,0)+v)%8
   need({m:v for m,v in P.items() if v}==({0:4} if h==0 else{}),'complete formal coefficient identity');symbols.append({'query':[a,c],'h':h+1,'all_coefficients':[[m,P.get(m,0)] for m in range(1<<16) if m.bit_count()<=2]})
 # Complete primality and endpoint-balanced checks on the fixed quotient.
 modules=[M for M in range(256) if 2<=M.bit_count()<8 and all(Q[v]&M in(0,M) for v in range(8) if not M>>v&1)];need(not modules,'fixed quotient prime')
 @lru_cache(None)
 def hp(M,u,v):
  if M==1<<v:return int(u==v)
  S=M^(1<<v);return sum(hp(S,u,w) for w in range(8) if S>>w&1 and Q[w]>>v&1)
 starts=[sum(hp(255,u,v) for v in range(8)) for u in range(8)];ends=[sum(hp(255,u,v) for u in range(8)) for v in range(8)];sm=sum((v%2)<<i for i,v in enumerate(starts));em=sum((v%2)<<i for i,v in enumerate(ends));need(sm==em==55,'endpoint-balanced obstruction independently confirmed')
 need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer receipts');pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','covers.cpp','physical.cpp','normal.json','optimized.json']}
 print(json.dumps({'status':'PASS','physical_models':516,'maximum_physical_order':14,'complete_physical_cover_values':values,'full_independent_kernels':kernel,'complete_formal_coefficients':symbols,'quotient_start_counts':starts,'quotient_end_counts':ends,'endpoint_parity_mask':sm,'source_pins':pins,'audit_source_sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['audit.py','interpolate.cpp']},'scope':'One specified prime8 quotient, arbitrary nontrivial children and fixed ambient; full ordered-cover certificate.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
