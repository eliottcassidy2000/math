"""Independent explicit-cover, run-word and physical interpolation audit."""
from pathlib import Path
from itertools import permutations,product,combinations
from functools import lru_cache
import hashlib,json,math,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'marked_quartic'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def opposite(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def induce(T,names):return tuple(sum(1<<j for j,w in enumerate(names) if T[u]>>w&1) for u in names)
def tournament(n,code):
 T=[0]*n
 for j,(u,v) in enumerate(combinations(range(n),2)):
  a,b=(u,v) if code>>j&1 else (v,u);T[a]|=1<<b
 return tuple(T)
@lru_cache(None)
def marked(T):
 n=len(T);q=[0]*(n+1);E=[[0]*(n+1) for _ in T];S=[[0]*(n+1) for _ in T]
 for word in permutations(range(n)):
  bad=sum(1<<j for j,(u,v) in enumerate(zip(word,word[1:])) if not T[u]>>v&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   k=1+cuts.bit_count();last=(cuts&-cuts).bit_length()-1 if cuts else n-1
   q[k]+=1;E[word[last]][k]+=1;S[word[0]][k]+=1
 return q,E,S
@lru_cache(None)
def runkernel(Q,caps,i,j):
 n=len(Q);cols=[{} for _ in Q];rows=[{} for _ in Q]
 def walk(first,last,counts):
  key=tuple(counts);rows[first][key]=rows[first].get(key,0)+1;cols[last][key]=cols[last].get(key,0)+1
  for nxt in range(n):
   if Q[last]>>nxt&1 and counts[nxt]<caps[nxt]:
    counts[nxt]+=1;walk(first,nxt,counts);counts[nxt]-=1
 for first in range(n):
  counts=[0]*n;counts[first]=1;walk(first,first,counts)
 K={}
 for sign,L,R in ((1,cols[i],rows[j]),(-1,cols[j],rows[i])):
  for a,x in L.items():
   for b,y in R.items():
    key=tuple(u+v for u,v in zip(a,b))
    if all(u<=v for u,v in zip(key,caps)):K[key]=K.get(key,0)+sign*x*y
 return {e:c for e,c in K.items() if c}
def evaluate(K,ps):return sum(c*math.prod(p[k] for p,k in zip(ps,e)) for e,c in K.items())
def inflate(Q,blocks):
 off=[sum(len(A) for A in blocks[:i]) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off

def endpoints(T):
 n=len(T);D=[[0]*n for _ in range(1<<n)];P=opposite(T)
 for u in range(n):D[1<<u][u]=1
 for M in range(1,1<<n):
  rem=M
  while rem:
   b=rem&-rem;rem-=b;u=b.bit_length()-1;S=M^b
   if S:D[M][u]=sum(D[S][v] for v in range(n) if S>>v&1 and P[u]>>v&1)
 return D
def Cmatrix(T):
 E=endpoints(T);S=endpoints(opposite(T));N=(1<<len(T))-1
 return [[sum(E[M][u]*S[N^M][v] for M in range(1,N)) for v in range(len(T))] for u in range(len(T))]
def prediction(Q,blocks,i,u,j,v):
 ps=[marked(tuple(A)) for A in blocks];K=runkernel(tuple(Q),tuple(map(len,blocks)),i,j);P=[x[0] for x in ps];P[i]=ps[i][1][u];P[j]=ps[j][2][v];a=evaluate(K,P);P[i]=ps[i][2][u];P[j]=ps[j][1][v];return a-evaluate(K,P)
def main():
 paths=[SRC/name for name in ('PROOF.md','verify.py','literal.cpp','literal.json','receipt.json')]+[BASE/'response'/'PROOF.md',BASE/'covered_response'/'verify.py',BASE/'covered_response'/'full_response.json'];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'receipt.json').read_text());profiles=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   T=tournament(n,code);q,E,S=marked(T)
   need(all(sum(r[k] for r in E)==sum(r[k] for r in S)==q[k] for k in range(n+1)),'explicit-cover normalization')
   if n>=2:
    for u in range(n):need(E[u][n-1]-S[u][n-1]==math.factorial(n-2)*(n-1-2*T[u].bit_count()),'explicit-cover top score formula')
   profiles+=1
 # Arbitrary quotient and child interiors, with all cross-child physical marks.
 childtypes=((0,),(2,0),(0,1),(2,4,1));factorizations=0;configurations=0
 for code in range(8):
  Q=tournament(3,code)
  for blocks in product(childtypes,repeat=3):
   T,off=inflate(Q,blocks);U,_=inflate(opposite(Q),blocks);C,D=Cmatrix(T),Cmatrix(U)
   for i,j in combinations(range(3),2):
    for u in range(len(blocks[i])):
     for v in range(len(blocks[j])):
      p,q=off[i]+u,off[j]+v;actual=C[p][q]+C[q][p]-D[p][q]-D[q][p];need(prediction(Q,blocks,i,u,j,v)==actual,'all-order factorization finite independent control');factorizations+=1
   configurations+=1
 cases=[];inputs=[];expected=[]
 literal=json.loads((SRC/'literal.json').read_text())
 for rec,lit in zip(data['exact_cases'],literal['cases']):
  Q=tuple(rec['quotient']);blocks=list(map(tuple,rec['child_rows']));T,off=inflate(Q,blocks);U,_=inflate(opposite(Q),blocks);need(T==tuple(rec['T']) and U==tuple(rec['U']),'physical case reconstruction');w,u,v=rec['deleted_and_marks'];child=lambda x:next(i for i,o in reversed(list(enumerate(off))) if x>=o);h,i,j=map(child,(w,u,v));names=[list(range(len(A))) for A in blocks];names[h].remove(w-off[h]);cards=[induce(A,ns) for A,ns in zip(blocks,names)];mu=names[i].index(u-off[i]);mv=names[j].index(v-off[j]);K=runkernel(Q,tuple(map(len,cards)),i,j)
  need([[*e,c] for e,c in sorted(K.items())]==rec['kernel_terms'],'independent complete run kernel');ps=[marked(A) for A in cards]
  need([x[0] for x in ps]==rec['child_card_profiles'],'explicit card q profiles');need([ps[i][1][mu],ps[i][2][mu],ps[j][1][mv],ps[j][2][mv]]==[rec[k] for k in ('first_end_profile','first_start_profile','second_end_profile','second_start_profile')],'explicit marked-card endpoint profiles');coefficient=prediction(Q,cards,i,mu,j,mv);need(coefficient==rec['shared_quartic_difference'],'named shared-label value');cases.append({'physical_marks':[w,u,v],'coefficient':coefficient,'kernel_terms':len(K)})
  need(len(lit['values'])==16,'all actual interpolation queries retained')
  for a,c,f,g in lit['values']:
   for R,z in ((T,f),(U,g)):inputs.append(str(len(R))+' '+' '.join(map(str,R))+f' {a} {c}');expected.append(z)
 with tempfile.TemporaryDirectory(prefix='marked-quartic-independent-') as tmp:
  exe=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True,capture_output=True);out=subprocess.check_output([str(exe)],input=(str(len(inputs))+'\n'+'\n'.join(inputs)+'\n').encode());counts=list(map(int,out.split()));need(counts==expected,'independent predecessor-DP physical query values')
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'receipt.json').read_bytes()==(SRC/'receipt.optimized.json').read_bytes(),'producer mode replay');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'frozen producer sources')
 print(json.dumps({'status':'PASS','explicit_permutation_cut_profiles':profiles,'arbitrary_quotient_child_configurations':configurations,'physical_cross_child_factorizations':factorizations,'independent_named_cases':cases,'independent_physical_counts':counts,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
