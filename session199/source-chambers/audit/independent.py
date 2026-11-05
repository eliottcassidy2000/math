"""Independent positive-coordinate source bank and chamber inverse audit."""
from pathlib import Path
from functools import lru_cache
from itertools import product
from hashlib import sha256
import heapq,json
ROOT=Path(__file__).resolve().parent.parent
A={'O':((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),
 'Z':((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),
 'J':((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))}
t=(0,0,0,1)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def mv(M,x):return tuple(sum(a*b for a,b in zip(row,x))for row in M)
def raw(v):
 a,f,e,p=v
 return a+2*e+p,e,a+f+2*e+2*p,e+p

def state(w):
 v=t
 for c in w:v=mv(A[c],v)
 return raw(v)

def U(v):return raw(v)[2]
def endpoint(D,V,d):
 if V>=2*D:return 0
 return d if V>D+d else D+d-V

def move(s,c):
 D,E,V,d=s
 if c=='O':return D,0,2*V,d
 if c=='Z':return D-E+V,d,2*V-d,d
 return D-E+V,D-E,D+V,D
cap=2048;bank={t:''};heap=[(2,t)];pred={}
while heap:
 _,v=heapq.heappop(heap)
 for c in A:
  w=mv(A[c],v)
  need(U(v)<U(w),'strict positive U growth')
  if U(w)>cap:continue
  pred.setdefault(raw(w),set()).add((raw(v),c))
  if w not in bank:bank[w]=bank[v]+c;heapq.heappush(heap,(U(w),w))
rawbank={raw(v):w for v,w in bank.items()}
triplebank={(D,V,d)for D,E,V,d in rawbank}
nonbalanced=balanced=0
for v,w in bank.items():
 D,E,V,d=raw(v);a,f,e,p=v
 need(endpoint(D,V,d)==E,'unique endpoint')
 need(f==0 or f>d,'strict f gap')
 need((E==0)==(V>=2*D),'weak separator')
 if E:
  need(E!=V-D,'strict strong separator')
  dd=E if E<V-D else V-D+E
  need(dd==d,'strong deletion recovery')
  if f==0:need(e>p,'strong f-zero separator')
 if V==D:
  balanced+=1
  if (D,d)!=(3,1):
   need(any(S[0]==d and S[2]==(D-d)//2 for S in rawbank),'canonical balanced predecessor')
 elif raw(v)!=(1,0,2,1):
  nonbalanced+=1
  need(len(pred[raw(v)])==1,'one actual letter and predecessor off balance')

def legal(s):
 D,E,V,d=s
 return D>=1 and d>=1 and D%2==d%2==1 and V>=2 and 0<=E<=d and D-E>=d and V>=D+d-E and V>=2*d

def lift(x):
 D,V,d=x
 return D,endpoint(D,V,d),V,d

@lru_cache(None)
def inverse(s):
 if not legal(s):return None
 if s==(1,0,2,1):return ''
 D,E,V,d=s
 candidates=[]
 if D==V:
  if (D,E,d)==(3,1,1):return 'Z'
  h=(D-d)//2
  for dd in range(1,min(d,h//2)+1,2):
   candidates.append((lift((d,h,dd)),'OJ'))
 elif E==0:
  if V%2:return None
  candidates=[(lift((D,V//2,d)),'O')]
 elif E==d:
  if (V+d)%2:return None
  h=(V+d)//2;D0=D-h
  if D0%2==0:D0+=d
  candidates=[(lift((D0,h,d)),'Z')]
 else:
  # Recover old deletion from its known strong (D,E,U), not the target-side selector.
  D0,E0,V0=d,d-E,V-d
  if E0==V0-D0:return None
  d0=E0 if E0<V0-D0 else V0-D0+E0
  candidates=[((D0,E0,V0,d0),'J')]
 for p,tail in candidates:
  if p[2]>=V:continue
  prev=inverse(p)
  if prev is None:continue
  q=p
  for c in tail:q=move(q,c)
  if q==s:return prev+tail
 return None
for s in rawbank:
 w=inverse(s)
 need(w is not None and state(w)==s,'independent reconstructed actual source')
box=guarded=accepted=0
for D in range(1,102,2):
 for V in range(2,161):
  for d in range(1,min(D,V//2)+1,2):
   s=lift((D,V,d));box+=1;guarded+=legal(s)
   w=inverse(s)
   need((w is not None)==(s in rawbank),'complete triple box exact membership')
   if w is not None:need(state(w)==s,'literal recipe at accepted triple');accepted+=1
need(state('OOJO')==(9,0,18,1)and state('ZOJO')==(9,0,18,3),'weak loss boundary')
need(move(state('OOJO'),'J')==move(state('ZOJO'),'J')==(27,9,27,9),'balanced J collision')
for k in range(1,10):
 s=state('Z'*k)
 need(s==(2**k+1,1,2**k+1,1),'zero-J balanced family')
report={'status':'PASS','proof_sha256':sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),
 'candidate_verifier_sha256':sha256((ROOT/'verify.py').read_bytes()).hexdigest(),
 'independent_verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
 'complete_cover_cap':cap,'physical_profiles':len(rawbank),'unique_nonbalanced_predecessors':nonbalanced,
 'balanced_profiles':balanced,'triple_box':{'D_max':101,'U_max':160,'triples':box,'guarded':guarded,'physical':accepted},
 'scope':'Exact nonlinear source language and deterministic nonbalanced inverse PASS; canonical reset does not retain minimum-J cost; no scalar coverage assertion.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
