"""Independent direct-permutation audit plus unchanged producer replays."""
import itertools,json,hashlib,subprocess,sys
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'module_decomposition'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,mask):
 R=[0]*n
 for k,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:R[u]|=1<<v
  else:R[v]|=1<<u
 return tuple(R)
def mod(R,M):return all((R[v]&M) in (0,M) for v in range(len(R)) if not M>>v&1)
def modules(R):return [S for S in range(1,1<<len(R)) if mod(R,S)]
def flip(R,F):
 U=list(R)
 for u,v in itertools.combinations(range(len(R)),2):
  if sum(bool(M>>u&1 and M>>v&1) for M in F)%2:U[u]^=1<<v;U[v]^=1<<u
 return tuple(U)
@lru_cache(None)
def tables(R):
 n=len(R);N=(1<<n)-1;H=[0]*(N+1);H[0]=1
 for S in range(1,N+1):
  names=[v for v in range(n) if S>>v&1]
  H[S]=sum(all(R[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(names))
 q=[[0]*(n+1) for _ in range(N+1)];q[0][0]=1
 for S in range(1,N+1):
  part=S
  while part:
   rest=S^part
   for k in range(1,n+1):q[S][k]+=H[part]*q[rest][k-1]
   part=(part-1)&S
 C={t for t in itertools.combinations(range(n),3) if H[sum(1<<v for v in t)]==3}
 return H,tuple(x for a in q for x in a),C

def Q(R,children):
 reps=[(M&-M).bit_length()-1 for M in children]
 return tuple(sum(1<<j for j,v in enumerate(reps) if R[u]>>v&1) for u in reps)
def transitive(R):return sorted(row.bit_count() for row in R)==list(range(len(R)))
def opp(R):return tuple(((1<<len(R))-1)^row^(1<<i) for i,row in enumerate(R))

def one(R,counts):
 n=len(R);Ms=modules(R);nontriv=[M for M in Ms if M.bit_count()>1];strong=[M for M in Ms if all(not(M&N) or M&N in (M,N) for N in Ms)];nodes={}
 for M in strong:
  if M.bit_count()<2:continue
  subs=[N for N in strong if N!=M and N&M==N]
  ch=[N for N in subs if not any(N!=Z and N&Z==N for Z in subs)];nodes[M]=(ch,Q(R,ch))
 H,q,C=tables(R);counts['sources']+=1
 for bits in range(1<<len(nontriv)):
  F=[M for i,M in enumerate(nontriv) if bits>>i&1];U=flip(R,F);uH,uq,uC=tables(U);bad=False
  for M in strong:need(mod(U,M),'persistent original strong block')
  for M,(ch,oldQ) in nodes.items():
   newQ=Q(U,ch)
   if transitive(oldQ):bad|=not transitive(newQ)
   else:need(newQ in (oldQ,opp(oldQ)),'prime whole-converse action')
  need(C<=uC,'old triangles survive')
  need(all(a<=b for a,b in zip(H,uH)),'direct-permutation hereditary Hamilton domination')
  need(all(a<=b for a,b in zip(q,uq)),'independent cover domination')
  eq=H[-1]==uH[-1];need(eq==(not bad)==(H==uH)==(q==uq)==(C==uC),'equality equivalences')
  counts['families']+=1;counts['scalar_equal']+=eq;counts['scalar_strict']+=not eq
  counts['H_coordinates']+=len(H);counts['cover_coordinates']+=len(q)

def run():
 source_hashes={p:hashlib.sha256((SOURCE/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py')}
 counts={'sources':0,'families':0,'scalar_equal':0,'scalar_strict':0,'H_coordinates':0,'cover_coordinates':0}
 for n in range(1,5):
  for mask in range(1<<(n*(n-1)//2)):one(rows(n,mask),counts)
 # All interval-XOR targets from TT5, plus two mixed hierarchical controls.
 one(tuple(sum(1<<j for j in range(i+1,5)) for i in range(5)),counts)
 # Explicit constructors avoid guessed rows.
 def substitute(Q,blocks):
  offsets=[0]
  for B in blocks:offsets.append(offsets[-1]+len(B))
  R=[0]*offsets[-1]
  for i,B in enumerate(blocks):
   for u,row in enumerate(B):
    R[offsets[i]+u]=row<<offsets[i]
    for j,C in enumerate(blocks):
     if Q[i]>>j&1:R[offsets[i]+u]|=((1<<len(C))-1)<<offsets[j]
  return tuple(R)
 one(substitute((2,4,1),((2,0),(2,0),(0,))),counts)
 one(substitute((2,4,1),((2,4,1),(2,0),(0,))),counts)
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized,'producer normal optimized agreement')
 (HERE/'producer.normal.stdout').write_bytes(normal);(HERE/'producer.optimized.stdout').write_bytes(optimized)
 need(source_hashes=={p:hashlib.sha256((SOURCE/p).read_bytes()).hexdigest() for p in source_hashes},'source stayed frozen during audit')
 result={'status':'PASS','independent_counts':counts,'independent_universe':'Every labelled source through order 4; all 1024 interval-XOR families from TT5; C3[TT2,TT2,1]; C3[C3,TT2,1]. H from literal permutations.','producer_result':json.loads(normal),'source_sha256':{p:hashlib.sha256((SOURCE/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py')},'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
