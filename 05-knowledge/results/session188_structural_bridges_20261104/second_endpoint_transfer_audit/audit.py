"""Independent literal quotient words and predecessor-cover physical controls."""
from pathlib import Path
from itertools import combinations,permutations
from math import comb,factorial
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'second_endpoint_transfer'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(n,code):
 R=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>bit&1 else(j,i);R[u]|=1<<v
 return R
def op(R):return [((1<<len(R))-1)^r^(1<<v) for v,r in enumerate(R)]
def endpoint(Q):
 n=len(Q);S=[0]*n;E=[0]*n;rep=[[0]*n for _ in Q]
 for w in permutations(range(n)):
  if all(Q[u]>>v&1 for u,v in zip(w,w[1:])):S[w[0]]+=1;E[w[-1]]+=1
 for j in range(n):
  for P in combinations(range(n+1),2):
   for v in permutations([i for i in range(n) if i!=j]):
    it=iter(v);w=tuple(j if t in P else next(it) for t in range(n+1))
    if all(Q[u]>>v&1 for u,v in zip(w,w[1:])):rep[j][w[0]]+=1
 return S,E,rep
def refined(Q):
 S,E,R=endpoint(Q);_,_,C=endpoint(op(Q));sigma=[(S[i]-E[i]+2*sum(r[i]-c[i] for r,c in zip(R,C)))%4 for i in range(len(Q))];return S,E,sigma
# Larger named exemplars use an independent memoized remaining-inventory endpoint count.
def refined_large(Q):
 from functools import lru_cache
 def starts(R,k):
  @lru_cache(None)
  def tail(last,left):
   if not any(left):return 1
   z=0
   for v,t in enumerate(left):
    if t and R[last]>>v&1:u=list(left);u[v]-=1;z+=tail(v,tuple(u))
   return z
  out=[]
  for i in range(len(R)):u=list(k);u[i]-=1;out.append(tail(i,tuple(u)))
  return out
 n=len(Q);P=op(Q);S=starts(Q,[1]*n);E=starts(P,[1]*n);R=[];C=[]
 for j in range(n):k=[1]*n;k[j]=2;R.append(starts(Q,k));C.append(starts(P,k))
 return S,E,[(S[i]-E[i]+2*sum(r[i]-c[i] for r,c in zip(R,C)))%4 for i in range(n)]
def tt(n):return [sum(1<<v for v in range(u+1,n)) for u in range(n)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 T=[0]*n
 for i,b in enumerate(B):
  for v,r in enumerate(b):
   T[off[i]+v]=r<<off[i]
   for j,c in enumerate(B):
    if Q[i]>>j&1:T[off[i]+v]|=((1<<len(c))-1)<<off[j]
 return T,off
def augmented(Q,B,i,K=None):
 T,o=blow(Q,B);a=((1<<len(B[i]))-1)<<o[i];c=0
 if K is not None:
  n=len(T);T,z=blow(K,[T]+[[0]]*(len(K)-1));out=sum(1<<z[v] for v in range(1,len(K)) if K[v]&1);a|=out;c|=out
 n=len(T);T+= [1<<(n+1),0]
 for v in range(n):
  for w,mask in((n,a),(n+1,c)):
   if mask>>v&1:T[v]|=1<<w
   else:T[w]|=1<<v
 return tuple(T)
def main():
 raw=(SOURCE/'normal.json').read_bytes();need(raw==(SOURCE/'optimized.json').read_bytes(),'paired producer bytes');j=json.loads(raw)
 for f,h in j['source_pins'].items():need(sha(SOURCE/f)==h,'frozen source pin')
 for f,h in j['dependency_pins'].items():need(sha(ROOT.parent/f)==h,'frozen dependency pin')
 allwords=0
 for n,code,p,ref,*sigma in j['complete_quotient_table']:
  Q=rows(n,code);S,E,sig=refined(Q);p0=sum(((s+e)%2)<<v for v,(s,e) in enumerate(zip(S,E)));ref0=sum((s==2 and not p0>>v&1)<<v for v,s in enumerate(sig));need(p==p0 and ref==ref0 and sigma==sig,'all quotient table entries');allwords+=n
 swaps=0
 for n,N in j['involution_word_counts']:
  count=0
  for j0 in range(n):
   for P in combinations(range(n+2),3):
    if any(b==a+1 for a,b in zip(P,P[1:])):continue
    for other in permutations([v for v in range(n) if v!=j0]):
     it=iter(other);w=tuple(j0 if v in P else next(it) for v in range(n+2));a,b,c=P;new=w[:a+1]+w[b+1:c+1]+w[a+1:b+1]+w[c+1:];need(new!=w and new[0]==w[0] and new[-1]==w[-1] and sorted(zip(w,w[1:]))==sorted(zip(new,new[1:])),'independent return swap');count+=1
  need(count==N==n*comb(n,3)*factorial(n-1),'complete swap count');swaps+=count
 graphs=[];expected=[];index={};occ=0
 def add(T,v):
  nonlocal occ
  occ+=1
  if T in index:need(expected[index[T]]==v,'duplicate graph values')
  else:index[T]=len(graphs);graphs.append(T);expected.append(v)
 types=[tt(2),tt(4),tt(3),[2,4,1]]
 for rec in j['complete_small_physical_controls']:
  Q=rows(rec['order'],rec['quotient_code']);B=[types[t] for t in rec['types']];i=rec['anchor'];S,E,sigma=refined(Q);need(not(S[i]+E[i])%2 and rec['sigma']==sigma[i],'small anchor');need([(a-b)%8 for a,b in zip(*rec['values'])]==[2*sigma[i]%8,0,0],'all cover residues')
  for R,v in zip((Q,op(Q)),rec['values']):add(augmented(R,B,i),v)
 for rec in j['special_controls']:
  Q=rec['quotient'];S,E,sig=refined_large(Q);need([S,E,sig]==[rec['start'],rec['end'],rec['sigma']],'large exact quotient vectors');n=len(Q);need(not any(all(Q[v]&M in(0,M) for v in range(n) if not M>>v&1) for M in range(1,(1<<n)-1) if M.bit_count()>1),'explicit prime modules');need([(a-b)%8 for a,b in zip(*rec['values'])]==[4,0,0],'prime ambient residues')
  for R,v in zip((Q,op(Q)),rec['values']):add(augmented(R,rec['blocks'],rec['anchor'],rec['ambient']),v)
 Q=[238,4,240,118,195,146,34,74]
 for R,v in zip((Q,op(Q)),j['singleton_boundary_values']):add(augmented(R,[[0]]*8,0),v)
 need([v[0] for v in j['singleton_boundary_values']]==[771,2403],'singleton boundary exact')
 with tempfile.TemporaryDirectory(prefix='second-endpoint-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T)) for T in graphs)+'\n';v=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need([v[i:i+3] for i in range(0,len(v),3)]==expected,'all complete actual q1 q2 q3 values')
 return dict(status='PASS',source_sha256={f:sha(SOURCE/f) for f in('PROOF.md','verify.py','physical.cpp','normal.json','optimized.json')},dependency_pins=j['dependency_pins'],checker_sha256=sha(Path(__file__)),counter_sha256=sha(ROOT/'count.cpp'),complete_quotient_cases=len(j['complete_quotient_table']),literal_repeated_label_cases=allwords,explicit_return_swap_words=swaps,complete_small_physical_controls=len(j['complete_small_physical_controls']),special_prime_ambient_controls=len(j['special_controls']),actual_graph_occurrences=occ,distinct_actual_graphs=len(graphs),complete_profile_coordinates_checked=3*occ,singleton_boundary_values=j['singleton_boundary_values'])
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
