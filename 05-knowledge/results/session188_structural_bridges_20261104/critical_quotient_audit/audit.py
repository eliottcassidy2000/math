"""Primary-convention and five-vertex classification controls, without producer imports."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'critical_quotient_synthesis'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def graph(n,arc):return tuple(sum(1<<v for v in range(n) if u!=v and arc(u,v)) for u in range(n))
def family(r,kind):
 n=2*r+1
 if kind==1:return graph(n,lambda u,v:0<(v-u)%n<=r)
 if kind==2:return graph(n,lambda u,v:((0<(v-u)%n<=r) != (u<r and v<r)))
 return graph(n,lambda u,v: u<v if u<2*r and v<2*r else (v%2==0 if u==2*r else u%2==1))
def source(r,kind):
 if kind==1:return family(r,1),list(range(2*r+1))
 if kind==2:
  names=[('a',j) for j in range(r+1)]+[('b',j) for j in range(1,r+1)]
  arcs=set()
  for i in range(1,r+1):
   for j in range(i+1,r+1):arcs.add((('b',i),('b',j)))
   for j in range(r+1):
    arcs.add((('b',i),('a',j)) if i<=j else (('a',j),('b',i)))
  for i in range(r+1):
   for j in range(i):arcs.add((('a',i),('a',j)))
  return graph(len(names),lambda u,v:(names[u],names[v]) in arcs),[2*r-j if k=='a' else r-j for k,j in names]
 names=[('b',0)]+[('a',j) for j in range(1,2*r+1)]
 arcs={(('a',i),('a',j)) for i in range(1,2*r+1) for j in range(i+1,2*r+1)}
 for j in range(1,2*r+1):arcs.add((('b',0),('a',j)) if j%2 else (('a',j),('b',0)))
 return graph(len(names),lambda u,v:(names[u],names[v]) in arcs),[2*r]+list(range(2*r))
def prime(T):
 n=len(T)
 return not any(all(T[u]&M in(0,M) for u in range(n) if not M>>u&1) for M in range(1,(1<<n)-1) if M.bit_count()>1)
def code(T):return sum((bool(T[u]>>v&1))<<j for j,(u,v) in enumerate(combinations(range(len(T)),2)))
def main():
 paths=[SRC/'PROOF.md',SRC/'verify.py']+[BASE/d/'PROOF.md' for d in ('ordered_module_skew','alternating_prime_quotient','staircase_prime_quotient','marked_pair_quotient')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};maps=[]
 for r in range(2,20):
  for k in (1,2,3):
   T,f=source(r,k);Q=family(r,k);n=len(Q);need(sorted(f)==list(range(n)),'bijective physical convention map');need(all(bool(T[u]>>v&1)==bool(Q[f[u]]>>f[v]&1) for u in range(n) for v in range(n)),'all primary edge classes under exact map');maps.append({'r':r,'kind':k,'map':f})
 classes=[];all_codes=set()
 for k in (1,2,3):
  Q=family(2,k);cs={code(graph(5,lambda u,v:Q[p[u]]>>p[v]&1)) for p in permutations(range(5))};need(not(all_codes&cs),'pairwise nonisomorphic three primary families')
  all_codes|=cs;classes.append({'kind':k,'labelled_count':len(cs),'codes':sorted(cs)})
 pairs=list(combinations(range(5),2));actual=set()
 for c in range(1<<len(pairs)):
  T=[0]*5
  for bit,(u,v) in enumerate(pairs):
   a,b=(u,v) if c>>bit&1 else (v,u);T[a]|=1<<b
  if prime(T):actual.add(c)
 need(actual==all_codes and len(actual)==264,'complete labelled prime-five base case')
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer replay')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'frozen sources')
 print(json.dumps({'status':'PASS','primary_url':'https://trotter.math.gatech.edu/papers/83.pdf','rendered_pages_reviewed':[198,203],'independent_convention_maps':maps,'labelled_prime_five_classes':classes,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
