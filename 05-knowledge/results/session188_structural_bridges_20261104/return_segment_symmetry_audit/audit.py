"""Independent forward-inventory and explicit orbit audit; no producer imports."""
from pathlib import Path
from itertools import combinations,product
from collections import defaultdict,Counter
from math import factorial
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'return_segment_symmetry'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def graph(Q,a,c):
 n=len(Q);T=list(Q)+[1<<(n+1),0]
 for u in range(n):
  for v,mask in((n,a),(n+1,c)):
   if mask>>u&1:T[u]|=1<<v
   else:T[v]|=1<<u
 return tuple(T)
cache={}
def count(T,k):
 key=T,tuple(k)
 if key in cache:return cache[key]
 cap=tuple(k)+(1,1);D={(tuple(0 for _ in cap),-1):1}
 for step in range(sum(cap)):
  new=defaultdict(int)
  for (used,last),z in D.items():
   for v,limit in enumerate(cap):
    if used[v]==limit or(last>=0 and not T[last]>>v&1):continue
    w=list(used);w[v]+=1;new[tuple(w),v]+=z
  D=new
 ans=sum(D.values());cache[key]=ans;return ans
def words(T,k):
 cap=tuple(k)+(1,1);stack=[((),cap)];ans=[]
 while stack:
  w,left=stack.pop()
  if not any(left):ans.append(w);continue
  for v,z in enumerate(left):
   if z and(not w or T[w[-1]]>>v&1):u=list(left);u[v]-=1;stack.append((w+(v,),tuple(u)))
 return ans
def returns(Q,k,i,R):
 cap=[z//R for z in k];cap[i]=0;stack=[((),tuple(cap))];out=[]
 while stack:
  w,left=stack.pop()
  if w and Q[w[-1]]>>i&1:out.append(w)
  for v,z in enumerate(left):
   if z and Q[i if not w else w[-1]]>>v&1:u=list(left);u[v]-=1;stack.append((w+(v,),tuple(u)))
 return sorted(out)
def contraction(Q,a,c,k,i):
 R=k[i]-1;T=graph(Q,a,c);terms=[]
 for W in returns(Q,k,i,R):
  reduced=list(k);reduced[i]=1
  for v in W:reduced[v]-=R
  terms.append([list(W),reduced,count(T,reduced)])
 lhs=count(T,k);rhs=sum(r[2] for r in terms);need((lhs-rhs)%2==0,'independent contraction')
 heavy=[v for v,z in enumerate(k) if v!=i and z>=R];tri=any(Q[i]>>u&1 and Q[u]>>v&1 and Q[v]>>i&1 for u in heavy for v in heavy);need(tri or lhs%2==0,'independent heavy triangle obstruction');return lhs,terms,tri

def main():
 nraw=(SOURCE/'normal.json').read_bytes();need(nraw==(SOURCE/'optimized.json').read_bytes(),'producer paired bytes');j=json.loads(nraw)
 for f,h in j['source_pins'].items():need(sha(SOURCE/f)==h,'producer source binding')
 need(sha(ROOT.parent/'coarse_mod8_affine'/'PROOF.md')==j['dependency_pin'],'inherited source binding')
 for Q,a,c,k,i,value,terms,tri in j['complete_parity_records']:
  v,t,b=contraction(Q,a,c,k,i);need(v==value and sorted(t)==sorted(terms) and b==tri,'all complete contraction records')
 total_words=total_orbits=0
 for Q,a,c,k,size,hist in j['orbit_controls']:
  ww=words(graph(Q,a,c),k);need(len(ww)==size==count(graph(Q,a,c),k),'independent explicit words');groups=defaultdict(set)
  for w in ww:
   at=[p for p,v in enumerate(w) if v==0];P=w[:at[0]];S=w[at[-1]+1:];L=tuple(w[a+1:b] for a,b in zip(at,at[1:]));need(all(L),'nonempty returns');groups[P,tuple(sorted(L)),S].add(w)
  hh=Counter()
  for (P,L,S),orbit in groups.items():
   den=1
   for v in Counter(L).values():den*=factorial(v)
   expected=factorial(len(L))//den;need(len(orbit)==expected,'full orbit size');need((len(orbit)%2==1)==(len(set(L))==1),'power-of-two odd orbit');hh[len(orbit)]+=1
  need(sorted(hh.items())==[tuple(x) for x in hist],'complete orbit histogram');total_words+=len(ww);total_orbits+=len(groups)
 for Q,a,c,k,v in j['factorial_controls']:
  need(count(graph(Q,a,c),k)==v and v%factorial(k[0]-1)==0,'all free factorial controls')
 for Q,k,v in j['sharp_factorial']:need(count(graph(Q,0,(1<<len(Q))-1),k)==v==factorial(k[0]-1),'sharp factorial')
 for Q,a,c,k,i,v,terms,tri in j['sharp_odd']:
  w,t,b=contraction(Q,a,c,k,i);need(w==v==1 and sorted(t)==sorted(terms) and b==tri,'sharp odd example')
 Q,a,c,k,i,v,terms,tri=j['triangle_not_sufficient'];need(tri and count(graph(Q,a,c),k)==v and v%2==0,'necessary-not-sufficient example')
 # Independent four-child controls outside the producer parity domain.
 extra=[]
 for code in range(64):
  Q=[0]*4
  for b,(u,v) in enumerate(combinations(range(4),2)):
   p,q=(u,v) if code>>b&1 else(v,u);Q[p]|=1<<q
  for a,c in((0,15),(3,5),(6,9),(15,0)):
   k=[3,2,2,1];v,t,tri=contraction(Q,a,c,k,0);extra.append([Q,a,c,k,v,t,tri])
 # Non-power-of-two return count: a mixed orbit can be odd (3!/(2!1!)=3).
 need(factorial(3)//(factorial(2)*factorial(1))==3,'power-of-two restriction')
 # Fresh producer replays are compared as full bytes, not only parsed summaries.
 for flags in([],['-O']):need(subprocess.check_output([sys.executable]+flags+[str(SOURCE/'verify.py')])==nraw,'fresh complete producer replay')
 return dict(status='PASS',source_sha256={f:sha(SOURCE/f) for f in('PROOF.md','verify.py','normal.json','optimized.json')},inherited_sha256=sha(ROOT.parent/'coarse_mod8_affine'/'PROOF.md'),checker_sha256=sha(Path(__file__)),complete_parity_cases=len(j['complete_parity_records']),complete_factorial_cases=len(j['factorial_controls']),explicit_word_sets=len(j['orbit_controls']),explicit_words=total_words,explicit_orbits=total_orbits,additional_four_child_controls=extra,distinct_forward_inventory_counts=len(cache),full_producer_normal_optimized_replays=True)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
