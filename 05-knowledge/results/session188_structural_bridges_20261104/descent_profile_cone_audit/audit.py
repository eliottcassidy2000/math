"""Independent word-derived kernel and simplex-vertex certificate audit."""
from pathlib import Path
from itertools import product,permutations,combinations,combinations_with_replacement
from functools import lru_cache
import hashlib,json,math,subprocess,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'descent_profile_cone'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def mul(P,Q,n):
 out={}
 for a,x in P.items():
  for b,y in Q.items():
   e=tuple(u+v for u,v in zip(a,b))
   if max(e)<=n:out[e]=out.get(e,0)+x*y
 return {e:c for e,c in out.items() if c}
def path_words(n,reverse):
 def arc(u,v):
  if u==v:return False
  if u<3 and v<3:return ((u-v if reverse else v-u)%3)==1
  if u==3:return v==4 or(v<3 and v!=0)
  if u==4:return v<3 and v!=1
  if v==3:return u==0
  if v==4:return u==1
  return False
 @lru_cache(None)
 def f(counts,mask,last):
  total=sum(counts)+mask.bit_count()
  if last<3:
   if not counts[last]:return 0
   q=list(counts);q[last]-=1;q=tuple(q);m=mask
  else:
   b=1<<(last-3)
   if not mask&b:return 0
   q=counts;m=mask^b
  if total==1:return 1
  return sum(f(q,m,prev) for prev in range(5) if arc(prev,last))
 return {c:sum(f(c,3,last) for last in range(5)) for c in product(range(n+1),repeat=3)}
def word_series(n):
 T=path_words(n,False);U=path_words(n,True);D={e:T[e]-U[e] for e in T};B={};pis=list(permutations(range(3)))
 for e in D:
  value=sum(D[tuple(e[j] for j in pi)] for pi in pis);need(value%2==0,'symmetrized literal-word kernel integral')
  if value:B[e]=value//2
 E={}
 for first in range(3):
  c=[0,0,0]
  for length in range(1,3*n+1):
   c[(first+length-1)%3]+=1
   if max(c)>n:break
   e=tuple(c);E[e]=E.get(e,0)+1
 # Independent all-support check of the previously proved rational kernel.
 K={(2,2,0):1,(2,0,2):1,(0,2,2):1,(2,1,0):1,(2,0,1):1,(1,2,0):1,(1,0,2):1,(0,2,1):1,(0,1,2):1,(2,1,1):-1,(1,2,1):-1,(1,1,2):-1,(1,1,1):-6};reference={}
 for e,c in K.items():
  for j in range(n-max(e)+1):
   a=tuple(x+j for x in e);reference[a]=reference.get(a,0)+(j+1)*c
 need(B=={e:c for e,c in reference.items() if c},'physical-word kernel equals K/(1-XYZ)^2')
 return B,E

def mix(F,qs):return sum(c*qs[0][e[0]]*qs[1][e[1]]*qs[2][e[2]] for e,c in F.items())
def gen(n,d):
 # Build from literal cut subsets of an abstract word whose backward positions
 # are the first d, rather than evaluating the producer binomial expression.
 out=[0]*(n+1)
 for e in {d,n-1-d}:
  mandatory=(1<<e)-1
  for cuts in range(1<<(n-1)):
   if cuts&mandatory==mandatory:out[cuts.bit_count()+1]+=1
 return out

def small():
 records=[]
 for n in range(1,6):
  census=[];pp=list(permutations(range(n)))
  for code in range(1<<(n*(n-1)//2)):
   rows=[0]*n
   for b,(u,v) in enumerate(combinations(range(n),2)):
    x,y=(u,v) if code>>b&1 else(v,u);rows[x]|=1<<y
   c=[0]*n;q=[0]*(n+1)
   for p in pp:
    bad=sum((not(rows[u]>>v&1))<<i for i,(u,v) in enumerate(zip(p,p[1:])))
    c[bad.bit_count()]+=1
    for cuts in range(1<<(n-1)):
     if cuts&bad==bad:q[cuts.bit_count()+1]+=1
   need(c==c[::-1] and sum(c)==math.factorial(n),'word reversal and full mass')
   gs=[gen(n,d) for d in range((n+1)//2)]
   need(q==[sum(c[d]*gs[d][k] for d in range(len(gs))) for k in range(n+1)],'explicit cuts and paired cone')
   census.append([code,c,q])
  records.append({'n':n,'count':len(census),'sha256':hashlib.sha256(json.dumps(census,separators=(',',':')).encode()).hexdigest()})
 return records

def main():
 src=json.loads((SRC/'normal.json').read_text());B,E=word_series(13);records=[];count=0
 for original in src['records']:
  n=original['n'];gs=[gen(n,d) for d in range((n+1)//2)];need(gs==original['generators'],'literal cut generators including central multiplicity')
  F={e:c for e,c in B.items() if max(e)<=n};ee={e:c for e,c in E.items() if max(e)<=n};byh={r['h']:r for r in original['coordinates']};rs=[]
  for h in range(1,3*n+3):
   if h>=2:
    cs=[]
    for t in combinations_with_replacement(range(len(gs)),3):
     mult=len(set(permutations(t)));num=h*mult*mix(F,[gs[i] for i in t]);need(num%3==0,'integral symmetrized cubic coefficient');v=num//3;need(v<=0,'all cone coefficients nonpositive');cs.append([list(t),v]);count+=1
    need(cs==byh[h]['coefficients'],'every frozen cubic coefficient independently matches')
    rs.append({'h':h,'count':len(cs),'minimum':min(v for _,v in cs),'maximum':max(v for _,v in cs),'sha256':hashlib.sha256(json.dumps(cs,separators=(',',':')).encode()).hexdigest()})
   F=mul(F,ee,n)
  need(not F,'higher support empty');records.append({'n':n,'coordinates':rs})
 gs=[gen(13,d) for d in range(7)];boundary=2*mix(mul(B,E,13),[gs[0],gs[0],gs[6]]);need(boundary==1930071780,'order13 coefficient boundary, not sign counterexample');need(count==7434,'complete coefficient count')
 controls=small();need([r['sha256'] for r in controls]==[r['complete_descent_cover_records_sha256'] for r in src['literal_transform_controls']],'all1099 independent word/cut records')
 pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','normal.json','optimized.json']}
 need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer receipts')
 print(json.dumps({'status':'PASS','method':'Independent actual run-word kernel, explicit cut-subset generators and coefficient contraction; no producer imports','coefficients_compared':count,'coordinates':records,'small_controls':controls,'order13_method_boundary':boundary,'source_pins':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
