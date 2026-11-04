"""Independent word-derived kernel and simplex-vertex certificate audit."""
from pathlib import Path
from itertools import product,permutations,combinations,combinations_with_replacement
from functools import lru_cache
import hashlib,json,math,subprocess,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'descent_amgm_ambient'
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

def main():
 from fractions import Fraction
 src=json.loads((SRC/'normal.json').read_text());B,E=word_series(16);records=[];total=0;used=[]
 for original in src['records']:
  n=original['n'];gs=[gen(n,d) for d in range((n+1)//2)];need(gs==original['generators'],'independent cut generators');F={e:c for e,c in B.items() if max(e)<=n};ee={e:c for e,c in E.items() if max(e)<=n};byh={r['h']:r for r in original['coordinates']};out=[]
  for h in range(1,3*n+3):
   if h>=2:
    cs={}
    for t in combinations_with_replacement(range(len(gs)),3):
     num=h*len(set(permutations(t)))*mix(F,[gs[i] for i in t]);need(num%3==0,'integer coefficients');cs[t]=num//3;total+=1
    need([[list(t),v] for t,v in cs.items()]==byh[h]['coefficients'],'complete independently reconstructed cubic')
    positive={t:v for t,v in cs.items() if v>0};rr={'h':h,'count':len(cs),'coefficients_sha256':hashlib.sha256(json.dumps([[list(t),v] for t,v in cs.items()],separators=(',',':')).encode()).hexdigest()}
    if positive:
     cert=byh[h]['amgm'];d=cert['scales'];need(len(d)==len(gs) and all(type(x)==int and x>0 for x in d),'positive integer scales');A=[-cs[i,i,i] for i in range(len(gs))];need(all(x>=0 for x in A),'no positive pure cube')
     # Independent rational loads directly from weighted AM-GM.
     load=[sum(Fraction(v*t.count(i)*math.prod(d[j] for j in t),3*d[i]**3) for t,v in positive.items()) for i in range(len(gs))]
     need(all(x<=a for x,a in zip(load,A)),'exact rational cube domination')
     capacity=[3*A[i]*d[i]**3 for i in range(len(gs))];usage=[int(load[i]*3*d[i]**3) for i in range(len(gs))];slack=[c-u for c,u in zip(capacity,usage)]
     need(A==cert['negative_cube_magnitudes'] and capacity==cert['capacity'] and usage==cert['usage'] and slack==cert['slack'],'all retained integer loads capacities slacks')
     need([[list(t),v] for t,v in positive.items()]==cert['positive_terms'],'all positive monomials allocated exactly')
     rr['scales']=d;rr['nonnegative_slack']=slack;rr['positive_terms']=len(positive);used.append([n,h])
    else:need(byh[h]['all_coefficients_nonpositive'],'negative-coefficient coordinate')
    out.append(rr)
   F=mul(F,ee,n)
  need(not F,'all higher support zero');records.append({'n':n,'coordinates':out})
 need(total==src['coefficient_count'] and used==src['positive_coefficient_coordinates'] and len(used)==9,'all coordinates and exactly nine allocated cubics')
 need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired modes')
 pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','normal.json','optimized.json']}
 print(json.dumps({'status':'PASS','source_pins':pins,'coefficients_compared':total,'coordinates':records,'positive_coordinates':used,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Entire nonnegative palindromic descent cone orders13..16, with exact weighted AM-GM. No conclusion above16.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
