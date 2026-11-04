"""Independent word-derived kernel and simplex-vertex certificate audit."""
from pathlib import Path
from itertools import product,permutations,combinations,combinations_with_replacement
from functools import lru_cache
import hashlib,json,math,subprocess,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'equal_ten_ambient'
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

def transitive(n,cap):
 # Ordered surjections, computed without the producer Stirling recurrence.
 return [sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1)) if 1<=k<=n else 0 for k in range(cap+1)]
def mix(F,qs):return sum(c*qs[0][e[0]]*qs[1][e[1]]*qs[2][e[2]] for e,c in F.items())
def finite_swap_controls():
 controls=[]
 for n in range(1,6):
  perms=list(permutations(range(n)));r=n//2;orbit=lambda p:tuple(tuple(sorted(p[2*j:2*j+2])) for j in range(r))+(tuple(p[2*r:]),);need(len({orbit(p) for p in perms})==math.factorial(n)//2**r,'free adjacent-swap orbits');maximum=0;count=0
  for code in range(1<<(n*(n-1)//2)):
   T=[0]*n
   for bit,(u,v) in enumerate(combinations(range(n),2)):
    a,b=(u,v) if code>>bit&1 else(v,u);T[a]|=1<<b
   seen=set();h=0
   for p in perms:
    if all(T[u]>>v&1 for u,v in zip(p,p[1:])):
     key=orbit(p);need(key not in seen,'at most one Hamilton ordering per orbit');seen.add(key);h+=1
   need(h<=math.factorial(n)//2**r,'finite Hamilton bound');maximum=max(maximum,h);count+=1
  controls.append({'order':n,'tournaments':count,'orbit_bound':math.factorial(n)//2**r,'largest_H_in_census':maximum})
 return controls

def complete_maxima(n):
 # Count each permutation in every orientation completing its path arcs.
 # This has no tournament endpoint/subset dynamic program.
 pairs=list(combinations(range(n),2));lookup={e:i for i,e in enumerate(pairs)};allbits=(1<<len(pairs))-1;counts=[0]*(allbits+1)
 for p in permutations(range(n)):
  fixed=ones=0
  for u,v in zip(p,p[1:]):
   bit=1<<lookup[min(u,v),max(u,v)];fixed|=bit
   if u<v:ones|=bit
  free=allbits^fixed;sub=free
  while True:
   counts[ones|sub]+=1
   if not sub:break
   sub=(sub-1)&free
 hist={}
 for h in counts:hist[h]=hist.get(h,0)+1
 maximum=max(counts)
 return {'n':n,'labelled_tournaments':len(counts),'maximum':maximum,'first_maximizer_code':counts.index(maximum),'histogram':[[h,c] for h,c in sorted(hist.items())]},hashlib.sha256(bytes(counts)).hexdigest()

def main():
 paths=[SRC/x for x in('PROOF.md','PROVENANCE.md','verify.py','maxima.cpp')]+[BASE/d/'PROOF.md' for d in('covered_ambient','covered_positivity','profile_cone_obstruction','equal_small_ambient','equal_nine_ambient')]
 pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};producer=json.loads((SRC/'normal.json').read_text());n=10;L=47249
 maxima=[];count_hashes=[]
 for m in(4,6):
  record,digest=complete_maxima(m);maxima.append(record);count_hashes.append({'n':m,'complete_orientation_count_bytes_sha256':digest})
 need(maxima==producer['complete_small_maximum_censuses'],'independent full small-max histograms and attaining codes')
 need([x['maximum'] for x in maxima]==[5,45] and math.comb(10,6)*5*45==L+1,'6+4 bound')
 B,E=word_series(n);top=transitive(n,n);orders=[2,4,6,8,10];vs=[[a+L*b for a,b in zip(top,transitive(s,n))] for s in orders[:-1]]+[top];expected=producer['records'][0]
 need(expected['generator_orders']==orders and expected['lower_weight_budget']==L,'exact order-ten domain');coords={c['h']:c for c in expected['coordinates']};F=B;rows=[];total=0
 for h in range(1,3*n+3):
  if h>=2:
   coefficients=[]
   for t in combinations_with_replacement(range(len(vs)),3):
    mult=math.factorial(3)//math.prod(math.factorial(t.count(i)) for i in set(t));numer=h*mult*mix(F,[vs[i] for i in t]);need(numer%3==0,'independent integral simplex coefficient');value=numer//3;need(value<=0,'complete coefficient sign');coefficients.append([list(t),value]);total+=1
   need(coefficients==coords[h]['coefficients'],'every complete order-ten coefficient matches')
   rows.append({'h':h,'coefficients':len(coefficients),'minimum':min(v for _,v in coefficients),'maximum':max(v for _,v in coefficients),'complete_list_sha256':hashlib.sha256(json.dumps(coefficients,separators=(',',':')).encode()).hexdigest()})
  F=mul(F,E,n)
 need(not F and total==1085 and len(rows)==31,'complete finite support and coefficient totals')
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer normal/optimized replay');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','scope':'Order ten only: complete continuous certificate H<=47250; earlier orders separately audited. No inference above ten.','independent_small_maximum_censuses':maxima,'complete_maximum_count_hashes':count_hashes,'split_bound':L+1,'budget':L,'simplex_vertices':vs,'coordinates':rows,'complete_coefficients_compared':total,'literal_word_kernel_exponent_box':n,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
