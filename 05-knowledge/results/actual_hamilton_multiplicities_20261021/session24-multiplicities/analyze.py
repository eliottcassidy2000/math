import pathlib,json,csv,collections,math,sympy as s
P=pathlib.Path(__file__).parent
D=json.load(open(P/'census.json'))
rows=list(csv.DictReader(open(P/'n8_classes.csv')))
c=collections.defaultdict(collections.Counter);atoms=collections.defaultdict(collections.Counter);centers=collections.defaultdict(collections.Counter);rh=collections.Counter()
for r in rows:
 r={k:int(v) for k,v in r.items()};H=r['H'];a=r['aut'];sc=r['self_converse'];assert H%a==0
 c[H]['ordinary']+=1;c[H]['self_converse']+=sc;c[H]['labelled']+=math.factorial(8)//a;rh[H//a]+=1
 if r['strong']:
  atoms[H][8]+=1
  if sc:centers[H][8]+=1
for H,ct in c.items():
 assert (ct['ordinary']+ct['self_converse'])%2==0
 ct['merged']=(ct['ordinary']+ct['self_converse'])//2
 assert ct['labelled']==D['labelled_n8'][str(H)]
assert len(rows)==6880 and sum(ct['self_converse'] for ct in c.values())==176
assert sum(ct['merged'] for ct in c.values())==3528
assert sum(sum(ct.values()) for ct in atoms.values())==6008
assert sum(h*k for h,k in rh.items())==2**21
D['counts_by_order']['8']={H:dict(ct) for H,ct in sorted(c.items())}
for key,new in [('strong_atom_inventory_through_order8',atoms),('self_converse_strong_inventory_through_order8',centers)]:
 old=D[key.replace('order8','order7')]
 for h,by_n in old.items():
  for n,num in by_n.items():new[int(h)][int(n)]+=num
 D[key]={h:dict(sorted(by_n.items())) for h,by_n in sorted(new.items())}
D['h_counts_n8']={h:k for h,k in sorted(rh.items())}
D['n8_validation']={'ordinary':6880,'self_converse':176,'merged':3528,'strong':6008,'labelled_histogram_matches_independent_frozen_census':'all 320 entries','sum_labelled':2**28,'sum_H_over_Aut':2**21}
(P/'extended_census.json').write_text(json.dumps(D,indent=2,sort_keys=True)+'\n')
z,n=s.symbols('z n');A={int(h):sum(k*z**int(m) for m,k in byn.items()) for h,byn in D['strong_atom_inventory_through_order8'].items() if int(h)<75};C={int(h):sum(k*z**int(m) for m,k in byn.items()) for h,byn in D['self_converse_strong_inventory_through_order8'].items() if int(h)<75};C[1]=1+z
from functools import lru_cache
@lru_cache(None)
def words(H):
 if H==1:return ((0,0,1),)
 out=collections.Counter()
 for q,pol in A.items():
  if H%q:continue
  for (m,),coef in s.Poly(pol,z).terms():
   for r,l,k in words(H//q):out[r+1,l+m]+=int(coef)*k
 return tuple((r,l,k) for (r,l),k in sorted(out.items()))
@lru_cache(None)
def F(H):return s.sympify(sum(k*z**m/(1-z)**(r+1) for r,m,k in words(H)))
def J(H):return sum(F(x).subs(z,z*z)*C.get(H//(x*x),0) for x in range(1,math.isqrt(H)+1) if H%(x*x)==0)
def pcount(H,N):return sum(k*math.comb(N-m+r,r) for r,m,k in words(H) if N>=m)
analyses={}
for H in range(1,75,2):
 f=s.cancel(F(H));j=s.cancel(J(H));mrg=s.cancel((f+j)/2)
 poly=s.expand(sum(k*s.prod(n-m+i for i in range(1,r+1))/math.factorial(r) for r,m,k in words(H)))
 start=max((m for r,m,k in words(H)),default=0)
 for N in range(3,9):
  ct=D['counts_by_order'][str(N)].get(str(H),D['counts_by_order'][str(N)].get(H,{}))
  assert pcount(H,N)==ct.get('ordinary',0),(H,N,pcount(H,N),ct)
  coef=int(s.expand(s.series(j,z,0,N+1).removeO()).coeff(z,N))
  assert coef==ct.get('self_converse',0),(H,N,coef,ct)
 for N in range(max(1,start),max(1,start)+5):assert int(poly.subs(n,N))==pcount(H,N)
 analyses[H]={'strong_atoms':str(A.get(H,0)),'strong_centers':str(C.get(H,0)),'ordinary_gf':str(f),'self_converse_gf':str(j),'merged_gf':str(mrg),'eventual_ordinary_polynomial':str(s.factor(poly)),'polynomial_valid_from_safe_bound':start,'nontrivial_scc_max':max((r for r,m,k in words(H)),default=None),'ordinary_denominator':str(s.factor(s.denom(f))),'merged_denominator':str(s.factor(s.denom(mrg))),'self_converse_denominator':str(s.factor(s.denom(j))),'core_word_terms':words(H)}
for H,r in analyses.items():
 B=0
 while 5**(B+1)<=H**3:B+=1
 r['universal_Busch_ordinary_onset']=B
 r['universal_Busch_parity_onset']=1+B
 r['polynomial_valid_from_core_term_bound']=max((m-r0 for r0,m,k in r['core_word_terms']),default=0)
with open(P/'multiplicities_n3_n8.csv','w') as out:
 w=csv.writer(out);w.writerow(['n','actual_H','ordinary_classes','self_converse_classes','converse_merged_classes','labelled_tournaments'])
 for n0 in range(3,9):
  hist=D['counts_by_order'][str(n0)]
  for h0 in sorted({int(x) for x in hist}|{7,21}):
   c0=hist.get(str(h0),hist.get(h0,{}))
   w.writerow([n0,h0,c0.get('ordinary',0),c0.get('self_converse',0),c0.get('merged',0),c0.get('labelled',0)])
(P/'low_H_formulas.json').write_text(json.dumps(analyses,indent=2,sort_keys=True)+'\n')
print('n8 checks:',D['n8_validation'])
for H in [3,5,9,15,19,21,25,27,33,35,39,41,49,63]:
 a=analyses[H];print(H,a['eventual_ordinary_polynomial'],'J=',a['self_converse_gf'],'A=',a['strong_atoms'],'C=',a['strong_centers'])
