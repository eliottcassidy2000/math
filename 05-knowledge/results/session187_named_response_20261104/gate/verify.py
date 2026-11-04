"""Standard-library, explicit-gate exact verification; no removable assertions."""
import itertools,json,math
from pathlib import Path
from base import SEEDS,inflate,reverse,packet,need
from literal import hp,extend

def endpoints(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if not m>>v&1:continue
   rem=m^(1<<v)
   if rem:d[m][v]=sum(d[rem][u] for u in range(n) if rem>>u&1 and a[u]>>v&1)
 return d[-1]

def skeleton(a,root,ignored=None,deleted=(),boundary=None):
 ext=[v for v in range(len(a)) if v!=root and v not in deleted];cs=[0]*(len(ext)+2)
 for p in itertools.permutations(ext):
  blocked=None
  if boundary:
   left,right=boundary;k=p.index(left)
   if k+1==len(p) or p[k+1]!=right:continue
   blocked=k+1
  good=[];bad=[]
  for g in range(len(p)+1):
   if g==blocked:continue
   u=None if not g else p[g-1];v=None if g==len(p) else p[g]
   if (u is None or u==ignored or a[u]>>root&1) and (v is None or root==ignored or a[root]>>v&1):good.append(g)
   if u is not None and v is not None and u!=ignored and not a[u]>>v&1:bad.append(g)
  if not set(bad)<=set(good):continue
  for k in range(len(good)+1):
   if 0<=k-len(bad)<=len(good)-len(bad):cs[k]+=math.comb(len(good)-len(bad),k-len(bad))
 return cs

B=(44,42,42);C=(-5,-1,3);D=(-4,-3,-3)
def score(j,r,k):return 18*3**(r-k)*4**(k-1)+B[j]*2**(r-k)*3**(k-1)+C[j]*2**(k-1)+D[j]*(k==r)
def global_gap(j,r):return [2*(4*3**r+6*2**r-6),2*(10*3**r+12*2**r-9),2*(10*3**r-3*2**r-9)][j]
def coefficient(t,u,am,cm):
 z=0;a=am
 while True:
  c=cm
  while True:
   z+=(-1)**((am^a).bit_count()+(cm^c).bit_count())*(hp(extend(t,a,c))-hp(extend(u,a,c)))
   if c==0:break
   c=(c-1)&cm
  if a==0:break
  a=(a-1)&am
 return z

def run():
 profiles=[];rows=[];literal_checks=0
 for j,s in enumerate(SEEDS):
  a=inflate(s,2);first=skeleton(a,6,5)
  expected=[(53,137,98,18),(56,137,96,18),(60,141,96,18)][j]
  need(tuple(first[:4])==expected and not any(first[4:]),('first skeleton',j,first))
  p=skeleton(a,6,deleted=(0,),boundary=(5,3));q=skeleton(a,6,deleted=(0,),boundary=(3,5))
  diff=[x-y for x,y in zip(p,q)]
  need(diff[:3]==[2,4,1] and not any(diff[3:]),('quartic skeleton',j,diff))
  left=skeleton(s,5,0);right=skeleton(reverse(s),5,0)
  need(tuple(left[:4])==[(15,45,39,8),(27,75,63,14),(18,54,51,14)][j] and not any(left[4:]),('original exterior profile',j,left))
  need(tuple(right[:4])==[(9,27,21,4),(12,30,21,4),(18,39,24,4)][j] and not any(right[4:]),('op exterior profile',j,right))
  need(tuple(p[:3])==[(3,6,2),(6,11,4),(6,11,4)][j] and not any(p[3:]),('forward quartic profile',j,p))
  need(tuple(q[:3])==[(1,2,1),(4,7,3),(4,7,3)][j] and not any(q[3:]),('backward quartic profile',j,q))
  profiles.append({'state':j,'first_score_profile':first,'exterior0_score_profile':left,'op_exterior0_score_profile':right,'quartic_forward_profile':p,'quartic_backward_profile':q,'quartic_difference':diff})
  for r in range(1,8):
   t=inflate(s,r);h,c=packet(t);e=endpoints(t);scores=[e[v]+sum(c[v]) for v in range(5,5+r)]
   expected=[score(j,r,k) for k in range(1,r+1)]
   need(scores==expected,('all scores',j,r,scores,expected))
   need(all(scores[k+1]-scores[k]==score(j,r-1,k+1) for k in range(r-1)),('adjacent recurrence',j,r))
   opp=reverse(t);hh,cc=packet(opp);ee=endpoints(opp)
   need(hh==h and 2*(e[0]+sum(c[0])-ee[0]-sum(cc[0]))==global_gap(j,r)>0,('global gap',j,r))
   # Literal physical-port checks independent of endpoint-cover derivation.
   for v in {0,5,4+r}:
    need(hp(extend(t,1<<v,0))==h+2*(e[v]+sum(c[v])),('literal score',j,r,v));literal_checks+=1
   quartic=None
   if r>=2:
    u=list(t);u[5]^=1<<6;u[6]^=1<<5
    quartic=coefficient(t,u,(1<<0)|(1<<3),(1<<0)|(1<<6));literal_checks+=32
    need(quartic==2**(r-2)+2-(r==2),('quartic',j,r,quartic))
   rows.append({'state':j,'r':r,'H':h,'module_scores':scores,'global_query_gap':global_gap(j,r),'pair_reversal_quartic':quartic})
 # Complete emitted finite polynomial data, with checks not trusting its status string.
 raw=json.loads(Path(__file__).with_name('full_response.json').read_text());cases=raw['cases'];battery_queries=0;coefficient_count=0
 need(len(cases)==39,'full table scope')
 expected_cases={(j,r,'global_converse') for j in range(3) for r in range(1,8)}|{(j,r,'pair_reversal') for j in range(3) for r in range(2,8)}
 need({(x['state'],x['r'],x['kind']) for x in cases}==expected_cases,'complete unique comparison universe')
 for case in cases:
  r=case['r'];j=case['state'];t=inflate(SEEDS[j],r);u=reverse(t)
  if case['kind']=='pair_reversal':u=list(t);u[5]^=1<<6;u[6]^=1<<5
  need(case['rows']==t and case['other_rows']==u,'physical rows')
  n=len(t);ms=[m for m in range(1<<n) if m.bit_count()<=2];need(ms==case['battery_masks'],'complete battery')
  vals={(a,c):pair[0]-pair[1] for (a,c),pair in zip(itertools.product(ms,repeat=2),case['battery'])}
  need(len(case['battery'])==len(ms)**2,'battery length');cs={}
  for a,c in vals:
   z=0;aa=a
   while True:
    cc=c
    while True:
     z+=(-1)**((a^aa).bit_count()+(c^cc).bit_count())*vals[aa,cc]
     if cc==0:break
     cc=(cc-1)&c
    if aa==0:break
    aa=(aa-1)&a
   if z:cs[a,c]=z
  need(cs=={(a,c):v for a,c,v in case['coefficients']},'all coefficients')
  degree=max((a.bit_count()+c.bit_count() for a,c in cs),default=0)
  need(degree==case['degree']==(3 if case['kind']=='global_converse' else 4),'exact degree')
  if case['kind']=='global_converse':
   need(all(not ((a==0 or c==0) and a.bit_count()+c.bit_count()==2) for a,c in cs),'pure quadratic')
   need(vals[1,0]==global_gap(j,r),'all-order fixed query')
  else:need(cs[9,65]==2**(r-2)+2-(r==2),'retained quartic')
  for a,c,expected in case['dense_controls']:
   need(sum(v for (x,y),v in cs.items() if x&a==x and y&c==y)==expected,'dense reconstruction')
  # One nonuniform higher-weight control independently recomputed per case.
  a,c,expected=case['dense_controls'][-1]
  need(hp(extend(t,a,c))-hp(extend(u,a,c))==expected,'independent dense physical DP');literal_checks+=2
  battery_queries+=2*len(vals);coefficient_count+=len(cs)
 return {'status':'PASS','profiles':profiles,'all_order_formula_controls':rows,'full_response_cases':len(cases),'retained_battery_counts':battery_queries,'retained_nonzero_coefficients':coefficient_count,'independent_literal_counts':literal_checks,'scope':'All 3 states r=1..7; complete named global-converse data and pair-reversal data for r>=2. All-order claims separately proved.'}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
