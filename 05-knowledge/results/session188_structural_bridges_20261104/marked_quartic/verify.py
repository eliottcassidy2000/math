"""Endpoint-marked child profiles and arbitrary-quotient shared quartic kernels."""
from pathlib import Path
import importlib.util,itertools,json,math,hashlib
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('covered',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def induced(rows,names):return [sum(1<<j for j,w in enumerate(names) if rows[u]>>w&1) for u in names]
def marked(a):
 n=len(a);N=1<<n;end=v.endpoints(a);start=v.endpoints(v.opposite(a));H=[sum(x) for x in end];H[0]=1;D=[[0]*(n+1) for _ in range(N)];D[0][0]=1
 for mask in range(1,N):
  sub=mask
  while sub:
   for k in range(1,n+1):D[mask][k]+=H[sub]*D[mask^sub][k-1]
   sub=(sub-1)&mask
 E=[[0]*(n+1) for _ in range(n)];S=[[0]*(n+1) for _ in range(n)];full=N-1
 for mask in range(1,N):
  for u in range(n):
   if mask>>u&1:
    for k in range(1,n+1):E[u][k]+=end[mask][u]*D[full^mask][k-1];S[u][k]+=start[mask][u]*D[full^mask][k-1]
 for k in range(1,n+1):v.need(sum(r[k] for r in E)==sum(r[k] for r in S)==D[full][k],'marked profile normalization')
 if n>=2:
  for u in range(n):v.need(E[u][n-1]-S[u][n-1]==math.factorial(n-2)*(n-1-2*a[u].bit_count()),'top marked imbalance score formula')
 return D[full],E,S
def add(*ps):
 d={}
 for p in ps:
  for e,c in p.items():d[e]=d.get(e,0)+c
 return {e:c for e,c in d.items() if c}
def multiply(p,q,caps):
 d={}
 for a,x in p.items():
  for b,y in q.items():
   c=tuple(u+w for u,w in zip(a,b))
   if all(u<=w for u,w in zip(c,caps)):d[c]=d.get(c,0)+x*y
 return {c:x for c,x in d.items() if x}
def run_matrix(Q,caps):
 q=len(Q);counts=sorted((c for c in itertools.product(*[range(m+1) for m in caps]) if any(c)),key=lambda c:(sum(c),c));G=[[{} for _ in Q] for _ in Q]
 for first in range(q):
  D={}
  for c in counts:
   z=[0]*q
   if sum(c)==1:
    if c[first]==1:z[first]=1
   else:
    for last in range(q):
     if not c[last]:continue
     p=list(c);p[last]-=1;before=D.get(tuple(p),[0]*q);z[last]=sum(before[u] for u in range(q) if Q[u]>>last&1)
   D[c]=z
   for last,x in enumerate(z):
    if x:G[first][last][c]=x
 return G
def kernel(Q,caps,i,j):
 G=run_matrix(Q,caps);col=lambda j:add(*[G[u][j] for u in range(len(Q))]);row=lambda i:add(*G[i]);A=multiply(col(i),row(j),caps);B=multiply(col(j),row(i),caps)
 return add(A,{e:-c for e,c in B.items()})
def evaluate(K,profiles):
 total=0
 for e,c in K.items():
  z=c
  for r,q in zip(e,profiles):z*=q[r]
  total+=z
 return total
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,r in enumerate(A):
   T[off[i]+u]=r<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return T,off
def csym(a,u,w):
 E=v.endpoints(a);S=v.endpoints(v.opposite(a));full=(1<<len(a))-1;value=0
 for mask in range(1,full):
  if mask>>u&1 and not mask>>w&1:value+=E[mask][u]*S[full^mask][w]
  if mask>>w&1 and not mask>>u&1:value+=E[mask][w]*S[full^mask][u]
 return value
def case(Q,blocks,deleted,first,second):
 T,off=blow(Q,blocks);U,_=blow(v.opposite(Q),blocks);h,localw=deleted;i,localu=first;j,localv=second;v.need(i!=j,'different marked children')
 names=[list(range(len(A))) for A in blocks];names[h].remove(localw);cards=[induced(A,ns) for A,ns in zip(blocks,names)];profiles=[marked(A) for A in cards];caps=list(map(len,cards));K=kernel(Q,caps,i,j);u=names[i].index(localu);w=names[j].index(localv)
 p=[x[0] for x in profiles];p[i]=profiles[i][1][u];p[j]=profiles[j][2][w];one=evaluate(K,p)
 p=[x[0] for x in profiles];p[i]=profiles[i][2][u];p[j]=profiles[j][1][w];predicted=one-evaluate(K,p)
 W=off[h]+localw;I=off[i]+localu;J=off[j]+localv;remain=[x for x in range(len(T)) if x!=W];tc=induced(T,remain);uc=induced(U,remain);actual=csym(tc,remain.index(I),remain.index(J))-csym(uc,remain.index(I),remain.index(J));v.need(actual==predicted,'independent full-card symmetric-cover quartic')
 return {'quotient':Q,'child_rows':blocks,'deleted_and_marks':[W,I,J],'a_mask':(1<<W)|(1<<I),'c_mask':(1<<W)|(1<<J),'T':T,'U':U,'kernel_terms':[[*e,c] for e,c in sorted(K.items())],'child_card_profiles':[x[0] for x in profiles],'first_end_profile':profiles[i][1][u],'first_start_profile':profiles[i][2][u],'second_end_profile':profiles[j][1][w],'second_start_profile':profiles[j][2][w],'shared_quartic_difference':actual}
if __name__=='__main__':
 count=0
 for n in range(1,6):
  for c in range(1<<(n*(n-1)//2)):marked(v.tournament(n,c));count+=1
 Q3=[2,4,1];Q5=[sum(1<<((i+j)%5) for j in (1,2)) for i in range(5)];C=[2,4,1]
 cases=[case(Q3,[C]*3,(0,0),(0,1),(1,0)),case(Q3,[C]*3,(2,0),(0,1),(1,0)),case(Q5,[C]*5,(0,0),(0,1),(1,0)),case(Q5,[C]*5,(2,0),(0,1),(1,0))]
 full=json.loads((ROOT.parent/'covered_response'/'full_response.json').read_text());lookup={(r[0],r[1]):r[5]-r[6] for r in full['coordinate_data']};v.need(cases[0]['shared_quartic_difference']==lookup[3,9]==-24,'pair-free nine-name retained coefficient');localized=0;shared=0;disjoint=0
 for (a,c),z in lookup.items():
  if a.bit_count()!=2 or c.bit_count()!=2 or not z:continue
  if not a&c:disjoint+=1;continue
  if (a&c).bit_count()!=1:continue
  w=(a&c).bit_length()-1;u=(a^(a&c)).bit_length()-1;t=(c^(a&c)).bit_length()-1;shared+=1
  if u//3!=t//3:v.need(w//3 in (u//3,t//3),'balanced child support localization');localized+=1
 v.need((shared,disjoint,localized)==(216,324,216),'complete retained quartic support counts')
 literal=json.loads((ROOT/'literal.json').read_text());v.need(literal['status']=='PASS' and len(literal['cases'])==len(cases),'complete physical interpolation controls');literal_count=0
 for case,record in zip(cases,literal['cases']):
  a,c=case['a_mask'],case['c_mask'];v.need(record['query']==[a,c] and len(record['values'])==16,'matching physical quartic query');delta=0
  for x,y,f,g in record['values']:
   v.need(x&a==x and y&c==y and 0<=f<=math.factorial(len(case['T'])+2) and 0<=g<=math.factorial(len(case['U'])+2),'exact literal count bounds');delta+=(-1)**((a^x).bit_count()+(c^y).bit_count())*(f-g)
   if len(case['T'])<=9:v.need([v.hp(v.extend(rows,x,y)) for rows in (case['T'],case['U'])]==[f,g],'independent Python physical query control')
   if x==y==0:v.need([v.hp(rows) for rows in (case['T'],case['U'])]==[f,g],'all-source port count equals old Hamilton count')
   literal_count+=2
  v.need(delta==record['coefficient_difference']==case['shared_quartic_difference'],'physical finite difference versus marked kernel')

 print(json.dumps({'status':'PASS','scope':'All-order arbitrary-quotient cross-child shared-quartic factorization and endpoint-balanced support localization; no general separation or collision theorem.','marked_tournament_profiles_checked':count,'literal_physical_counts_retained':literal_count,'exact_cases':cases,'pair_free_nine_full_quartic_support':{'shared':shared,'disjoint':disjoint,'localized_cross_child_shared':localized},'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT/'literal.cpp',ROOT/'literal.json',ROOT.parent/'covered_response'/'verify.py',ROOT.parent/'response'/'PROOF.md']}},sort_keys=True,indent=2))
