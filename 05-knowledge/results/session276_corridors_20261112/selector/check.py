from pathlib import Path
import sys,json,itertools
sys.path.insert(0,str(Path(__file__).resolve().parent))
# Load helper definitions without exploratory work at import.
ns={"__file__":str(Path(__file__).parent/"probe.py")};exec((Path(__file__).parent/'probe.py').read_text().split('for n in (3,5,7,9):')[0],ns)
hp,profile,apex,cyclic,covers=[ns[x] for x in ('hp','profile','apex','cyclic','covers')]

def substitute(B,v,G):
 outer=[i for i in range(len(B)) if i!=v];n=len(G); rows=list(G)+[0]*len(outer)
 for i,a in enumerate(outer,n):
  for j,b in enumerate(outer,n):
   if B[a]>>b&1:rows[i]|=1<<j
  for j in range(n):
   if B[a]>>v&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
 return tuple(rows)
checks=[]
for n,sets,want in [(5,({0,1},{0,2}),[41,43]),(7,({0,2},{0,1}),[559,561])]:
 ts=[apex(cyclic(n),a) for a in sets]; ps=[profile(t) for t in ts]; cs=[covers(t) for t in ts]
 assert [p['H'] for p in ps]==want
 assert {k:v for k,v in ps[0].items() if k!='H'}=={k:v for k,v in ps[1].items() if k!='H'}
 delta=[b-a for a,b in zip(*cs)]
 outer=(2,4,1); actual=[hp(substitute(outer,0,t))[-1] for t in ts]
 assert actual==[3*c[1]+2*c[2] for c in cs]
 fan=[]
 for m in range(4):
  # outer vertex0 -> transitive1..m -> vertex m+1 -> vertex0.
  B=[0]*(m+2)
  B[0]=sum(1<<i for i in range(1,m+1));B[m+1]=1
  for i in range(1,m+1):B[i]=sum(1<<j for j in range(i+1,m+2))
  vals=[hp(substitute(B,0,t))[-1] for t in ts]
  predicted=[((1<<m)+1)*c[1]+2*((1<<m)-1)*c[2] for c in cs]
  # m0 outer is w->v: one-run coefficient1, not2.
  if m==0:predicted=[c[1] for c in cs]
  assert vals==predicted
  fan.append(dict(m=m,counts=vals))
 checks.append(dict(fan=fan,order=n+1,graphs=ts,profiles=ps,covers=cs,cover_difference=delta,triangle_outputs=actual))
# Complete labelled outer contexts through order3; independent literal DP.
count=0
for k in range(1,4):
 edges=list(itertools.combinations(range(k),2))
 for word in range(1<<len(edges)):
  B=[0]*k
  for i,(a,b) in enumerate(edges):
   if word>>i&1:B[a]|=1<<b
   else:B[b]|=1<<a
  for v in range(k):
   ts=[apex(cyclic(5),a) for a in ({0,1},{0,2})]
   values=[hp(substitute(B,v,t))[-1] for t in ts];d=values[1]-values[0];h=hp(tuple(B))[-1]
   assert d>=2*h and (d==2)==(h==1);count+=1
out=dict(status='PASS',literal_sources=checks,complete_outer_context_checks=count)
print(json.dumps(out,indent=2))
