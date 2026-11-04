"""Exact SCC-local response factorization; no external modules or assertions."""
import itertools,json
from pathlib import Path

def need(x,m):
 if not x:raise RuntimeError(m)
def hp(a):
 n=len(a);full=(1<<n)-1;d=[[0]*n for _ in range(full+1)]
 for u in range(n):d[1<<u][u]=1
 for m in range(1,full+1):
  for u in range(n):
   t=d[m][u]
   if not t:continue
   avail=a[u]&~m
   while avail:
    bit=avail&-avail;avail-=bit;v=bit.bit_length()-1
    d[m|bit][v]+=t
 return sum(d[-1])
def induced(a,vs):return [sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs]
def conv(a):return [((1<<len(a))-1)^r^(1<<u) for u,r in enumerate(a)]
def tournament(n,code):
 a=[0]*n
 for k,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1:a[u]|=1<<v
  else:a[v]|=1<<u
 return a
def sccs(a):
 n=len(a);reach=[r|(1<<u) for u,r in enumerate(a)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 remain=set(range(n));out=[]
 while remain:
  u=min(remain);s=sorted(v for v in remain if reach[u]>>v&1 and reach[v]>>u&1);out.append(s);remain.difference_update(s)
 parts=list(out)
 out.sort(key=lambda s:-sum(a[s[0]]>>t[0]&1 for t in parts if t!=s))
 return out
def substitute(q,m):
 # root 0 gets physical names 0..len(m)-1; old q vertex i>0 becomes len(m)+i-1.
 n=len(m);out=list(m)+[0]*(len(q)-1)
 for i in range(len(q)):
  for j in range(len(q)):
   if not(q[i]>>j&1):continue
   I=range(n) if i==0 else [n+i-1];J=range(n) if j==0 else [n+j-1]
   for u in I:
    for v in J:out[u]|=1<<v
 return out
def ext(a,x,y):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for u in range(n):
  if x>>u&1:b[u]|=1<<n
  else:b[n]|=1<<u
  if y>>u&1:b[u]|=1<<(n+1)
  else:b[n+1]|=1<<u
 return b
def setup(q,m):
 parts=sccs(q);K=next(s for s in parts if 0 in s);qK=induced(q,K);a=substitute(q,m);b=substitute(qK,m)
 keep=list(range(len(m)))+[len(m)+v-1 for v in K if v!=0]
 factor=1
 for s in parts:
  if s!=K:factor*=hp(induced(q,s))
 incoming=sum(1<<(len(m)+v-1) for v in range(1,len(q)) if v not in K and q[v]&1)
 def lift(mask):return incoming|sum(1<<keep[u] for u in range(len(b)) if mask>>u&1)
 return a,b,factor,lift,parts

def run():
 records=[];queries=0
 # Every labelled rooted quotient of order<=4, fixed physical root0.
 # Relabelling covers arbitrary roots. Two nonisomorphic cores cover a
 # reversal-invariant and reversal-sensitive scalar profile independently.
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   q=tournament(n,code)
   for m in ([2,4,1],[6,4,0]):
    a,b,f,lift,parts=setup(q,m);basis=[s for s in range(1<<len(b)) if s.bit_count()<=2]
    # Full bidegree-(2,2) basis, including quartics; no cubic truncation.
    for x in basis:
     for y in basis:
      need(hp(ext(a,lift(x),lift(y)))==f*hp(ext(b,x,y)),'SCC factorization')
      queries+=1
    for x,y in [(2**len(b)-1,2**len(b)-1),(2**len(b)-1,0),(0,2**len(b)-1)]:
     need(hp(ext(a,lift(x),lift(y)))==f*hp(ext(b,x,y)),'high mask factorization');queries+=1
    records.append({'quotient_rows':q,'core_rows':m,'sccs':parts,'active_rows':b,'ambient_rows':a,'factor':f,'basis_queries':len(basis)**2})
 return {'status':'PASS','labelled_rooted_quotients':75,'core_choices':2,'cases':len(records),'exact_query_checks':queries,'scope':'Complete named bidegree-(2,2) interpolation per listed active block plus three high-mask controls; not a reconstruction census.','records':records}
if __name__=='__main__':print(json.dumps(run(),sort_keys=True,indent=2))
