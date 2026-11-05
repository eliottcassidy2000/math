"""Finite-arity physical construction inventory; no unbounded source search."""
import json,hashlib
from pathlib import Path
from functools import lru_cache
from array import array

EXPECTED=[7,21,67,231,247,533,821,1861]
def need(x,m):
 if not x:raise RuntimeError(m)
@lru_cache(None)
def constants(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=constants(L//2);return g,Q,2*B
 l=(L+1)//2;g,Q,B=constants(l)
 return g+l-l%2,Q+B,2*B

def independent_constants(L):
 n=L-1;k=n.bit_length();B=1<<k;R=sum(((n>>i)&1)<<(k-1-i)for i in range(k))
 prefix=0;f=1
 for b in bin(n)[2:] if n else '':
  if b=='0':f+=prefix+prefix%2
  prefix=2*prefix+int(b)
 return f,B-R,B

def T(v,L,independent=False):
 D,E,U=v;g,Q,B=(independent_constants if independent else constants)(L)
 return g*D-(L%2)*E+L*U,D-int(L==1)*E,Q*D+B*U

def F(v):return v[0]-v[1]+v[2]
def O(v,r=1):return v[0],0,(1<<r)*v[2]
def pair(L,K,kind,independent=False):
 const=independent_constants if independent else constants
 D,E,U,d=const(L)[0],L%2,L,1
 G,e,V,f=const(K)[0],K%2,K,1
 if kind=='N':return (D-E+U)*G-(D-E+d)*e+(D+3*U-d)*V-(U+E)*f
 return (D+2*U)*G-(U+d)*e+(2*D-E+6*U-d)*V-(U+E)*f

def construct_inventory(N,depth=4):
 states=[(1,1,1)];parents=array('i',[-1]);labels=array('i',[-1]);levels=[[0]];seen={states[0]:0};support={}
 for k in range(depth+1):
  nxt=[]
  for i in levels[k]:
   v=states[i];D,E,U=v
   h=F(v)
   if 3<=h<=N:support.setdefault(h,('C1',i,0))
   r=1
   while D+(1<<r)*U<=N:
    support.setdefault(D+(1<<r)*U,('C1',i,r));r+=1
   if k<depth:
    for L in range(1,N//U+1):
     w=T(v,L)
     if F(w)<=N and w not in seen:
      j=len(states);seen[w]=j;states.append(w);parents.append(i);labels.append(L);nxt.append(j)
  if k<depth:levels.append(nxt)
 c1=len(support);tested=0
 for L in range(2,N+1):
  if 3*L*L+1>N:break
  for K in range(L,(N-1)//(3*L)+1):
   tested+=1
   for kind in('N','M'):
    n=pair(L,K,kind)
    if n<=N:support.setdefault(n,(kind,L,K))
 def record(h):
  rec=support[h]
  if rec[0]!='C1':return rec
  _,i,r=rec;ls=[]
  while parents[i]>=0:ls.append(labels[i]);i=parents[i]
  return 'C1',tuple(reversed(ls)),r
 records={h:record(h)for h in sorted(support)}
 return records,dict(normalized_prefixes=len(states),prefix_levels=list(map(len,levels)),C1_values=c1,one_apex_unordered_pairs=tested)

def literal(rows):
 n=len(rows);size=1<<n
 need(all(not rows[v]>>v&1 for v in range(n)),'loops')
 need(all(((rows[u]>>v)&1)+((rows[v]>>u)&1)==1 for u in range(n)for v in range(u+1,n)),'tournament')
 end=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
 for v in range(n):end[1<<v][v]=1
 start[1][0]=1
 for mask in range(1,size):
  for v in range(n):
   if not mask>>v&1:continue
   old=mask^(1<<v)
   if old:
    end[mask][v]=sum(end[old][u]for u in range(n)if old>>u&1 and rows[u]>>v&1)
    start[mask][v]=sum(start[old][u]for u in range(n)if old>>u&1 and rows[u]>>v&1)
 return sum(end[-1]),sum(start[-1][u]for u in range(n)if rows[u]&1)

def main():
 controls=0;N=100001
 for L in range(1,N+1):need(constants(L)==independent_constants(L),'gadget constants')
 # The sink-absorption identity is linear, checked on a basis including signed coordinates.
 for j in range(3):
  v=tuple(int(i==j)for i in range(3))
  for L in range(1,65):
   for r in range(1,9):
    need(T(O(v,r),L)==T(v,(1<<r)*L),'sink absorption');controls+=1
 records,meta=construct_inventory(N)
 missing=[h for h in range(3,N+1,2)if h not in records]
 need(missing==EXPECTED,'restricted family gaps')
 digest=hashlib.sha256();counts={'C1':0,'N':0,'M':0};max_order=0
 for h,rec in records.items():
  counts[rec[0]]+=1
  if rec[0]=='C1':
   _,ls,r=rec;need(len(ls)<=4,'finite arity');v=(1,1,1)
   for L in ls:v=T(v,L,True)
   actual=T(v,1<<r,True)[0]
   order=len(ls)+2+sum((L-1).bit_length()for L in ls)+r
  else:
   actual=pair(rec[1],rec[2],rec[0],True)
   order=(3 if rec[0]=='N'else 4)+(rec[1]-1).bit_length()+(rec[2]-1).bit_length()
  need(actual==h,'independent certificate replay')
  need(order<=h.bit_length()-1+10,'uniform graph order bound');max_order=max(max_order,order)
  digest.update(json.dumps((h,rec),separators=(',',':')).encode());digest.update(b'\n')
 repairs=json.loads((Path(__file__).parent/'REPAIRS.json').read_text())
 need({x['H']for x in repairs}==set(EXPECTED)-{7,21},'repair set')
 for x in repairs:
  need(literal(x['rows'])==(x['H'],x['C']),'literal repair')
  need(len(x['rows'])<=x['H'].bit_length()-1+10,'repair order')
 print(json.dumps(dict(status='PASS',bound=N,gadget_constant_checks=N,restricted_family_missing=missing,covered_with_repairs_except=[7,21],basis_absorption_checks=controls,independent_positive_replays=len(records),certificate_sha256=digest.hexdigest(),selected_certificate_types=counts,literal_repairs=len(repairs),maximum_selected_graph_order=max_order,**meta),indent=2,sort_keys=True))
if __name__=='__main__':main()
