"""Literal full-response interpolation and independent ordered-cover contrasts."""
import itertools,json,math
from pathlib import Path

def need(ok,why):
 if not ok: raise RuntimeError(why)
def endpoint(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not mask>>v&1:continue
   rem=mask^(1<<v)
   if not rem:continue
   d[mask][v]=sum(d[rem][u] for u in range(n) if rem>>u&1 and a[u]>>v&1)
 return d
def converse(a):
 n=len(a);return [((1<<n)-1)^row^(1<<u) for u,row in enumerate(a)]
def hp(a):return sum(endpoint(a)[-1])
def covers(a):
 n=len(a);e=endpoint(a);s=endpoint(converse(a));full=(1<<n)-1
 return [[0 if u==v else sum(e[b][u]*s[full^b][v] for b in range(1,full) if b>>u&1 and not b>>v&1) for v in range(n)] for u in range(n)]
def extend(a,b,c):
 n=len(a);o=list(a)+[1<<(n+1),0]
 for i in range(n):
  if b>>i&1:o[i]|=1<<n
  else:o[n]|=1<<i
  if c>>i&1:o[i]|=1<<(n+1)
  else:o[n+1]|=1<<i
 return o
def inflate(a,r):
 b=[row&31 for row in a[:5]]+[0]*r
 for i in range(5):
  for v in range(5,5+r):
   if a[i]>>5&1:b[i]|=1<<v
   else:b[v]|=1<<i
 for u in range(5,5+r):
  for v in range(u+1,5+r):b[u]|=1<<v
 return b
def reverse_pair(a,u,v):
 b=list(a);b[u]^=1<<v;b[v]^=1<<u;return b
def polynomial_difference(a,b):
 n=len(a);small=[m for m in range(1<<n) if m.bit_count()<=2]
 vals={(x,y):hp(extend(a,x,y))-hp(extend(b,x,y)) for x in small for y in small}
 coeff=[]
 for x in small:
  for y in small:
   v=0;xx=x
   while True:
    yy=y
    while True:
     v+=(-1)**((x^xx).bit_count()+(y^yy).bit_count())*vals[xx,yy]
     if not yy:break
     yy=(yy-1)&y
    if not xx:break
    xx=(xx-1)&x
   if v:coeff.append([x,y,v])
 # Non-interpolation assignments verify use of actual named masks, not sorted invariants.
 for x,y in [(2**n-1,2**n-1),(2**n-1,0),(0,2**n-1),(0b1010101&((1<<n)-1),0b1101010&((1<<n)-1))]:
  got=sum(v for i,j,v in coeff if i&x==i and j&y==j)
  need(got==hp(extend(a,x,y))-hp(extend(b,x,y)),'full interpolation control')
 return {'basis':'a_mask,c_mask,integer_coefficient; all omitted coefficients zero; actual names unchanged','query_count':len(vals),'nonzero':coeff,'degree':max((x.bit_count()+y.bit_count() for x,y,v in coeff),default=0),'first_query_witness':next(([x,y,v] for (x,y),v in vals.items() if v),None)}
def run():
 seeds=[(16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3)]
 out=[]
 for j,q in enumerate(seeds):
  a=inflate(q,2);b=reverse_pair(a,5,6);C=covers(a);D=covers(b);Q=covers(q)
  signs=[]
  for w in range(5):
   L=C[6][w]-C[5][w]+C[w][5]-C[w][6]
   R=D[6][w]-D[5][w]+D[w][5]-D[w][6]
   need(L==Q[5][w]+Q[w][5]>0,'twin sign formula');need(R==-L,'reversed twin sign')
   signs.append([w,L,R])
  need(hp(a)==4*3**2+9*2**2-3+2*j,'gate count')
  ambient=polynomial_difference(a,b)
  cubic=polynomial_difference(a,converse(a))
  need(cubic['degree']<=3,'converse cubic degree')
  need(all(not((x==0 and y.bit_count()==2) or (y==0 and x.bit_count()==2)) for x,y,v in cubic['nonzero']),'pure quadratic zero')
  out.append({'state':j,'rows':a,'module_reversed_rows':b,'global_converse_rows':converse(a),'H':hp(a),'full_named_C':C,'reversed_full_named_C':D,'exterior_L_signs':signs,'full_named_module_reversal_polynomial':ambient,'full_named_global_converse_cubic':cubic})
 return {'status':'all explicit checks pass; all-order exclusion inherited from C2','cases':out}
if __name__=='__main__':
 result=run();print(json.dumps(result,indent=2,sort_keys=True))
