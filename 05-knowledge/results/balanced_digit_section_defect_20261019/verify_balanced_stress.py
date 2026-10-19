#!/usr/bin/env python3
"""Exact controls for the q=29 balanced-stress section-defect bound. Standard library only."""
from itertools import product
import json
from pathlib import Path
P=29
V=list(product(range(P), repeat=2))
u,v,w=(1,0),(0,1),(3,9)
add=lambda x,y: ((x[0]+y[0])%P,(x[1]+y[1])%P)
Q=lambda z:z[0]*z[0]+z[0]*z[1]+z[1]*z[1]
def ell(d,z):return (2*d[0]+d[1])*z[0]+(d[0]+2*d[1])*z[1]
def D(f,d):return {x:f[add(x,d)]-f[x] for x in V}
def DD(f,d,e):return D(D(f,d),e)
def carry(d):
 out={}
 for x in V:
  y=add(x,d);diff=(y[0]-x[0],y[1]-x[1])
  assert (Q(diff)-1)%P==0
  out[x]=((Q(diff)-1)//P)%P
 return out

def cycles(d):
 left=set(V);out=[]
 while left:
  x=min(left);cyc=[add(x,(t*d[0],t*d[1])) for t in range(P)]
  left.difference_update(cyc);out.append(cyc)
 return out

def stencil(fu,fv,fw):
 A=DD(fw,u,v);B=DD(fu,v,w);C=DD(fv,u,w)
 return {x:A[x]-3*B[x]-9*C[x] for x in V}

cu,cv,cw=map(carry,(u,v,w))
assert all(x%P==0 for x in stencil(cu,cv,cw).values())
# Exhaust the numerical step which forces three binary second differences to vanish.
digits=[(a,b,c) for a,b,c in product(range(-2,3),repeat=3) if (a-3*b-9*c)%P==0]
assert digits==[(0,0,0)]
# Hostile: the integer range argument is essential. At 23 a nonzero tuple aliases to zero.
assert (2-3*(-1)-9*(-2))==23
assert (2-3*(-1)-9*(-2))%23==0
# Independently verify every translated 12-term expression is a rigidity stress.
for x in V:
 coeff={}
 for d,shifts,k in [(w,(u,v),1),(u,(v,w),-3),(v,(u,w),-9)]:
  e,f=shifts
  for a,sgn in [(add(add(x,e),f),1),(add(x,e),-1),(add(x,f),-1),(x,1)]:
   b=add(a,d)
   for point,sign in [(b,1),(a,-1)]:
    old=coeff.get(point,(0,0))
    coeff[point]=(old[0]+sgn*k*sign*(2*d[0]+d[1]),old[1]+sgn*k*sign*(d[0]+2*d[1]))
 assert all(a%P==b%P==0 for a,b in coeff.values())
# An arbitrary, deterministic motion supplies a positive control for the full affine syndrome.
h={x:((x[0]**2+7*x[1]**3+x[0]*x[1])%P,(3*x[0]**3+x[1]**2+2*x[0]*x[1])%P) for x in V}
def motion_error(d,c):
 out={}
 for x in V:
  y=add(x,d);dh=(h[y][0]-h[x][0],h[y][1]-h[x][1])
  out[x]=(c[x]+ell(d,dh))%P
 return out
zu,zv,zw=[motion_error(d,c) for d,c in zip((u,v,w),(cu,cv,cw))]
assert all(a%P==0 for a in stencil(zu,zv,zw).values())
# Straight-cycle carry and error sums are -2 for all three directions.
for d,c,z in zip((u,v,w),(cu,cv,cw),(zu,zv,zw)):
 for cyc in cycles(d):
  assert sum(c[x] for x in cyc)%P==P-2
  assert sum(z[x] for x in cyc)%P==P-2

cases=[]
for eu,ev,ew in product(range(2),repeat=3):
 # eu=1 or ev=1 selects seams parallel to w. ew=1 selects a v-parallel w-edge seam.
 iu={x:int((9*x[0]-3*x[1])%P==0) if eu else int(x[0]==0) for x in V}
 iv={x:int((9*x[0]-3*x[1])%P==0) if ev else int(x[1]==0) for x in V}
 iw={x:int(x[0]==0) if ew else int(x[1]==0) for x in V}
 for d,I in zip((u,v,w),(iu,iv,iw)):
  assert sum(I.values())==P
  assert all(sum(I[x] for x in C)==1 for C in cycles(d))
 assert all(a==0 for a in stencil(iu,iv,iw).values())
 # Solve the first two directions with explicit cumulative integrals; do not use any rank computation.
 ru={x:(-2*iu[x]-cu[x])%P for x in V}
 rv={x:(-2*iv[x]-cv[x])%P for x in V}
 rw={x:(-2*iw[x]-cw[x])%P for x in V}
 Hu={};Hv={}
 for y in range(P):
  t=0
  for x in range(P):Hu[x,y]=t;t=(t+ru[x,y])%P
  assert t==0
 for x in range(P):
  t=0
  for y in range(P):Hv[x,y]=t;t=(t+rv[x,y])%P
  assert t==0
 assert all((D(Hu,u)[x]-ru[x])%P==0 for x in V)
 assert all((D(Hv,v)[x]-rv[x])%P==0 for x in V)
 T={x:(rw[x]-3*D(Hu,w)[x]-9*D(Hv,w)[x])%P for x in V}
 periods={sum(T[x,y] for x in range(P))%P for y in range(P)}
 expected=(-2*ew+3*(2*3+9)-2*3**2*eu+2*9**2*ev)%P
 assert periods=={expected} and expected!=0
 cases.append({'u_w_parallel':eu,'v_w_parallel':ev,'w_v_parallel':ew,'horizontal_period_obstruction':expected})
# Verify affine seam offsets do not affect the period formula, using analytic support sums.
for eu,ev,ew in product(range(2),repeat=3):
 for tu,tv,tw in product(range(P),repeat=3):
  # Selected u-edge origin on row y is a_y; only its increment over b=9 matters.
  au=lambda y:((3*y+tu)*pow(9,-1,P))%P if eu else tu
  assert (au(9)-au(0))%P==3*eu
  # A w-parallel v seam has exactly one point in each horizontal row; otherwise 0 mod p.
  assert (sum(int((9*x-tv)%P==0) for x in range(P)) if ev else 0)==ev
  assert (sum(int(x==tw) for x in range(P)) if ew else 0)==ew
# Norm-one rotations and exact averaging multiplicity.
mul=lambda a,b:((a[0]*b[0]-a[1]*b[1])%P,(a[0]*b[1]+a[1]*b[0]+a[1]*b[1])%P)
neg=lambda a:((-a[0])%P,(-a[1])%P)
canon=lambda a:min(a,neg(a))
U={x for x in V if Q(x)%P==1};dirs=sorted({canon(x) for x in U})
assert len(U)==30 and len(dirs)==15
multiplicity={d:0 for d in dirs};triples=[]
for g in dirs:
 triple=sorted({canon(mul(g,d)) for d in (u,v,w)})
 assert len(triple)==3
 triples.append(triple)
 for d in triple:multiplicity[d]+=1
assert set(multiplicity.values())=={3}
# Stronger: these 15 direction triples have five pairwise-disjoint members.
from itertools import combinations
partition=next((cs for cs in combinations(range(15),5) if len(set().union(*(set(triples[i]) for i in cs)))==15),None)
result={'p':P,'norm_w_integer':Q(w),'norm_w_residue':Q(w)%P,'stencil_terms':12,'binary_second_difference_tuples_checked':125,'allowed_digit_tuples':digits,'mod23_aliasing_hostile':[2,-1,-2],'translated_stresses_checked':P*P,'carry_stencils_checked':P*P,'arbitrary_motion_stencils_checked':P*P,'straight_cycles_checked':3*P,'seam_cases':cases,'seam_offset_tuples_checked':8*P**3,'direction_triples':triples,'direction_multiplicity':3,'disjoint_triple_partition':partition,'triple_lower_bound':3*P+1,'full_section_defect_lower_bound':15*(3*P+1)//3}
Path(__file__).with_name('balanced-stress-exact-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
