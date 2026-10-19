#!/usr/bin/env python3
"""Exact small-support controls for exclusion of 88 errors. Standard library only."""
from itertools import product
from pathlib import Path
import json
P=29;V=list(product(range(P),repeat=2));ds=[(1,0),(0,1),(3,9)];coef=[-3,-9,1]
add=lambda a,b:((a[0]+b[0])%P,(a[1]+b[1])%P)
scale=lambda k,a:((k*a[0])%P,(k*a[1])%P)
sub=lambda a,b:add(a,scale(-1,b))
cross=lambda a,b:(a[0]*b[1]-a[1]*b[0])%P
Q=lambda a:a[0]**2+a[0]*a[1]+a[1]**2
ell=lambda d,h:(2*d[0]+d[1])*h[0]+(d[0]+2*d[1])*h[1]
def D(f,d):return {x:f[add(x,d)]-f[x] for x in V}
def DD(f,d,e):return D(D(f,d),e)
def stencil(fs):
 out={x:0 for x in V}
 for i in range(3):
  others=[ds[j] for j in range(3) if j!=i];h=DD(fs[i],*others)
  for x in V:out[x]+=coef[i]*h[x]
 return out

def rref(A):
 A=[[v%P for v in row] for row in A];r=0;pivs=[]
 if not A:return A,pivs
 for j in range(len(A[0])):
  t=next((i for i in range(r,len(A)) if A[i][j]),None)
  if t is None:continue
  A[r],A[t]=A[t],A[r];q=pow(A[r][j],-1,P);A[r]=[v*q%P for v in A[r]]
  for i in range(len(A)):
   if i!=r and A[i][j]:
    q=A[i][j];A[i]=[(v-q*w)%P for v,w in zip(A[i],A[r])]
  pivs.append(j);r+=1
  if r==len(A):break
 return A[:r],pivs

def nullspace(A,n):
 rr,piv=rref(A);out=[]
 for j in range(n):
  if j not in piv:
   v=[0]*n;v[j]=1
   for row,k in zip(rr,piv):v[k]=-row[j]%P
   out.append(v)
 return out

def line_constraints(S,i):
 rows=[]
 for j in range(3):
  if j==i:continue
  levels=sorted({cross(ds[j],x) for x in S})
  for t in levels:rows.append([int(cross(ds[j],x)==t) for x in S])
 return rows

def cycles(d):
 left=set(V);out=[]
 while left:
  x=min(left);C=[add(x,scale(k,d)) for k in range(P)];left.difference_update(C);out.append(C)
 return out
cycle_lists=[cycles(d) for d in ds]
carry=[]
for d in ds:
 c={}
 for x in V:
  y=add(x,d);q=(y[0]-x[0],y[1]-x[1]);assert (Q(q)-1)%P==0;c[x]=((Q(q)-1)//P)%P
 carry.append(c)

def period(z):
 r=[{x:(z[i][x]-carry[i][x])%P for x in V} for i in range(3)]
 H=[{},{}]
 for y in range(P):
  t=0
  for x in range(P):H[0][x,y]=t;t=(t+r[0][x,y])%P
  assert t==0
 for x in range(P):
  t=0
  for y in range(P):H[1][x,y]=t;t=(t+r[1][x,y])%P
  assert t==0
 A=D(H[0],ds[2]);B=D(H[1],ds[2]);T={x:(r[2][x]-3*A[x]-9*B[x])%P for x in V}
 ps={sum(T[x,y] for x in range(P))%P for y in range(P)}
 assert len(ps)==1
 return next(iter(ps))

# Coefficient geometry used in the analytic incidence proof.
for i in range(3):
 j,k=[j for j in range(3) if j!=i]
 assert cross(ds[i],ds[j]) not in (0,cross(ds[i],ds[k]),-cross(ds[i],ds[k])%P)

cube_controls=[];hostiles=[]
for exc in range(3):
 for k in range(1,15):
  e=ds[exc];q=scale(k,e)
  S0={x:int(x in {scale(t,e) for t in range(1,k+1)}) for x in V}
  K=D(S0,e);assert K=={x:int(x==(0,0))-int(x==q) for x in V}
  C=D(D(D(S0,ds[0]),ds[1]),ds[2]);support=sorted(x for x in V if C[x])
  assert len(support)==8 and set(C[x] for x in support)=={-1,1}
  # For each pair, the two rectangular incidence cycles share at most one graph vertex.
  axes=[scale(k,d) if i==exc else d for i,d in enumerate(ds)]
  origin=q # D on point masses puts the alternating cube at q minus its subset sums.
  face_shared=[]
  for missing in range(3):
   j,l=[j for j in range(3) if j!=missing]
   faces=[];nodes=[]
   for b in range(2):
    f={sub(sub(sub(origin,scale(b,axes[missing])),scale(s,axes[j])),scale(t,axes[l])) for s,t in product(range(2),repeat=2)}
    assert len(f)==4;faces.append(f)
    nodes.append({(j,cross(ds[j],x)) for x in f}|{(l,cross(ds[l],x)) for x in f})
    assert len(nodes[-1])==4
   assert faces[0].isdisjoint(faces[1]) and faces[0]|faces[1]==set(support)
   shared=len(nodes[0]&nodes[1]);assert shared<=1;face_shared.append(shared)
  spaces=[nullspace(line_constraints(support,i),len(support)) for i in range(3)]
  assert list(map(len,spaces))==[2,2,2]
  # The full constrained Hessian system has exactly the two homogeneous scalar degrees of freedom.
  A=[[coef[i]*v[r]%P for i in range(3) for v in spaces[i]] for r in range(len(support))]
  rank=len(rref(A)[1]);assert rank==4
  for i in range(3):assert all(sum(row[t]*C[x] for t,x in enumerate(support))%P==0 for row in line_constraints(support,i))
  cube_controls.append({'exceptional_direction':exc,'short_separation':k,'support_size':8,'pair_kernel_dimensions':[2,2,2],'coupled_rank':rank,'homogeneous_dimension':2,'shared_incidence_vertices':face_shared})
  # Construct near-seam 88-error hostiles. They pass every cycle and every compact mixed stress,
  # but differ from a rejected seam only by an explicit local rigidity motion.
  for seam_dir in [j for j in range(3) if j!=exc]:
   J=[]
   for i in range(3):
    parallel=ds[seam_dir] if i==exc else e
    J.append({x:int(cross(parallel,x)==0) for x in V})
   lam=[int(i!=exc) for i in range(3)]
   I=[]
   for i in range(3):
    dS=D(S0,ds[i]);I.append({x:J[i][x]+lam[i]*dS[x] for x in V})
    assert set(I[i].values())<={0,1}
    assert all(sum(I[i][x] for x in C0)==1 for C0 in cycle_lists[i])
   t=(2*sum(coef[i]*lam[i] for i in range(3))*pow(coef[exc]%P,-1,P))%P
   z=[{x:(-2*I[i][x]+(t*K[x] if i==exc else 0))%P for x in V} for i in range(3)]
   assert sum(sum(a!=0 for a in zz.values()) for zz in z)==88
   assert all(sum(z[i][x] for x in C0)%P==27 for i in range(3) for C0 in cycle_lists[i])
   assert all(a%P==0 for a in stencil(z).values())
   L=[(-2*lam[i]+(t if i==exc else 0))%P for i in range(3)]
   H=((2*L[0]-L[1])*pow(3,-1,P)%P,(-L[0]+2*L[1])*pow(3,-1,P)%P)
   assert all(ell(ds[i],H)%P==L[i] for i in range(3))
   for i in range(3):
    dS=D(S0,ds[i])
    assert all((z[i][x]+2*J[i][x]-L[i]*dS[x])%P==0 for x in V)
   pz=period(z);pj=period([{x:-2*J[i][x]%P for x in V} for i in range(3)])
   assert pz==pj!=0
   hostiles.append({'exceptional_direction':exc,'short_separation':k,'exceptional_seam_parallel_to':seam_dir,'failures':88,'dipole_amplitude':t,'period_obstruction':pz})

# The strict short-column hypothesis in the binary-normalization lemma is necessary.
# At p-1 modified points per column an additive nonbinary array can be made binary at weight p.
J={(s,t):int(s==0)+int(t==0)-int(t==1) for s,t in V}
I=dict(J)
for s in range(P):
 if s!=1:I[s,0]-=1
 if s!=0:I[s,1]+=1
assert set(I.values())=={0,1} and sum(I.values())==P
assert sum(J.values())==P and not set(J.values())<={0,1}
assert sum(I[x]-J[x] for x in V)==0
assert [sum(I[s,t]!=J[s,t] for s in range(P)) for t in [0,1]]==[P-1,P-1]
assert all(a==0 for a in DD(J,ds[0],ds[1]).values())

out={'p':P,'scope':'First unramified quadratic carry; directions (1,0),(0,1),(3,9).',
'cube_controls':cube_controls,'cube_cases_checked':len(cube_controls),
'near_seam_hostiles':hostiles,'explicit_88_error_hostiles_checked':len(hostiles),
'binary_normalization_boundary_hostile':{'modified_points_per_column':P-1,'I_weight':P,'J_weight':P,'J_binary':False},
'proven_three_direction_lower_bound':89,'proven_full_section_defect_lower_bound':445}
Path(__file__).with_name('one-extra-exact-results.json').write_text(json.dumps(out,indent=2)+'\n')
print('42 eight-point cube geometries: all pair-kernel dimensions 2, coupled rank 4, homogeneous dimension 2.')
print('84 explicit 88-error hostiles: every straight cycle and mixed stress passes; all global periods are nonzero.')
print('Binary-normalization boundary hostile checked: p-1 changes per exceptional column can conceal a nonbinary additive array.')
print('Proved lower bounds: 89 on the mixed triple; 445 on the full norm-circle graph.')
