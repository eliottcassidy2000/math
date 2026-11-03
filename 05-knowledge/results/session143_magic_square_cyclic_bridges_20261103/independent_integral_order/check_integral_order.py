#!/usr/bin/env python3
"""Independent exact audit of Z[P+P^2] inside Z[C5]. Standard library only.
Uses cyclic convolution, Leibniz minors, determinantal divisors and adjugate.
"""
from itertools import permutations,combinations,product
from math import gcd
import argparse,json
p=argparse.ArgumentParser();p.add_argument('--mutate',choices=['index','growing-cokernel','cyclotomic-inverse','hadamard-repair']);args=p.parse_args()
checks=0
def need(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,k):return [k*x for x in a]
def mul(a,b):
 c=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[(i+j)%5]+=x*y
 return c
I=[1,0,0,0,0];P=[0,1,0,0,0];A=add(P,mul(P,P));powers=[I]
for i in range(1,6):powers.append(mul(powers[-1],A))
M=[[powers[j][i] for j in range(5)] for i in range(5)]
def det(M):
 n=len(M);ans=0
 for s in permutations(range(n)):
  t=(-1)**sum(s[i]>s[j] for i in range(n) for j in range(i+1,n))
  for i in range(n):t*=M[i][s[i]]
  ans+=t
 return ans
D=det(M);need(D==(-10 if args.mutate=='index' else -11),'power-basis determinant minus eleven')
# Every k-minor is evaluated independently; gcds certify the Smith invariants.
divisors=[1];minor_counts=[]
for k in range(1,6):
 g=0;n=0
 for rs in combinations(range(5),k):
  for cs in combinations(range(5),k):
   g=gcd(g,det([[M[r][c] for c in cs] for r in rs]));n+=1
 divisors.append(abs(g));minor_counts.append(n)
need(divisors==[1,1,1,1,1,11],'determinantal divisors')
smith=[divisors[k+1]//divisors[k] for k in range(5)]
need(smith==[1,1,1,1,11],'Smith invariants')
# Monic degree-five relation shows ALL polynomial values are in this lattice.
need(powers[5]==add(add(scale(powers[2],5),scale(powers[1],5)),scale(I,2)),'monic annihilating relation')
inverse_num=add(add(add(add(scale(powers[4],5),scale(powers[3],-1)),scale(powers[2],-2)),scale(A,-18)),scale(I,-17))
need(inverse_num==scale(P,11),'constructive 11P inverse')
# Check exact adjugate too, rather than infer matrix invertibility by float.
adj=[[(-1)**(i+j)*det([[M[r][c] for c in range(5) if c!=i] for r in range(5) if r!=j]) for j in range(5)] for i in range(5)]
for i,j in product(range(5),repeat=2):need(sum(M[i][k]*adj[k][j] for k in range(5))==D*int(i==j),'adjugate identity')

def cyc(v):
 # Quotient by Phi5: subtract v4*(1+z+z²+z³+z⁴).
 return [v[i]-v[4] for i in range(4)]
lam=powers[1]
cinv=add(add(add(scale(powers[3],-1),scale(powers[2],-2)),scale(lam,-3)),scale(I,-2))
if args.mutate=='cyclotomic-inverse':cinv[0]+=1
need(cyc(cinv)==cyc(P),'integral cyclotomic inverse')
need(cinv==add(P,[-5]*5),'full-ring correction is minus five times J')
# Ring multiplication by z-1 and z+2 on the cyclotomic quotient.
def normO(v):
 cols=[]
 for i in range(4):
  e=[0]*5;e[i]=1;cols.append(cyc(mul(v,e)))
 return det([[cols[j][i] for j in range(4)] for i in range(4)])
zminus1=add(P,scale(I,-1));zplus2=add(P,scale(I,2))
need(mul(zminus1,zplus2)==add(A,scale(I,-2)),'lambda-minus-two factorization')
need(normO(zminus1)==5 and normO(zplus2)==11 and normO(add(A,scale(I,-2)))==55,'gluing norms five eleven fiftyfive')
# Trace discriminants of these explicitly based etale orders.
def transpose(M):return list(map(list,zip(*M)))
def matmul(U,V):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*V)] for row in U]
G5=[[5 if (i+j)%5==0 else 0 for j in range(5)] for i in range(5)]
G4=[[4 if (i+j)%5==0 else -1 for j in range(4)] for i in range(4)]
F4=transpose([cyc(v) for v in powers[:4]])
need(abs(det(F4))==1,'cyclotomic power basis unimodular')
discs=[det(G5),det(matmul(matmul(transpose(M),G5),M)),det(G4)]
need(discs==[3125,378125,125],'trace discriminants')
need(discs[1]==121*discs[0],'discriminant index-square law')
need(det(matmul(matmul(transpose(F4),G4),F4))==125,'unchanged cyclotomic discriminant')
# Full labeled adjacency has more operations than ordinary polynomial calculus.
repair_moduli=[1,2,3,4,5,7,11,22,55,121,1331]
for step in range(1,5):
 Q=[[int((j-i)%5==step) for j in range(5)] for i in range(5)]
 Q2=matmul(Q,Q);B=[[Q[i][j]+Q2[i][j] for j in range(5)] for i in range(5)]
 B2=matmul(B,B);B3=matmul(B2,B)
 J=[[int(i==j)+B[i][j]+B[j][i] for j in range(5)] for i in range(5)]
 need(J==[[1]*5 for _ in range(5)],'converse sidecar recovers all-ones matrix')
 repair=[[5*J[i][j]-B3[i][j]-2*B2[i][j]-3*B[i][j]-2*int(i==j) for j in range(5)] for i in range(5)]
 need(repair==Q,'sidecar repairs every oriented generator integrally')
 sign=1 if args.mutate=='hadamard-repair' else -1
 had=[[B[i][j]+sign*B[i][j]*B2[i][j] for j in range(5)] for i in range(5)]
 need(had==Q,'entrywise Hadamard repair recovers labeled generator')
 for m in repair_moduli:
  need(all((had[i][j]-Q[i][j])%m==0 for i in range(5) for j in range(5)),'repair survives coefficient reduction')
# Evaluation difference is linear and surjective, not a ring homomorphism.
def ev(v,x,m):return sum(c*pow(x,i,m) for i,c in enumerate(v))%m
need(all(ev(v,1,11)==ev(v,9,11) for v in powers),'all A powers satisfy equalizer condition')
need((ev(P,9,11)-ev(P,1,11))%11==8,'evaluation difference surjective on additive group')
need(ev(P,1,11)!=ev(P,9,11),'P not in index-eleven subring')
# Full split algebra at 11: exactly one collision.
roots=[r for r in range(11) if pow(r,5,11)==1]
images={r:(r+r*r)%11 for r in roots}
need(roots==[1,3,4,5,9] and images=={1:2,3:1,4:9,5:8,9:2},'split mod11 roots and collision')
def polymul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
cp=[1]
for v in images.values():cp=polymul(cp,[-v,1]);cp=[x%11 for x in cp]
need(cp==[x%11 for x in [-2,-5,-5,0,0,1]],'mod11 characteristic factorization')
# Smith reduction: additive cokernel is Z/gcd(m,11), never 11^v_11(m).
moduli=[1,2,3,5,7,11,13,17,22,25,55,121,1331,14641]
mods=[]
for m in moduli:
 cokernel=gcd(m,11)
 claimed=(m if args.mutate=='growing-cokernel' and m==121 else cokernel)
 need(claimed==cokernel,'bounded module cokernel does not grow eleven-adically')
 if cokernel==1 and m>1:
  inv11=pow(11,-1,m)
  need([inv11*x%m for x in inverse_num]==[x%m for x in P],'modular polynomial inverse')
 mods.append({'modulus':m,'additive_cokernel_order':cokernel,'full_algebra_generated':cokernel==1})
# Characteristic-five boundary: spectrum coalesces but generator still works.
q=[1]
for _ in range(5):q=polymul(q,[-2,1])
need([x%5 for x in q]==[x%5 for x in [-2,-5,-5,0,0,1]],'mod5 coalesced characteristic polynomial')
need([x%5 for x in inverse_num]==[x%5 for x in P],'mod5 full nonsemisimple generation')
need([x%5 for x in mul(add(P,scale(I,-1)),add(P,scale(I,-1)))]!=[0]*5,'mod5 nonzero nilpotence boundary')
print(json.dumps({'status':'PASS','checks':checks,'power_basis_matrix':M,'determinant':D,'minor_counts':minor_counts,
 'determinantal_divisors':divisors,'smith_invariants':smith,'cyclotomic_inverse_coefficients':[-2,-3,-2,-1],
 'adjacency_repair_family_size':4,'adjacency_repair_moduli':repair_moduli,'trace_discriminants':discs,'gluing_norms':[5,11,55],'roots_mod11':roots,'A_values_mod11':images,'characteristic_mod11_ascending':cp,
 'modular_cokernels':mods,'scope':'Additive order/module defect and exact algebra generation. Cokernel is not a quotient ring. No FRACTRAN state/transition transport or Baby Monster embedding inferred.'},indent=2,sort_keys=True))
