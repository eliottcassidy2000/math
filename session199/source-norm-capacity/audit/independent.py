from pathlib import Path
from fractions import Fraction
from itertools import product
import json,hashlib,math
G=[((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))]
b=(0,0,0,1);u=(1,1,2,2)
P={'N':((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13)),'M':((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29)),'A':((14,10,21,20),(10,7,16,15),(21,16,33,31),(20,15,31,29)),'B':((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),'C':((23,18,44,40),(18,14,34,31),(44,34,80,73),(40,31,73,67)),'J0':((12,12,24,24),(12,12,21,21),(15,12,30,27),(15,12,30,27)),'J1':((12,12,21,21),(12,12,18,18),(18,15,32,30),(18,15,31,29))}
for x,y in [('T','B'),('J2','J1'),('J3','J0')]:P[x]=tuple(zip(*P[y]))
def need(c,s):
 if not c:raise RuntimeError(s)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def rm(r,A):return tuple(dot(r,c) for c in zip(*A))
def mv(A,v):return tuple(dot(r,v) for r in A)
def minimum(S):
 R=[]
 for v in sorted(set(S),key=lambda x:sum(x)):
  if not any(all(a<=b for a,b in zip(w,v)) for w in R):R.append(v)
 return set(R)
def old(v):
 a,f,e,p=v;return a+2*e+p,e,a+f+2*e+2*p,e+p
R=[{u}]
for k in range(40):R.append(minimum(rm(v,A) for v in R[-1] for A in G))
need(R[7]=={tuple(5*x for x in v) for v in R[4]},'full recurrence-state equality')
for k in range(4,38):need(R[k+3]=={tuple(5*x for x in v) for v in R[k]},'state periodicity')
S={b};mins=[]
for k in range(11):
 exact=min(dot(u,v) for v in S);need(exact==min(r[3] for r in R[k]),'unpruned source state minima')
 mins.append(exact);S={mv(A,v) for v in S for A in G}
for k in range(41):
 m=k+2
 expected=2 if m==2 else (3,5,9)[(m-3)%3]*5**((m-3)//3)
 need(min(r[3] for r in R[k])==expected,'all residue-class minima')
 need(expected**3>=5**(m-1),'cover growth inequality')
# Independent physical complete bank, with U increasing under each step.
bank={b};queue=[b]
for v in queue:
 for A in G:
  w=mv(A,v)
  need(dot(u,w)>dot(u,v),'physical U increase')
  if dot(u,w)<=256 and w not in bank:bank.add(w);queue.append(w)
for v in bank:
 D,E,U,d=old(v);need((U+E)%2==0,'physical half-integrality')
 x,y,z=D-E+U,(U+E)//2,U-d
 need(2*y<x<4*y and x%2==1 and y>=1 and z>=1,'physical norm chamber')
 need(dot(v,mv(P['N'],v))==x*x+4*y*z,'diagonal norm')
 X,Y,Z=x-2*y,x-y+z,y
 need(X>0 and Y>0 and Z>0 and X*X+4*Y*Z==x*x+4*y*z and X<2*Y,'norm move hostile')
for name,A in P.items():need(all(A[i][j]>=3*u[i]*u[j] for i in range(4) for j in range(4)),'all product bounds')
# Exact quadratic/bilinear identity on rational matrices.
T=((2,1,3,3),(Fraction(1,2),Fraction(1,2),Fraction(3,2),1),(1,1,1,1))
Q=((1,0,0),(0,0,2),(0,2,0))
for i in range(4):
 for j in range(4):need(sum(T[a][i]*Q[a][c]*T[c][j] for a in range(3) for c in range(3))==P['N'][i][j],'norm bilinear matrix identity')
path=Path(__file__).resolve().parents[1]
report={'status':'PASS','proof_sha256':hashlib.sha256((path/'PROOF.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256((path/'verify.py').read_bytes()).hexdigest(),'full_Pareto_state_scaling_identity':True,'Pareto_recurrence_checked_through':40,'unpruned_source_state_minima_through_word_length10':mins,'physical_U_cap':256,'complete_profiles_in_bank':len(bank),'all_ten_product_matrices_verified':True,'exact_rational_bilinear_norm_identity':True,'manual_proof_audit':'Sharp all-order minima, attainability, rooted-class exponent, norm chamber/involution, common-source density, sublogarithmic thickening, and fixed-bank multiplicative closure are valid. No mathematical correction needed.','scope':'Independent-source spectrum remains open.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
