"""Exact aggregate checks; no optional packages and no removable asserts."""
from fractions import Fraction as F
import json
K=[[172,-24,408],[24,-4,68],[408,-68,706]]
J=[[1,0,0],[0,-1,0],[0,0,1]]
def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def tr(A):return list(map(list,zip(*A)))
def det(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2]
 return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def require(ok,label):
 if not ok:raise RuntimeError(label)
v=[[1],[0],[0]];cols=[];hs=[]
for r in range(44):
 hs.append(v[0][0])
 if r<3:cols.append([x[0] for x in v])
 v=mm(K,v)
B=tr(cols)
require(mm(J,K)==mm(tr(K),J),'adjoint')
require(det(K)==50400 and det(B)==-4406400,'determinants')
qv=[[cols[2][i]-874*cols[1][i]-45144*cols[0][i]] for i in range(3)]
require(qv==[[0],[10800],[0]],'negative polynomial')
require(mm(mm(tr(qv),J),qv)==[[-116640000]],'negative norm')
for r in range(3,len(hs)):
 require(hs[r]==874*hs[r-1]+43344*hs[r-2]+50400*hs[r-3],f'recurrence {r}')
for s in range(40):
 require(det([[hs[s+i+j] for j in range(3)] for i in range(3)])==-19416360960000*50400**s,f'shifted determinant {s}')
def p(x):return x**3-874*x*x-43344*x-50400
def dp(x):return 3*x*x-1748*x-43344
L=F('921.115395614954');U=F('921.115395614956')
require(p(-46)<0<p(-45) and p(-2)>0>p(-1),'negative roots')
require(p(L)<0<p(U),'dominant root')
cl=(L*L-702*L+1800)/dp(U);cu=(U*U-702*U+1800)/dp(L)
print(json.dumps({'status':'PASS','K':K,'J':J,'initial_moments':hs[:6], 'krylov_determinant':det(B),'negative_polynomial':[-45144,-874,1], 'negative_square_moment':-116640000, 'shifted_hankel_constant':-19416360960000,'shift_multiplier':50400, 'checked_shifts':40,'rho_interval':[str(L),str(U)],'asymptotic_coefficient_interval':[str(cl),str(cu)],'rho_over_18_interval':[str(L/18),str(U/18)]},indent=2))
