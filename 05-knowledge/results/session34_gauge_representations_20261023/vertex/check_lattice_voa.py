#!/usr/bin/env python3
import json
from fractions import Fraction
import sympy as s
import mpmath as mp
N=30
p=[1]+[0]*N
for k in range(1,N+1):
 for j in range(k,N+1):p[j]+=p[j-k]
# Independent Euler pentagonal recurrence.
pe=[1]+[0]*N
for j in range(1,N+1):
 for k in range(1,N+1):
  for g in (k*(3*k-1)//2,k*(3*k+1)//2):
   if g<=j: pe[j]+=(-1)**(k+1)*pe[j-g]
assert p==pe
co=[]
for a in (0,1):
 ts=[0]*(N+1)
 for n in range(-N,N+1):
  j=n*n if a==0 else n*(n+1)
  if j<=N:ts[j]+=1
 co.append([sum(ts[k]*p[j-k] for k in range(j+1)) for j in range(N+1)])
assert co[0][:5]==[1,3,4,7,13]
assert co[1][:5]==[2,2,6,8,14]
I=s.I;S=s.Matrix([[1,1],[1,-1]])/s.sqrt(2);theta=s.diag(1,I)
t=(s.sqrt(6)+s.sqrt(2)-I*(s.sqrt(6)-s.sqrt(2)))/4
T=t*theta
simp=lambda M:M.applyfunc(s.simplify)
assert simp(S*S)==s.eye(2)
assert simp((S*T)**3)==s.eye(2)
assert simp((S*theta)**3)==(1+I)/s.sqrt(2)*s.eye(2)
# Verlinde formula, exact.
fusion=[]
for a in range(2):
 for b in range(2):
  for c in range(2):
   z=s.simplify(sum(S[a,j]*S[b,j]*S[c,j]/S[0,j] for j in range(2)))
   assert z==int(c==(a^b));fusion.append([a,b,c,int(z)])
mp.mp.dps=80
pi=mp.pi
ex=lambda x:mp.exp(2*pi*mp.j*x)
def chars(tau):
 eta=ex(tau/24)*mp.fprod(1-ex(tau*n) for n in range(1,401))
 return [sum(ex(tau*(mp.mpf(n)+a)**2) for n in range(-60,61))/eta for a in (mp.mpf(0),mp.mpf('0.5'))]
num=[]
for re,im in [('0.17','0.83'),('0.31','1.2'),('-0.2','0.7')]:
 tau=mp.mpc(re,im);v=chars(tau);sv=chars(-1/tau);tv=chars(tau+1)
 se=max(abs(sv[0]-(v[0]+v[1])/mp.sqrt(2)),abs(sv[1]-(v[0]-v[1])/mp.sqrt(2)))
 te=max(abs(tv[a]-ex(Fraction(-1,24).numerator/mp.mpf(Fraction(-1,24).denominator)+mp.mpf(a)/4)*v[a]) for a in range(2))
 assert se<mp.mpf('1e-65') and te<mp.mpf('1e-65')
 num.append({'tau':f'{re}+({im})i','S_absolute_residual':mp.nstr(se,8),'T_absolute_residual':mp.nstr(te,8)})
print(json.dumps({'status':'exact algebra and coefficient checks; numerical analytic controls, not certified interval bounds','coefficient_order':N,'vacuum_exponent':'-1/24','nonvacuum_exponent':'5/24','vacuum_coefficients':co[0],'nonvacuum_coefficients':co[1],'partition_methods_agree':True,'S_squared':'I','character_ST_cubed':'I','ribbon_ST_cubed':'exp(pi*i/4) I','Verlinde_fusion':fusion,'modular_numeric_checks':num},indent=2))
