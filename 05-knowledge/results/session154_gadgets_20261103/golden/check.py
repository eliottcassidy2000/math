def require(x):
 if not x: raise AssertionError("7-adic check failed")
# Integer-only certificate for A^48, orbit, derivatives, and small lifts.
A=((11,18),(4,7))
def mul(P,Q,mod):
 return tuple(tuple(sum(P[i][k]*Q[k][j] for k in range(2))%mod for j in range(2)) for i in range(2))
def power(P,n,mod):
 R=((1,0),(0,1))
 while n:
  if n&1:R=mul(R,P,mod)
  P=mul(P,P,mod);n//=2
 return R
def row(v,P,mod):return tuple(sum(v[i]*P[i][j] for i in range(2))%mod for j in range(2))
require(power(A,48,49)==((43,42),(42,1)))
v=(25%7,34%7); orbit=[]
for k in range(48):
 require(v not in orbit and v!=(0,0))
 orbit.append(v);v=row(v,A,7)
require(v==orbit[0] and len(orbit)==48)
for y in range(7):
 require(sum(Y==y and (Y+X)%7!=0 for Y,X in orbit)==6)
for m in (1,2,3):
 mod=7**m; v=(25%mod,34%mod); counts=[0]*mod
 for k in range(48*7**(m-1)):
  counts[(v[0]+6)%mod]+=1;v=row(v,A,mod)
 require(min(counts)>=6)
 print('m=',m,'period bound=',48*7**(m-1),'min/max fibres=',min(counts),max(counts))
print('PASS')
